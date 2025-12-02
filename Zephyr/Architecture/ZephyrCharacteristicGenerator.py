import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj
import bisect

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "RP1"

TARGET_THRUST_LBF = 5000.0
CHAMBER_PRESSURES = [400]
OF_TARGET = 2.2
OF_RANGE = np.linspace(0.5, 3.0, 50)

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
LSTAR = 1.1

# --- Wall / material config  ---
C_bartz = 0.026             # empirical Bartz coefficient used in h estimate
WALL_THICKNESS_M = 0.01     # 1 cm total wall thickness (m)
WALL_NODES = 24            # number of 1D finite-difference nodes through thickness

# GRCop-42 nominal approximate values (temperature-dependent k used)
WALL_MATERIAL = {
    'name': 'GRCop-42',
    'rho': 8790.0,   # kg/m3 (NASA data ~8.79 g/cc). 
    'cp': 385.0      # J/kg-K (approx typical copper-alloy cp)
}

# outer convection coefficient (ambient to outer wall) for transient conduction model
H_EXT_AMBIENT = 10.0  # W/m^2-K, natural convection approx (user can change)

# =====================================================
# === CONSTANTS =======================================
# =====================================================

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822
R_UNIVERSAL = 8314.4621
KG_TO_LB = 2.20462
GAL_PER_LITER = 0.264172
L_PER_M3 = 1000.0
M3_TO_GAL = GAL_PER_LITER * L_PER_M3

DENSITY_RP1 = 810.0
DENSITY_LOX = 1140.0

# unit conversion constants used for CEA outputs
BTU_TO_J = 1055.05585
LBM_TO_KG = 0.45359237
LBM_PER_CUFT_TO_KG_PER_M3 = 16.01846337396
FT_TO_M = 0.3048
MILLIPOISE_TO_PA_S = 1e-4
SIGMA_SB = 5.670374419e-8  # Stefan-Boltzmann (not used but left for extension)

# =====================================================
# === MATERIAL PROPERTY HELPERS (from NASA figures) ==
# =====================================================

# Thermal conductivity vs temperature (digitized from NASA Gradl et al. 2019 Fig.15)
# Temperatures in °C then converted to K
# Source: Gradl et al., "GRCop-42 Development..." (2019) Fig.15 — values approximated from figure.
_GRCOP42_T_C = np.array([0., 100., 200., 300., 400., 500., 600., 700.])
_GRCOP42_K   = np.array([275, 276, 278, 279, 280, 275, 265, 255])  # W/m·K approx
_GRCOP42_T_K = _GRCOP42_T_C + 273.15

def thermal_conductivity_GRCop42(T_K):
    """Interpolate thermal conductivity for GRCop-42 (T in K). Returns k in W/mK.
    Source: Gradl et al. 2019, Fig.15 (values digitized / approximated).
    """
    if T_K <= _GRCOP42_T_K[0]:
        return float(_GRCOP42_K[0])
    if T_K >= _GRCOP42_T_K[-1]:
        return float(_GRCOP42_K[-1])
    idx = bisect.bisect_left(_GRCOP42_T_K, T_K)
    T1, T2 = _GRCOP42_T_K[idx-1], _GRCOP42_T_K[idx]
    k1, k2 = _GRCOP42_K[idx-1], _GRCOP42_K[idx]
    k = k1 + (k2 - k1) * (T_K - T1) / (T2 - T1)
    return float(k)

# Yield-strength helper remains available for hoop check (digitized from Gradl et al. 2019 Fig.6)
_YIELD_T_C = np.array([0., 100., 200., 300., 400., 500., 600., 700., 800.])
_YIELD_GRCOP42_MPA = np.array([200., 189., 172., 158., 140., 125., 100., 70., 30.])  # MPa approx
_YIELD_T_K = _YIELD_T_C + 273.15

def yield_strength_GRCop42(T_K):
    """Return 0.2% offset yield strength (MPa) interpolated from NASA Gradl et al. 2019 Fig.6."""
    if T_K <= _YIELD_T_K[0]:
        return float(_YIELD_GRCOP42_MPA[0])
    if T_K >= _YIELD_T_K[-1]:
        return float(_YIELD_GRCOP42_MPA[-1])
    idx = bisect.bisect_left(_YIELD_T_K, T_K)
    T1, T2 = _YIELD_T_K[idx-1], _YIELD_T_K[idx]
    s1, s2 = _YIELD_GRCOP42_MPA[idx-1], _YIELD_GRCOP42_MPA[idx]
    s = s1 + (s2 - s1) * (T_K - T1) / (T2 - T1)
    return float(s)

# =====================================================
# === HELPER FUNCTIONS ================================
# =====================================================

def cstar_from_T_gamma_R(Tc, gamma, R_specific):
    term = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return math.sqrt(R_specific * Tc) / (gamma ** 0.5 * term)

def pe_over_pc_from_M(M, gamma):
    return (1.0 + 0.5 * (gamma - 1.0) * M**2) ** (-gamma / (gamma - 1.0))

def area_ratio_from_M(M, gamma):
    term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M**2)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / M) * (term ** exponent)

def Mach_from_pe_pc(pe_pc, gamma):
    if pe_pc <= 0 or pe_pc >= 1:
        raise ValueError("pe_pc must be between 0 and 1 (exclusive).")
    term = (pe_pc) ** (-(gamma - 1.0) / gamma)
    M = math.sqrt((2.0 / (gamma - 1.0)) * (term - 1.0))
    return M

def compute_Cf_ideal(gamma, Pe_Pc, Pa_Pc, Ae_over_At):
    term1 = (2.0 * gamma**2) / (gamma - 1.0)
    exponent = (gamma + 1.0) / (gamma - 1.0)
    bracket = (2.0 / (gamma + 1.0)) ** exponent
    momentum = math.sqrt(term1 * bracket * (1.0 - Pe_Pc ** ((gamma - 1.0) / gamma)))
    pressure_term = (Pe_Pc - Pa_Pc) * Ae_over_At
    return momentum + pressure_term

# =====================================================
# === CORE CEA WRAPPER ================================
# =====================================================

def get_numbers_extended(OF_ratio, pc_psi, ox, fuel):
    """
    Wrap CEA_Obj calls and return essential thermochemical and transport properties.
    Minimal additional CEA queries (cp, density, transport, sonic velocity).
    Conversions to SI are applied here.
    """
    cea = CEA_Obj(oxName=ox, fuelName=fuel)

    # temperatures (degR), convert primary chamber temperature to K
    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0  # Rankine -> Kelvin

    # molecular weight and gamma (MolWt returned in lbm/lbmole)
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt, gamma = mw_gamma

    # R_specific
    R_specific = R_UNIVERSAL / MolWt

    # cstar
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)

    # --- sparing additional CEA queries (units converted to SI) ---
    # cp: get_Chamber_Cp returns BTU/lbm-degR -> J/kg-K
    cp_btu_per_lbm_degR = cea.get_Chamber_Cp(Pc=pc_psi, MR=OF_ratio, eps=1)
    cp = cp_btu_per_lbm_degR * BTU_TO_J / LBM_TO_KG * (9.0 / 5.0)

    # density: get_Chamber_Density returns lbm/ft^3 -> kg/m^3
    rho_lbm_per_cuft = cea.get_Chamber_Density(Pc=pc_psi, MR=OF_ratio, eps=1)
    rho = rho_lbm_per_cuft * LBM_PER_CUFT_TO_KG_PER_M3

    # transport: get_Chamber_Transport returns [heatcap, viscosity (mP), thermal_cond (?), Pr]
    transport = cea.get_Chamber_Transport(Pc=pc_psi, MR=OF_ratio, eps=1, frozen=0)
    mu_millipoise = transport[1]
    mu = mu_millipoise * MILLIPOISE_TO_PA_S

    # sonic velocity (ft/s) -> convert to m/s
    sonic_vel_ft_s = cea.get_Chamber_SonicVel(Pc=pc_psi, MR=OF_ratio, eps=1)
    sonic_vel = sonic_vel_ft_s * FT_TO_M

    return {
        'cstar': cstar,
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
        'cp': cp,
        'mu': mu,
        'rho_chamber': rho,
        'sonic_vel': sonic_vel
    }

# =====================================================
# === TRANSIENT 1D WALL MODEL (throat & cylinder) =====
# =====================================================

def bartz_h(Tg, Tw, mu, cp, G, rt):
    """
    Bartz-like empirical convective coefficient with Tg/Tw correction.
    h = C_bartz * mu^0.2 * cp^0.6 * G^0.8 * rt^-0.1 * (Tg/Tw)^0.68
    Units:
      mu [Pa.s], cp [J/kg-K], G [kg/(m^2 s)], rt [m], Tg/Tw dimensionless
    Returns h in W/m^2-K
    """
    # safeguard
    mu = max(mu, 1e-8)
    cp = max(cp, 1e-8)
    G = max(G, 1e-12)
    rt = max(rt, 1e-8)
    Tw_safe = max(Tw, 1.0)
    h0 = C_bartz * (mu ** 0.2) * (cp ** 0.6) * (G ** 0.8) * (rt ** -0.1)
    h = h0 * ((Tg / Tw_safe) ** 0.68)
    return h

def transient_1d_wall(Tg, cp_g, mu, G, rt, wall_props, t_total, nodes=WALL_NODES, h_ext=H_EXT_AMBIENT):
    """
    Solve 1D transient conduction through wall thickness with inner convection (Bartz) and outer convection to ambient.
    - Tg: gas temperature (K)
    - cp_g: gas cp (J/kg-K) [used in bartz_h]
    - mu: gas viscosity (Pa.s)
    - G: mass flux (kg/m2-s) through throat area
    - rt: local radius used by Bartz (m)
    - wall_props: dict with 'rho', 'cp' (material)
    - t_total: total simulation time (s)
    - nodes: number of FD nodes
    - h_ext: outer convective coefficient (W/m2-K)
    Returns:
      Tw_inner_time_series (list), times (list), q_inner_initial (W/m2), h_inner_initial (W/m2-K)
    Notes: explicit scheme used; dt chosen for stability using local thermal diffusivity.
    """
    rho_w = wall_props['rho']
    cp_w = wall_props['cp']
    L = WALL_THICKNESS_M
    N = int(nodes)
    dx = L / N

    # initial temperature field (uniform ambient ~300 K)
    T_ambient = 300.0
    T = np.ones(N) * T_ambient

    # initial Tw for bartz calculation (inner surface temperature)
    Tw = T[0]

    # compute initial h
    h_initial = bartz_h(Tg, Tw, mu, cp_g, G, rt)
    q_initial = h_initial * (Tg - Tw)

    # thermal diffusivity based on initial k
    k_init = thermal_conductivity_GRCop42(Tw)
    alpha = k_init / (rho_w * cp_w)
    # stable explicit dt constraint: dt <= dx^2/(2*alpha)
    if alpha <= 0:
        dt_stable = 1e-3
    else:
        dt_stable = dx * dx / (2.0 * alpha)
    # choose dt conservatively but not too small
    dt = min(0.05, max(1e-4, 0.4 * dt_stable))

    nsteps = max(1, int(t_total / dt))
    times = []
    Tw_inner_series = []

    # loop
    for step in range(nsteps):
        # inner surface temperature used for h calc
        Tw = T[0]
        # recompute h using current Tw
        h_inner = bartz_h(Tg, Tw, mu, cp_g, G, rt)
        # inner heat flux (W/m2) into wall -> q = h*(Tg - Tw)
        q_inner = h_inner * (Tg - Tw)

        # finite difference explicit update (nodes 0..N-1)
        T_new = np.copy(T)
        for i in range(N):
            # local thermal conductivity at node temperature
            k_local = thermal_conductivity_GRCop42(T[i])
            # For internal nodes use central difference with neighbor k averaged
            if i == 0:
                # inner boundary node: convection from gas into node
                # use ghost node formulation: q_inner = -k * (T_ghost - T0)/dx  -> approximate T_ghost
                # Simpler explicit: heat in per unit volume = q_inner / (rho*cp*dx)
                dT_inner = (q_inner) / (rho_w * cp_w * dx)
                # conduction to neighbor:
                k_right = 0.5 * (k_local + thermal_conductivity_GRCop42(T[i+1]))
                conduction = k_right * (T[i+1] - T[i]) / (dx * dx) / (rho_w * cp_w)
                T_new[i] = T[i] + dt * (dT_inner + conduction)
            elif i == N-1:
                # outer boundary node: convection to ambient
                k_left = 0.5 * (k_local + thermal_conductivity_GRCop42(T[i-1]))
                conduction = k_left * (T[i-1] - T[i]) / (dx * dx) / (rho_w * cp_w)
                conv_out = h_ext * (T_ambient - T[i]) / (rho_w * cp_w * dx)
                T_new[i] = T[i] + dt * (conduction + conv_out)
            else:
                # interior node: conduction from both sides
                k_left = 0.5 * (k_local + thermal_conductivity_GRCop42(T[i-1]))
                k_right = 0.5 * (k_local + thermal_conductivity_GRCop42(T[i+1]))
                term_left = k_left * (T[i-1] - T[i])
                term_right = k_right * (T[i+1] - T[i])
                T_new[i] = T[i] + dt * ( (term_left + term_right) / (dx * dx) ) / (rho_w * cp_w)

        T = T_new
        times.append((step + 1) * dt)
        Tw_inner_series.append(T[0])

    # return inner temperature time history, times, initial q and h
    return Tw_inner_series, times, q_initial, h_initial

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
    cp = data.get('cp', None)
    mu = data.get('mu', None)
    rho_chamber = data.get('rho_chamber', None)
    sonic_vel = data.get('sonic_vel', None)

    Pc_Pa = pc_psi * PSI_TO_PA
    Pa_Pa = AMBIENT_P_PSI * PSI_TO_PA
    Pa_over_Pc = Pa_Pa / Pc_Pa

    # === Nozzle parameters ===
    if Pa_over_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_over_Pc, gamma)
        eps = area_ratio_from_M(Me, gamma)
        Pe_Pc = pe_over_pc_from_M(Me, gamma)
    else:
        Me, eps, Pe_Pc = 1.0, 1.0, 1.0

    Cf = compute_Cf_ideal(gamma, Pe_Pc, Pa_over_Pc, eps)
    At = F_newtons / (Pc_Pa * Cf)
    Dt = math.sqrt(4.0 * At / math.pi)

    cstar_eff = cstar * EFFICIENCY_FACTOR
    mdot = Pc_Pa * At / cstar_eff

    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Isp = (Cf * cstar) / G0
    Ve = Isp * G0 * EFFICIENCY_FACTOR

    # === Chamber geometry ===
    D_chamber = 2.0 * Dt
    V_chamber = LSTAR * At
    L_chamber = V_chamber / (math.pi/4 * D_chamber**2)

    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    V_converge = (1/3) * math.pi * L_converge * ((D_chamber/2)**2 + (D_chamber/2)*(Dt/2) + (Dt/2)**2)

    V_total = V_chamber + V_converge
    Lstar_with_converge = V_total / At

    R_throat = 0.382 * Dt  # Throat radius-of-curvature approximation (not geometric throat radius)
    L_throat = R_throat + 0.5 * Dt
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    # Surface areas
    A_cyl_wall = math.pi * D_chamber * L_chamber
    r1 = D_chamber / 2.0
    r2 = Dt / 2.0
    slant_converge = math.sqrt((r1 - r2)**2 + L_converge**2)
    A_converge_wall = math.pi * (r1 + r2) * slant_converge
    A_total_wall = A_cyl_wall + A_converge_wall

    # assemble base result dict (includes transport props returned from CEA)
    result = {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber': L_chamber*100.0,
        'L_chamber_cyl': L_chamber, 'L_converge': L_converge,
        'L_chamber_full': L_chamber + L_converge, 'V_chamber': V_chamber,
        'Lstar': LSTAR, 'V_total': V_total, 'Lstar_with_converge': Lstar_with_converge,
        'L_throat': L_throat*100.0, 'L_nozzle': L_nozzle*100.0,
        'A_cyl_wall': A_cyl_wall * 10000,
        'A_converge_wall': A_converge_wall * 10000,
        'A_total_wall': A_total_wall * 10000,
        'cp': cp, 'mu': mu, 'rho_chamber': rho_chamber,
        'sonic_vel': sonic_vel, 'R_specific': R_specific
    }

    # --- TRANSIENT: compute throat inner-surface transient using 1D conduction + Bartz h ---
    # Mass flux G = mdot / At (kg / (m^2 s))
    G = max(result['mdot'] / result['At'], 1e-8)
    # Use throat curvature radius for Bartz
    rt_throat = max(R_throat, 1e-6)

    Tw_series_throat, times_throat, q0_throat, h0_throat = transient_1d_wall(
        Tg=result['Tc'],
        cp_g=result['cp'],
        mu=result['mu'],
        G=G,
        rt=rt_throat,
        wall_props=WALL_MATERIAL,
        t_total=HOTFIRE_SECONDS,
        nodes=WALL_NODES,
        h_ext=H_EXT_AMBIENT
    )

    # Peak throat inner surface temperature and associated time
    if len(Tw_series_throat) > 0:
        Tw_peak_throat = max(Tw_series_throat)
        t_at_peak_throat = times_throat[Tw_series_throat.index(Tw_peak_throat)]
    else:
        Tw_peak_throat = float('nan')
        t_at_peak_throat = 0.0

    # --- SIMPLE ESTIMATE FOR CYLINDRICAL CHAMBER WALL ---
    # Chamber sees lower convective heating than throat. Approximate inner convection coefficient as fraction of throat.
    # We reuse the Bartz baseline but scale mass-flux effect by larger hydraulic radius. Simpler: use h_chamber = 0.2 * h_throat
    h_chamber_est = 0.2 * h0_throat
    q0_chamber = h_chamber_est * (result['Tc'] - 300.0)  # initial approx heat flux into chamber wall

    # run a transient wall model for chamber inner surface using a fixed h_inner = h_chamber_est
    # we adapt transient_1d_wall by temporarily using a near-constant h_inner via simplified call: emulate by setting G small and rt very large
    Tw_series_chamber, times_chamber, q0_ch, h0_ch = transient_1d_wall(
        Tg=result['Tc'],
        cp_g=result['cp'],
        mu=result['mu'],
        G=max(G*0.05, 1e-6),   # much smaller effective mass flux for chamber region
        rt=max(result['D_chamber']/2.0, 1e-6),
        wall_props=WALL_MATERIAL,
        t_total=HOTFIRE_SECONDS,
        nodes=WALL_NODES,
        h_ext=H_EXT_AMBIENT
    )
    if len(Tw_series_chamber) > 0:
        Tw_peak_chamber = max(Tw_series_chamber)
        t_at_peak_chamber = times_chamber[Tw_series_chamber.index(Tw_peak_chamber)]
    else:
        Tw_peak_chamber = float('nan')
        t_at_peak_chamber = 0.0

    # add transient results to result dict
    result['Tw_series_throat'] = Tw_series_throat
    result['Tw_times_throat'] = times_throat
    result['Tw_peak_throat'] = Tw_peak_throat
    result['Tw_peak_time_throat'] = t_at_peak_throat
    result['q0_throat'] = q0_throat
    result['h0_throat'] = h0_throat

    result['Tw_series_chamber'] = Tw_series_chamber
    result['Tw_times_chamber'] = times_chamber
    result['Tw_peak_chamber'] = Tw_peak_chamber
    result['Tw_peak_time_chamber'] = t_at_peak_chamber
    result['q0_chamber_est'] = q0_chamber
    result['h0_chamber_est'] = h_chamber_est

    # Hoop stresses (thin-wall approximation) computed once and included
    Pc_Pa = pc_psi * PSI_TO_PA
    r_chamber = result['D_chamber'] / 2.0
    t = WALL_THICKNESS_M
    sigma_hoop_chamber = Pc_Pa * r_chamber / t
    result['sigma_hoop_chamber_Pa'] = sigma_hoop_chamber

    r_throat_geom = result['Dt'] / 2.0
    sigma_hoop_throat = Pc_Pa * r_throat_geom / t
    result['sigma_hoop_throat_Pa'] = sigma_hoop_throat

    # yield values for checks (based on peak throat temp)
    result['yield_at_Tw_peak_throat_MPa'] = yield_strength_GRCop42(result['Tw_peak_throat']) if not math.isnan(result['Tw_peak_throat']) else None
    result['yield_at_Tc_MPa'] = yield_strength_GRCop42(result['Tc'])

    return result

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*60)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" gamma = {r['gamma']:.4f} | eps = {r['eps']:.3f}")
    print(f" Chamber Temp (gas) = {r['Tc']:.1f} K")
    print(f" Isp = {r['Isp']:.2f} s | Ve = {r['Ve']:.1f} m/s | c* = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Chamber L = {r['L_chamber']:.2f} cm")
    print(f" Converging length = {r['L_converge']*100:.2f} cm")
    print(f" Full chamber length= {r['L_chamber_full']*100:.2f} cm")
    print(f" Chamber volume (with converging section) = {r['V_chamber']*1e6:.2f} cm³")
    print(f" L* (just cylinder) = {r['Lstar']:.3f} m | L* (with converge) = {r['Lstar_with_converge']:.3f} m")
    print(f" Throat length = {r['L_throat']:.2f} cm | Nozzle length = {r['L_nozzle']:.2f} cm")
    print(f" SA Cyl = {r['A_cyl_wall']:.2f} cm^2 | SA Cyl+converge = {r['A_total_wall']:.2f} cm^2")

    # Throat transient results
    Tw_peak_throat = r.get('Tw_peak_throat', None)
    t_peak_throat = r.get('Tw_peak_time_throat', None)
    q0_th = r.get('q0_throat', None)
    h0_th = r.get('h0_throat', None)
    if Tw_peak_throat is not None:
        print(f"\nPeak inner wall temperature (throat) = {Tw_peak_throat:.1f} K at t = {t_peak_throat:.2f} s")
        if q0_th is not None:
            print(f"Initial throat heat flux (approx) = {q0_th:.1f} W/m^2")
        if h0_th is not None:
            print(f"Initial throat convective coeff (approx) = {h0_th:.1f} W/m^2K")

    # Chamber inner surface estimate
    Tw_peak_ch = r.get('Tw_peak_chamber', None)
    t_peak_ch = r.get('Tw_peak_time_chamber', None)
    q0_ch = r.get('q0_chamber_est', None)
    h0_ch = r.get('h0_chamber_est', None)
    if Tw_peak_ch is not None:
        print(f"\nPeak inner wall temperature (cyl chamber) = {Tw_peak_ch:.1f} K at t = {t_peak_ch:.2f} s")
        if q0_ch is not None:
            print(f"Initial chamber heat flux (est) = {q0_ch:.1f} W/m^2")
        if h0_ch is not None:
            print(f"Initial chamber convective coeff (est) = {h0_ch:.1f} W/m^2K")

    # Hoop stress & yield checks
    sigma_ch_Pa = r.get('sigma_hoop_chamber_Pa', None)
    sigma_th_Pa = r.get('sigma_hoop_throat_Pa', None)
    yield_th_peak = r.get('yield_at_Tw_peak_throat_MPa', None)
    if sigma_ch_Pa is not None:
        sigma_ch_MPa = sigma_ch_Pa / 1e6
        print(f"\nApprox chamber hoop stress (thin-wall) = {sigma_ch_MPa:.2f} MPa (at Pc = {r['pc_psi']} psi)")
    if sigma_th_Pa is not None:
        sigma_th_MPa = sigma_th_Pa / 1e6
        print(f"Approx throat hoop stress (thin-wall)   = {sigma_th_MPa:.2f} MPa")

    if yield_th_peak is not None:
        print(f"0.2% yield strength at peak throat wall T = {yield_th_peak:.1f} MPa")
        if sigma_th_MPa >= yield_th_peak:
            print(" **WARNING**: Throat hoop stress >= yield strength at predicted peak wall temperature (plastic risk).")
        elif sigma_th_MPa >= 0.8 * yield_th_peak:
            print(" **CAUTION**: Throat hoop stress >= 80% of yield strength at predicted peak wall temperature.")
        else:
            print("Throat hoop stress is below yield at predicted peak wall temperature.")

    of, mdot, t_fire = r['of'], r['mdot'], HOTFIRE_SECONDS
    mdot_fuel = mdot / (1.0 + of)
    mdot_ox = mdot - mdot_fuel
    m_fuel_total = mdot_fuel * t_fire
    m_ox_total = mdot_ox * t_fire

    V_fuel_gal = (m_fuel_total / DENSITY_RP1) * M3_TO_GAL
    V_ox_gal = (m_ox_total / DENSITY_LOX) * M3_TO_GAL

    print("\n")
    print(f"--- Hotfire {HOTFIRE_SECONDS:.1f}s ---")
    print(f" Fuel: {m_fuel_total:.2f} kg ({V_fuel_gal:.2f} gal)")
    print(f" Oxidizer: {m_ox_total:.2f} kg ({V_ox_gal:.2f} gal)")
    print(f" Total: {(m_fuel_total+m_ox_total):.2f} kg ({(V_fuel_gal+V_ox_gal):.2f} gal)")
    print("="*60 + "\n")

# =====================================================
# === MAIN =============================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (auto-optimal chamber geometry) ===\n")
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N
    show_opt = input("Also show optimal O/F results? (y/n): ").strip().lower()

    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print(r)

        Ve_list, mdot_list, cstar_list, gamma_list, eps_list = [], [], [], [], []
        At_list, Ae_list, Dt_list, De_list, Tc_list, Isp_list = [], [], [], [], [], []

        for of in OF_RANGE:
            r2 = compute_sizing(target_thrust_N, pc, of)
            Ve_list.append(r2['Ve'])
            mdot_list.append(r2['mdot'])
            cstar_list.append(r2['cstar'])
            gamma_list.append(r2['gamma'])
            eps_list.append(r2['eps'])
            At_list.append(r2['At'])
            Ae_list.append(r2['Ae'])
            Dt_list.append(r2['Dt']*100)
            De_list.append(r2['De']*100)
            Tc_list.append(r2['Tc'])
            Isp_list.append(r2['Isp'])

        Isp_arr = np.array(Isp_list)
        idx_isp = np.nanargmax(Isp_arr)
        opt_Isp_OF = OF_RANGE[idx_isp]

        if show_opt == "y":
            print("\n=== OPTIMAL RESULTS SUMMARY ===")
            print(f"At Pc = {pc} psi | Optimal O/F = {opt_Isp_OF:.3f}")
            r_opt_isp = compute_sizing(target_thrust_N, pc, opt_Isp_OF)
            pretty_print(r_opt_isp)

        # ====== Performance Plots ======
        plt.figure(figsize=(10,8))
        plt.suptitle(f"Performance at Pc={pc} psi")

        plt.subplot(2,2,1)
        plt.plot(OF_RANGE, Ve_list, label="Ve (m/s)")
        plt.title("Exhaust Velocity"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2,2,2)
        plt.plot(OF_RANGE, mdot_list, label="mdot (kg/s)")
        plt.title("Mass Flow Rate"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2,2,3)
        plt.plot(OF_RANGE, cstar_list, label="c* (m/s)")
        plt.title("Characteristic Velocity"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2,2,4)
        plt.plot(OF_RANGE, gamma_list, label="Gamma (γ)")
        plt.title("Gamma vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.tight_layout(); plt.show()

        # ====== Geometry & Temp Plots ======
        plt.figure(figsize=(10,8))
        plt.suptitle(f"Geometry & Temperature at Pc={pc} psi")

        plt.subplot(2,2,1)
        plt.plot(OF_RANGE, np.array(At_list)*1e4, label="At (cm²)")
        plt.plot(OF_RANGE, np.array(Ae_list)*1e4, '--', label="Ae (cm²)")
        plt.title("Areas vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2,2,2)
        plt.plot(OF_RANGE, Dt_list, label="Dt (cm)")
        plt.plot(OF_RANGE, De_list, '--', label="De (cm)")
        plt.title("Diameters vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2,2,3)
        plt.plot(OF_RANGE, eps_list, label="Expansion Ratio (ε)")
        plt.title("Expansion Ratio vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        ax1 = plt.subplot(2,2,4)
        ax1.plot(OF_RANGE, Tc_list, color='red', label="Tc (K)")
        ax1.set_xlabel("O/F"); ax1.set_ylabel("Tc (K)", color='red')
        ax1.tick_params(axis='y', labelcolor='red'); ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(OF_RANGE, Isp_list, color='blue', label="Isp (s)")
        ax2.set_ylabel("Isp (s)", color='blue')
        ax2.tick_params(axis='y', labelcolor='blue')

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

        plt.tight_layout(); plt.show()
