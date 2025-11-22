import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "RP1"

TARGET_THRUST_LBF = 5000.0
CHAMBER_PRESSURES = [600]
OF_TARGET = 1.5
OF_RANGE = np.linspace(0.5, 3.0, 50)

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
LSTAR = 1.1

# --- Wall / material config  ---
C_bartz = 0.026
WALL_THICKNESS_M = 0.01  # 1 cm
# GRCop-42 properties (approx): density kg/m3, cp J/kg-K, k W/m-K
WALL_MATERIAL = {
    'name': 'GRCop-42',
    'rho': 8900.0,
    'cp': 385.0,
    'k': 344.0
}

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
    Minimal additional CEA calls are made here (cp, density, transport, sonic velocity).
    All conversions to SI are done here so other code uses consistent units.
    """
    cea = CEA_Obj(oxName=ox, fuelName=fuel)

    # temperatures (degR), convert primary chamber temperature to K
    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0  # Rankine -> Kelvin

    # molecular weight and gamma (MolWt returned in lbm/lbmole)
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt, gamma = mw_gamma

    # R_specific: we keep the same approach as before (R_universal / MolWt)
    # (This follows the original code style; rocketcea uses English units internally).
    R_specific = R_UNIVERSAL / MolWt

    # cstar (as before)
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)

    # --- minimal, sparing additional CEA queries (units converted to SI) ---
    # cp: get_Chamber_Cp returns BTU/lbm-degR (default). Convert -> J/kg-K
    cp_btu_per_lbm_degR = cea.get_Chamber_Cp(Pc=pc_psi, MR=OF_ratio, eps=1)
    # Convert: BTU/lbm-degR -> J/kg-K
    # 1 BTU = 1055.05585 J; 1 lbm = 0.45359237 kg; degR->K factor: 1 K = 9/5 degR so multiply by 9/5
    cp = cp_btu_per_lbm_degR * BTU_TO_J / LBM_TO_KG * (9.0 / 5.0)

    # density: get_Chamber_Density returns lbm/ft^3 -> convert to kg/m^3
    rho_lbm_per_cuft = cea.get_Chamber_Density(Pc=pc_psi, MR=OF_ratio, eps=1)
    rho = rho_lbm_per_cuft * LBM_PER_CUFT_TO_KG_PER_M3

    # transport: get_Chamber_Transport returns [heatcap, viscosity, thermal_cond, Pr] in default units
    transport = cea.get_Chamber_Transport(Pc=pc_psi, MR=OF_ratio, eps=1, frozen=0)
    # transport[1] is viscosity in millipoise -> Pa.s (1 mP = 1e-4 Pa.s)
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
# === TRANSIENT WALL TEMPERATURE MODEL (new) ==========
# =====================================================

def compute_wall_temperatures_from_sizing(r):
    """
    Given the sizing dictionary 'r' (returned by compute_sizing),
    compute a Bartz-like convective coefficient at the throat and
    integrate a lumped-capacitance wall temperature rise over HOTFIRE_SECONDS.

    Returns:
      Tw_max (K), t_at_Tw_max (s), q_dot_initial (W/m^2)
    Uses only values already present in r (no new CEA calls).
    """

    # extract required values from r
    Tc = r['Tc']                # gas temperature (K)
    cp = r.get('cp', None)      # J/kg-K
    mu = r.get('mu', None)      # Pa.s
    rho = r.get('rho_chamber', None)  # kg/m3
    At = r['At']                # m2
    Dt = r['Dt']                # m
    mdot = r['mdot']            # kg/s
    sonic_vel = r.get('sonic_vel', None)  # m/s
    gamma = r['gamma']

    # Safety: if any required property missing, return NaNs
    if cp is None or mu is None or rho is None or At is None or Dt is None or mdot is None:
        return float('nan'), float('nan'), float('nan')

    # throat radius (m)
    rt = Dt / 2.0

    # characteristic mass flux (G = mdot / At) [kg/(s m^2)]
    G = mdot / At

    # Approximate throat velocity: use sonic_vel if present (choked flow approx)
    if sonic_vel is not None and sonic_vel > 0:
        Vt = sonic_vel
    else:
        # fallback: ideal sonic estimate (may be inconsistent if R_specific units mismatch)
        try:
            Vt = math.sqrt(gamma * r['R_specific'] * Tc)
        except Exception:
            Vt = max(1.0, mdot / (rho * At))

    # --- Bartz-like heat transfer coefficient (empirical form) ---
    # h = C * mu^0.2 * cp^0.6 * G^0.8 * rt^-0.1
    mu = max(mu, 1e-8)
    cp = max(cp, 1e-8)
    G = max(G, 1e-8)
    rt = max(rt, 1e-8)

    h = C_bartz * (mu ** 0.2) * (cp ** 0.6) * (G ** 0.8) * (rt ** -0.1)

    # initial wall temperature (assume ambient ~300 K)
    Tw = 300.0

    # initial heat flux [W/m^2]
    q0 = h * (Tc - Tw)

    # wall material properties from user config
    rho_w = WALL_MATERIAL['rho']   # kg/m3
    cp_w = WALL_MATERIAL['cp']     # J/kg-K
    t_w = WALL_THICKNESS_M         # m

    # Lumped-capacitance check: compute Biot number quickly (use convection h and k)
    k_w = WALL_MATERIAL['k']
    Bi = h * t_w / k_w if k_w > 0 else 1.0

    # We'll still integrate transient using lumped mass (reasonable if Bi < ~0.1).
    # Use small time step for stability
    dt = 0.1
    nsteps = max(1, int(HOTFIRE_SECONDS / dt))
    times = []
    temps = []

    for i in range(nsteps):
        # convective heat flux at current wall temp
        q = h * (Tc - Tw)  # W/m^2

        # rate of temperature rise (simple lumped capacitance)
        dTdt = q / (rho_w * cp_w * t_w)
        Tw = Tw + dTdt * dt

        times.append((i+1) * dt)
        temps.append(Tw)

    # find peak wall temp and time
    if len(temps) > 0:
        Tw_max = max(temps)
        t_at_Tw_max = times[temps.index(Tw_max)]
    else:
        Tw_max = Tw
        t_at_Tw_max = 0.0

    return Tw_max, t_at_Tw_max, q0

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

    R_throat = 0.382 * Dt  # Throat radius of curvature (not the geometric throat radius)
    L_throat = R_throat + 0.5 * Dt
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    # Surface areas
    A_cyl_wall = math.pi * D_chamber * L_chamber

    # Converging section is a frustum: lateral area = pi * (r1 + r2) * s
    # where s = slant length = sqrt((r1 - r2)^2 + axial_length^2)
    r1 = D_chamber / 2.0
    r2 = Dt / 2.0
    slant_converge = math.sqrt((r1 - r2)**2 + L_converge**2)
    A_converge_wall = math.pi * (r1 + r2) * slant_converge

    A_total_wall = A_cyl_wall + A_converge_wall

    # --- new: assemble base result dict and include the CEA-derived transport properties ---
    result = {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber': L_chamber*100.0,
        'L_chamber_cyl': L_chamber, 'L_converge': L_converge,
        'L_chamber_full': L_chamber + L_converge, 'V_chamber': V_chamber,
        'Lstar': LSTAR, 'V_total': V_total, 'Lstar_with_converge': Lstar_with_converge,
        'L_throat': L_throat*100.0, 'L_nozzle': L_nozzle*100.0,

        # surface areas (converted to cm^2 to match prior code's choice)
        'A_cyl_wall': A_cyl_wall * 10000,  # m² -> cm²
        'A_converge_wall': A_converge_wall * 10000,
        'A_total_wall': A_total_wall * 10000,

        # pass-through of additional properties (SI units)
        'cp': cp,
        'mu': mu,
        'rho_chamber': rho_chamber,
        'sonic_vel': sonic_vel,
        'R_specific': R_specific
    }

    # --- compute wall temperature transient and return peaks (use only values in 'result') ---
    Tw_peak, Tw_peak_time, q0 = compute_wall_temperatures_from_sizing(result)

    # add these outputs to result
    result['Tw_peak'] = Tw_peak
    result['Tw_peak_time'] = Tw_peak_time
    result['q_dot_initial'] = q0

    return result

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*60)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" gamma = {r['gamma']:.4f} | eps = {r['eps']:.3f}")
    print(f" Chamber Temp = {r['Tc']:.1f} K")
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

    # NEW: print peak wall (throat) temperature and time
    Tw_peak = r.get('Tw_peak', None)
    Tw_peak_time = r.get('Tw_peak_time', None)
    q0 = r.get('q_dot_initial', None)
    if Tw_peak is not None:
        print(f" Peak wall temperature (throat) = {Tw_peak:.1f} K at t = {Tw_peak_time:.2f} s")
        if q0 is not None:
            print(f" Initial throat heat flux (approx) = {q0:.1f} W/m^2")

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
    print(f" Total: {(m_fuel_total+m_ox_total):.2f} kg ({V_fuel_gal+V_ox_gal:.2f} gal)")
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
            r = compute_sizing(target_thrust_N, pc, of)
            Ve_list.append(r['Ve']); mdot_list.append(r['mdot'])
            cstar_list.append(r['cstar']); gamma_list.append(r['gamma'])
            eps_list.append(r['eps']); At_list.append(r['At'])
            Ae_list.append(r['Ae']); Dt_list.append(r['Dt']*100)
            De_list.append(r['De']*100); Tc_list.append(r['Tc'])
            Isp_list.append(r['Isp'])

        Isp_arr = np.array(Isp_list)
        idx_isp_opt = np.nanargmax(Isp_arr)
        opt_Isp_OF = OF_RANGE[idx_isp_opt]

        if show_opt == "y":
            print("\n=== OPTIMAL RESULTS SUMMARY ===")
            print(f"At Pc = {pc} psi | Optimal O/F = {opt_Isp_OF:.3f}")
            r_opt_isp = compute_sizing(target_thrust_N, pc, opt_Isp_OF)
            pretty_print(r_opt_isp)

        # ====== Performance Plots ======
        plt.figure(figsize=(10, 8))
        plt.suptitle(f"Performance at Pc={pc} psi")

        plt.subplot(2, 2, 1)
        plt.plot(OF_RANGE, Ve_list, label="Ve (m/s)")
        plt.title("Exhaust Velocity"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 2)
        plt.plot(OF_RANGE, mdot_list, label="mdot (kg/s)")
        plt.title("Mass Flow Rate"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 3)
        plt.plot(OF_RANGE, cstar_list, label="c* (m/s)")
        plt.title("Characteristic Velocity"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 4)
        plt.plot(OF_RANGE, gamma_list, label="Gamma (γ)")
        plt.title("Gamma vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.tight_layout(); plt.show()

        # ====== Geometry & Temp Plots ======
        plt.figure(figsize=(10, 8))
        plt.suptitle(f"Geometry & Temperature at Pc={pc} psi")

        plt.subplot(2, 2, 1)
        plt.plot(OF_RANGE, np.array(At_list)*1e4, label="At (cm²)")
        plt.plot(OF_RANGE, np.array(Ae_list)*1e4, '--', label="Ae (cm²)")
        plt.title("Areas vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 2)
        plt.plot(OF_RANGE, Dt_list, label="Dt (cm)")
        plt.plot(OF_RANGE, De_list, '--', label="De (cm)")
        plt.title("Diameters vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 3)
        plt.plot(OF_RANGE, eps_list, label="Expansion Ratio (ε)")
        plt.title("Expansion Ratio vs O/F"); plt.xlabel("O/F"); plt.grid(True); plt.legend()

        ax1 = plt.subplot(2, 2, 4)
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
