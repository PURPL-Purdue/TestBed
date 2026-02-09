import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj
import bisect

from matplotlib.widgets import Slider


# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "RP1"

TARGET_THRUST_LBF = 500.0
CHAMBER_PRESSURES = [320.0]
OF_TARGET = 2.0
OF_RANGE = np.linspace(0.5, 3.0, 50)

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
LSTAR = 1.0922


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

    return {
        'cstar': cstar,
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
    }

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


    return result

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*60)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    #print(f" gamma = {r['gamma']:.4f} | eps = {r['eps']:.3f}")
    print(f" Chamber Temp (gas) = {r['Tc']:.1f} K")
    #print(f" Isp = {r['Isp']:.2f} s | Ve = {r['Ve']:.1f} m/s | c* = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Chamber L = {r['L_chamber']:.2f} cm")
    #print(f" Converging length = {r['L_converge']*100:.2f} cm")
    print(f" Full chamber length= {r['L_chamber_full']*100:.2f} cm")
    print(f" Chamber volume (with converging section) = {r['V_chamber']*1e6:.2f} cm³")
    #print(f" L* (just cylinder) = {r['Lstar']:.3f} m | L* (with converge) = {r['Lstar_with_converge']:.3f} m")
    print(f" Throat length = {r['L_throat']:.2f} cm | Nozzle length = {r['L_nozzle']:.2f} cm")
    #print(f" SA Cyl = {r['A_cyl_wall']:.2f} cm^2 | SA Cyl+converge = {r['A_total_wall']:.2f} cm^2")

    # of, mdot, t_fire = r['of'], r['mdot'], HOTFIRE_SECONDS
    # mdot_fuel = mdot / (1.0 + of)
    # mdot_ox = mdot - mdot_fuel
    # m_fuel_total = mdot_fuel * t_fire
    # m_ox_total = mdot_ox * t_fire

    # V_fuel_gal = (m_fuel_total / DENSITY_RP1) * M3_TO_GAL
    # V_ox_gal = (m_ox_total / DENSITY_LOX) * M3_TO_GAL

    # print("\n")
    # print(f"--- Hotfire {HOTFIRE_SECONDS:.1f}s ---")
    # print(f" Fuel: {m_fuel_total:.2f} kg ({V_fuel_gal:.2f} gal)")
    # print(f" Oxidizer: {m_ox_total:.2f} kg ({V_ox_gal:.2f} gal)")
    # print(f" Total: {(m_fuel_total+m_ox_total):.2f} kg ({(V_fuel_gal+V_ox_gal):.2f} gal)")
    # print("="*60 + "\n")

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


# =====================================================
# === INTERACTIVE SLIDER UI (DYNAMIC MAX VALUES) =======
# =====================================================

def interactive_engine_dashboard():
    # ---- Slider bounds ----
    THRUST_MIN, THRUST_MAX = 1.0, 1000.0
    OF_MIN, OF_MAX = 0.1, 4.0
    PC_MIN, PC_MAX = 25.0, 500.0

    init_thrust_lbf = TARGET_THRUST_LBF
    init_of = OF_TARGET
    init_pc = CHAMBER_PRESSURES[0]

    # ---- Compute extrema over slider domain ----
    def compute_slider_maxima():
        temps, de, dt, vol = [], [], [], []

        for thrust in [THRUST_MIN, THRUST_MAX]:
            for of in [OF_MIN, OF_MAX]:
                for pc in [PC_MIN, PC_MAX]:
                    r = compute_sizing(thrust * LBF_TO_N, pc, of)
                    temps.append(r['Tc'])
                    de.append(r['De'] * 100)       # <- Exit diameter
                    dt.append(r['Dt'] * 100)
                    vol.append(r['V_chamber'] * 1e6)

        return {
            "Tc": max(temps),
            "De": max(de),
            "Dt": max(dt),
            "V": max(vol),
        }

    max_vals = compute_slider_maxima()

    # ---- Initial state ----
    r0 = compute_sizing(init_thrust_lbf * LBF_TO_N, init_pc, init_of)

    # ---- Figure setup ----
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    plt.subplots_adjust(bottom=0.32)
    fig.suptitle("Rocket Engine Interactive Sizing")

    # ---- Chamber Temperature ----
    temp_line, = axs[0,0].plot([0], [r0['Tc']], marker='o')
    axs[0,0].set_title("Chamber Temperature")
    axs[0,0].set_ylabel("K")
    axs[0,0].set_ylim(0, max_vals["Tc"] * 1.05)
    axs[0,0].set_xticks([])
    axs[0,0].grid(True)

    temp_text = axs[0,0].text(
        0.02, 0.95,
        "",
        transform=axs[0,0].transAxes,
        verticalalignment='top'
    )

    # ---- Exit Diameter ----
    de_line, = axs[0,1].plot([0], [r0['De'] * 100], marker='o')   # <-- changed
    axs[0,1].set_title("Exit Diameter")
    axs[0,1].set_ylabel("cm")
    axs[0,1].set_ylim(0, max_vals["De"] * 1.05)
    axs[0,1].set_xticks([])
    axs[0,1].grid(True)

    de_text = axs[0,1].text(
        0.02, 0.95,
        "",
        transform=axs[0,1].transAxes,
        verticalalignment='top'
    )

    # ---- Throat Diameter ----
    dt_line, = axs[1,0].plot([0], [r0['Dt'] * 100], marker='o')
    axs[1,0].set_title("Throat Diameter")
    axs[1,0].set_ylabel("cm")
    axs[1,0].set_ylim(0, max_vals["Dt"] * 1.05)
    axs[1,0].set_xticks([])
    axs[1,0].grid(True)

    dt_text = axs[1,0].text(
        0.02, 0.95,
        "",
        transform=axs[1,0].transAxes,
        verticalalignment='top'
    )

    # ---- Chamber Volume ----
    vol_line, = axs[1,1].plot([0], [r0['V_chamber'] * 1e6], marker='o')
    axs[1,1].set_title("Chamber Volume")
    axs[1,1].set_ylabel("cm³")
    axs[1,1].set_ylim(0, max_vals["V"] * 1.05)
    axs[1,1].set_xticks([])
    axs[1,1].grid(True)

    vol_text = axs[1,1].text(
        0.02, 0.95,
        "",
        transform=axs[1,1].transAxes,
        verticalalignment='top'
    )

    # ---- Sliders ----
    ax_thrust = plt.axes([0.15, 0.20, 0.65, 0.03])
    ax_of     = plt.axes([0.15, 0.15, 0.65, 0.03])
    ax_pc     = plt.axes([0.15, 0.10, 0.65, 0.03])

    s_thrust = Slider(ax_thrust, "Thrust (lbf)", THRUST_MIN, THRUST_MAX, valinit=init_thrust_lbf)
    s_of     = Slider(ax_of, "O/F", OF_MIN, OF_MAX, valinit=init_of)
    s_pc     = Slider(ax_pc, "Pc (psi)", PC_MIN, PC_MAX, valinit=init_pc)

    # ---- Update ----
    def update(val):
        r = compute_sizing(s_thrust.val * LBF_TO_N, s_pc.val, s_of.val)

        Tc = r['Tc']
        De = r['De'] * 100         # <-- exit diameter
        Dt = r['Dt'] * 100
        Vc = r['V_chamber'] * 1e6

        temp_line.set_ydata([Tc])
        de_line.set_ydata([De])
        dt_line.set_ydata([Dt])
        vol_line.set_ydata([Vc])

        temp_text.set_text(f"Tc = {Tc:.1f} K\nMax = {max_vals['Tc']:.1f} K")
        de_text.set_text(f"De = {De:.2f} cm\nMax = {max_vals['De']:.2f} cm")
        dt_text.set_text(f"Dt = {Dt:.2f} cm\nMax = {max_vals['Dt']:.2f} cm")
        vol_text.set_text(f"V = {Vc:.1f} cm³\nMax = {max_vals['V']:.1f} cm³")

        fig.canvas.draw_idle()

    s_thrust.on_changed(update)
    s_of.on_changed(update)
    s_pc.on_changed(update)

    update(None)
    plt.show()


if __name__ == "__main__":
    interactive_engine_dashboard()
