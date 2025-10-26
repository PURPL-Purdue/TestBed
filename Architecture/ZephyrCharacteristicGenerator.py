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
EFFICIENCY_FACTOR = 1
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
THROAT_LENGTH_FACTOR = 0.5

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
    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt, gamma = mw_gamma
    R_specific = R_UNIVERSAL / MolWt
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)
    return {'cstar': cstar, 'Tc': Tc, 'gamma': gamma, 'MolWt': MolWt, 'R_specific': R_specific}

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
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
    mdot = Pc_Pa * At / cstar
    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Isp = (Cf * cstar) / G0
    Ve = Isp * G0 * EFFICIENCY_FACTOR

    # === Chamber geometry iterative L* ===
    D_chamber = 2.5 * Dt
    L_chamber = 3.0 * Dt
    for _ in range(20):
        V_chamber = (math.pi / 4.0) * D_chamber**2 * L_chamber
        Lstar = V_chamber / At
        if Lstar < 1.0:
            L_chamber *= 1.1
            D_chamber *= 1.05
        elif Lstar > 1.5:
            L_chamber *= 0.9
            D_chamber *= 0.95
        else:
            break

    V_chamber = (math.pi / 4.0) * D_chamber**2 * L_chamber
    Lstar_cylinder = V_chamber / At

    # Converging section
    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    V_converge = (1/3) * math.pi * L_converge * ((D_chamber/2)**2 + (D_chamber/2)*(Dt/2) + (Dt/2)**2)
    V_total = V_chamber + V_converge
    Lstar_with_converge = V_total / At

    # Other geometry
    R_throat = 0.382 * Dt
    L_throat = R_throat + 0.5 * Dt
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    return {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber': L_chamber * 100.0,
        'L_chamber_cyl': L_chamber, 'L_converge': L_converge,
        'L_chamber_full': L_chamber + L_converge, 'V_chamber': V_chamber,
        'Lstar': Lstar_cylinder, 'V_total': V_total, 'Lstar_with_converge': Lstar_with_converge,
        'L_throat': L_throat * 100.0, 'L_nozzle': L_nozzle * 100.0
    }

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
