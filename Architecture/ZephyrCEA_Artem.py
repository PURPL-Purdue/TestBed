import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

# Propellants
OXIDIZER = "LOX"
FUEL = "RP1"

# Operating conditions
TARGET_THRUST_LBF = 5000.0       # Thrust [lbf]
CHAMBER_PRESSURES = [600]        # List of chamber pressures [psia]
OF_TARGET = 1.5                  # Single O/F point for diagnostics
OF_RANGE = np.linspace(0.5, 3.0, 50)  # Range for sweep

# Environment and efficiency
AMBIENT_P_PSI = 14.7             # Ambient pressure [psi]
EFFICIENCY_FACTOR = 1            # Efficiency multiplier for Ve (user scaling)

# Geometry and design heuristics
NOZZLE_HALF_ANGLE_DEG = 15.0     # Conical nozzle half-angle
THROAT_LENGTH_FACTOR = 0.5       # Throat length ≈ 0.5 × Dt

# =====================================================
# === CONSTANTS =======================================
# =====================================================

G0 = 9.80665                     # m/s²
PSI_TO_PA = 6894.76              # Pa per psi
LBF_TO_N = 4.44822               # N per lbf
R_UNIVERSAL = 8314.4621          # J / (kmol K)

# =====================================================
# === HELPER FUNCTIONS ================================
# =====================================================

def cstar_from_T_gamma_R(Tc, gamma, R_specific):
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    term = ((gamma + 1.0) / 2.0) ** exponent
    return math.sqrt(R_specific * Tc / gamma) * term

def pe_over_pc_from_M(M, gamma):
    return (1.0 + 0.5 * (gamma - 1.0) * M**2) ** (-gamma / (gamma - 1.0))

def area_ratio_from_M(M, gamma):
    term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M**2)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / M) * (term ** exponent)

def Mach_from_pe_pc(pe_pc, gamma):
    M =  math.sqrt( (2/(gamma-1)) * ((1/pe_pc) ** ((gamma -1)/gamma))-1)
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
    Tc = temps_rankine[0] * 5.0/9.0
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt = mw_gamma[0]
    gamma = mw_gamma[1]
    R_specific = R_UNIVERSAL / MolWt
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)
    return {'cstar': cstar, 'Tc': Tc, 'gamma': gamma, 'MolWt': MolWt, 'R_specific': R_specific}

# =====================================================
# === GEOMETRY & PERFORMANCE COMPUTATION ==============
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
    Pc_Pa = pc_psi * PSI_TO_PA
    Pa_Pa = AMBIENT_P_PSI * PSI_TO_PA
    Pa_over_Pc = Pa_Pa / Pc_Pa

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

    # =====================================================
    # === CHAMBER CALCULATIONS =============================
    # =====================================================

    # Iterate on D_chamber until L* converges between 1.0–1.5 m
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

    # Converging section approximation
    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    V_converge = (1/3) * math.pi * L_converge * ((D_chamber/2)**2 + (D_chamber/2)*(Dt/2) + (Dt/2)**2)
    V_total = V_chamber + V_converge
    Lstar_with_converge = V_total / At

    L_throat = THROAT_LENGTH_FACTOR * Dt
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    return {
        'pc_psi': pc_psi,
        'of': of,
        'gamma': gamma,
        'Tc': Tc,
        'cstar': cstar,
        'Isp': Isp,
        'Ve': Ve,
        'Cf': Cf,
        'eps': eps,
        'mdot': mdot,
        'At': At,
        'Dt': Dt,
        'Ae': Ae,
        'De': De,
        'D_chamber': D_chamber,
        'L_chamber': L_chamber * 100.0,
        'V_chamber': V_chamber,
        'Lstar': Lstar_cylinder,
        'V_total': V_total,
        'Lstar_with_converge': Lstar_with_converge,
        'L_throat': L_throat * 100.0,
        'L_nozzle': L_nozzle * 100.0
    }

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*70)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" gamma = {r['gamma']:.4f} | eps = {r['eps']:.3f}")
    print(f" Tc = {r['Tc']:.1f} K")
    print(f" Isp = {r['Isp']:.2f} s | Ve = {r['Ve']:.1f} m/s | c* = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dt = {r['Dt']*100:.2f} cm | De = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | L = {r['L_chamber']:.2f} cm")
    print(f" Chamber volume = {r['V_chamber']*1e6:.2f} cm³")
    print(f" L* (cylinder only) = {r['Lstar']:.3f} m")
    print(f" L* (with converge) = {r['Lstar_with_converge']:.3f} m")
    print(f" Throat length = {r['L_throat']:.2f} cm | Nozzle length = {r['L_nozzle']:.2f} cm")
    print("="*70 + "\n")

# =====================================================
# === MAIN EXECUTION ==================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (auto-calculated L*) ===\n")
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N

    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print(r)

        Ve_list, mdot_list, cstar_list, gamma_list, eps_list = [], [], [], [], []
        At_list, Ae_list, Dt_list, De_list, Tc_list, Isp_list = [], [], [], [], [], []

        for of in OF_RANGE:
            r = compute_sizing(target_thrust_N, pc, of)
            Ve_list.append(r['Ve'])
            mdot_list.append(r['mdot'])
            cstar_list.append(r['cstar'])
            gamma_list.append(r['gamma'])
            eps_list.append(r['eps'])
            At_list.append(r['At'])
            Ae_list.append(r['Ae'])
            Dt_list.append(r['Dt']*100)
            De_list.append(r['De']*100)
            Tc_list.append(r['Tc'])
            Isp_list.append(r['Isp'])

        Isp_arr, mdot_arr, Tc_arr = np.array(Isp_list), np.array(mdot_list), np.array(Tc_list)
        idx_isp_opt = np.nanargmax(Isp_arr)
        opt_Isp_OF = OF_RANGE[idx_isp_opt]

        print("\n=== OPTIMAL RESULTS SUMMARY ===")
        print(f"At Pc = {pc} psi")
        print("-" * 50)
        print(f"Optimal Isp at O/F = {opt_Isp_OF:.3f}")
        r_opt_isp = compute_sizing(target_thrust_N, pc, opt_Isp_OF)
        pretty_print(r_opt_isp)

        # ====== Performance Plots ======
        plt.figure(figsize=(10, 8))
        plt.suptitle(f"Performance at Pc={pc} psi")

        plt.subplot(2, 2, 1)
        plt.plot(OF_RANGE, Ve_list, label="Ve (m/s)")
        plt.title("Exhaust Velocity")
        plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 2)
        plt.plot(OF_RANGE, mdot_list, label="mdot (kg/s)")
        plt.title("Mass Flow Rate")
        plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 3)
        plt.plot(OF_RANGE, cstar_list, label="c* (m/s)")
        plt.title("Characteristic Velocity")
        plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.subplot(2, 2, 4)
        plt.plot(OF_RANGE, gamma_list, label="Gamma (γ)")
        plt.title("Gamma vs O/F")
        plt.xlabel("O/F"); plt.grid(True); plt.legend()

        plt.tight_layout()
        plt.show()

        # ====== Geometry & Temp Plots ======
        plt.figure(figsize=(10, 8))
        plt.suptitle(f"Geometry & Temp at Pc={pc} psi")

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
        ax1.set_xlabel("O/F"); ax1.set_ylabel("Chamber Temperature (K)", color='red')
        ax1.tick_params(axis='y', labelcolor='red'); ax1.grid(True)

        ax2 = ax1.twinx()
        ax2.plot(OF_RANGE, Isp_list, color='blue', label="Isp (s)")
        ax2.set_ylabel("Specific Impulse (s)", color='blue')
        ax2.tick_params(axis='y', labelcolor='blue')

        ax1.set_title("Temperature & Isp vs O/F")
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc="best")

        plt.tight_layout()
        plt.show()
