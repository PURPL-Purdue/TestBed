import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "RP1"

TARGET_THRUST_LBF = 1000.0
CHAMBER_PRESSURES = [200]          # Single-Pc runs still supported

OF_TARGET = 1.5
OF_RANGE = np.linspace(0.5, 3.0, 50)

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = 1
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
THROAT_LENGTH_FACTOR = 0.5

# === Chamber pressure sweep configuration ===
PC_SWEEP_START_PSI = 100.0
PC_SWEEP_END_PSI   = 600.0
PC_SWEEP_POINTS    = 6

PC_SWEEP_RANGE = np.linspace(
    PC_SWEEP_START_PSI,
    PC_SWEEP_END_PSI,
    PC_SWEEP_POINTS
)

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
    term = (pe_pc) ** (-(gamma - 1.0) / gamma)
    return math.sqrt((2.0 / (gamma - 1.0)) * (term - 1.0))

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
    MolWt, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    R_specific = R_UNIVERSAL / MolWt
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)
    return {'cstar': cstar, 'Tc': Tc, 'gamma': gamma}

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc = data['gamma'], data['cstar'], data['Tc']

    Pc_Pa = pc_psi * PSI_TO_PA
    Pa_Pa = AMBIENT_P_PSI * PSI_TO_PA
    Pa_over_Pc = Pa_Pa / Pc_Pa

    Me = Mach_from_pe_pc(Pa_over_Pc, gamma)
    eps = area_ratio_from_M(Me, gamma)
    Pe_Pc = pe_over_pc_from_M(Me, gamma)

    Cf = compute_Cf_ideal(gamma, Pe_Pc, Pa_over_Pc, eps)
    At = F_newtons / (Pc_Pa * Cf)
    mdot = Pc_Pa * At / cstar

    Isp = (Cf * cstar) / G0
    Ve = Isp * G0

    Ae = At * eps
    Dt = math.sqrt(4 * At / math.pi)
    De = math.sqrt(4 * Ae / math.pi)

    return {
        'pc_psi': pc_psi, 'of': of,
        'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp,
        'Ve': Ve, 'Cf': Cf,
        'eps': eps, 'mdot': mdot,
        'At': At, 'Ae': Ae,
        'Dt': Dt, 'De': De
    }

# =====================================================
# === MAIN =============================================
# =====================================================

if __name__ == "__main__":
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N

    # === Pressure sweep storage ===
    pc_sweep_Isp = {}
    pc_sweep_Ve = {}
    pc_sweep_mdot = {}
    pc_sweep_eps = {}

    # === Chamber pressure sweep ===
    for pc in PC_SWEEP_RANGE:
        Ve_list, mdot_list, Isp_list, eps_list = [], [], [], []
        for of in OF_RANGE:
            r = compute_sizing(target_thrust_N, pc, of)
            Ve_list.append(r['Ve'])
            mdot_list.append(r['mdot'])
            Isp_list.append(r['Isp'])
            eps_list.append(r['eps'])

        pc_sweep_Ve[pc] = Ve_list
        pc_sweep_mdot[pc] = mdot_list
        pc_sweep_Isp[pc] = Isp_list
        pc_sweep_eps[pc] = eps_list

    # =====================================================
    # === SWEPT CHAMBER PRESSURE PLOTS ====================
    # =====================================================

    plt.figure(figsize=(10, 8))
    plt.suptitle("Effect of Chamber Pressure on Performance")

    plt.subplot(2, 2, 1)
    for pc in PC_SWEEP_RANGE:
        plt.plot(OF_RANGE, pc_sweep_Isp[pc], label=f"{pc:.0f} psi")
    plt.title("Isp vs O/F")
    plt.xlabel("O/F"); plt.ylabel("Isp (s)")
    plt.grid(True); plt.legend()

    plt.subplot(2, 2, 2)
    for pc in PC_SWEEP_RANGE:
        plt.plot(OF_RANGE, pc_sweep_Ve[pc], label=f"{pc:.0f} psi")
    plt.title("Exhaust Velocity vs O/F")
    plt.xlabel("O/F"); plt.ylabel("Ve (m/s)")
    plt.grid(True); plt.legend()

    plt.subplot(2, 2, 3)
    for pc in PC_SWEEP_RANGE:
        plt.plot(OF_RANGE, pc_sweep_mdot[pc], label=f"{pc:.0f} psi")
    plt.title("Mass Flow vs O/F")
    plt.xlabel("O/F"); plt.ylabel("ṁ (kg/s)")
    plt.grid(True); plt.legend()

    plt.subplot(2, 2, 4)
    for pc in PC_SWEEP_RANGE:
        plt.plot(OF_RANGE, pc_sweep_eps[pc], label=f"{pc:.0f} psi")
    plt.title("Expansion Ratio vs O/F")
    plt.xlabel("O/F"); plt.ylabel("ε")
    plt.grid(True); plt.legend()

    plt.tight_layout()
    plt.show()
