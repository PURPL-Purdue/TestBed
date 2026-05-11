import numpy as np
import matplotlib.pyplot as plt
import math
# Using cea_obj_si to get SI units directly (m, kg, s, K, J, etc.)
from rocketcea.cea_obj import CEA_Obj, cea_obj_si 
import bisect
from matplotlib.widgets import Slider

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "ETHANOL"

TARGET_THRUST_LBF = 900.0
CHAMBER_PRESSURES = [200.0]
OF_TARGET = 1.0
OF_RANGE = np.linspace(0.5, 3.0, 50)
ExpansionRatio = 2.704

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = 1
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
    Wrap cea_obj_si calls to get SI values directly.
    """
    # Using SI version of the object
    cea = cea_obj_si.CEA_Obj(oxName=ox, fuelName=fuel)

    # Temperatures (K)
    Tc = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)[0]

    # MolWt (kg/kmol) and gamma
    MolWt, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    R_specific = R_UNIVERSAL / MolWt

    # Transport Properties (at Throat index 1)
    # Returns (Cp, Visc, Cond, Pr)
    trans_props = cea.get_Transport_Ppts(Pc=pc_psi, MR=OF_ratio, eps=1)
    cp = trans_props[0][1]      # J/kg-K
    mu = trans_props[1][1]      # Pa-s (kg/m-s)
    k_gas = trans_props[2][1]   # W/m-K
    prandtl = trans_props[3][1] # Dimensionless

    # Density and Velocity (at Throat index 1)
    rho_throat = cea.get_Densities(Pc=pc_psi, MR=OF_ratio, eps=1)[1] # kg/m3
    sonic_vel_throat = cea.get_SonicVelocities(Pc=pc_psi, MR=OF_ratio, eps=1)[1] # m/s
    u_throat = sonic_vel_throat * 1.0 # Mach 1 at throat

    # cstar
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)

    return {
        'cstar': cstar,
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
        'cp': cp,
        'mu': mu,
        'k_gas': k_gas,
        'prandtl': prandtl,
        'rho_throat': rho_throat,
        'u_throat': u_throat
    }

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
    
    # Heat transfer components
    cp = data['cp']
    mu = data['mu']
    k_gas = data['k_gas']
    prandtl = data['prandtl']
    rho_throat = data['rho_throat']
    u_throat = data['u_throat']

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
    V_total = V_chamber + ((1/3) * math.pi * L_converge * ((D_chamber/2)**2 + (D_chamber/2)*(Dt/2) + (Dt/2)**2))
    Lstar_with_converge = V_total / At

    result = {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber': L_chamber*100.0,
        'L_chamber_full': (L_chamber + L_converge)*100.0, 'V_chamber': V_chamber,
        'L_throat': (0.382 * Dt + 0.5 * Dt)*100.0, 'L_nozzle': ((De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG))))*100.0,
        'cp': cp, 'mu': mu, 'k_gas': k_gas, 'prandtl': prandtl, 
        'rho_throat': rho_throat, 'u_throat': u_throat, 'R_specific': R_specific
    }
    return result

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*60)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" Chamber Temp (gas) = {r['Tc']:.1f} K")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Full chamber L = {r['L_chamber_full']:.2f} cm")
    print(f" --- Heat Transfer Gas Props (Throat) ---")
    print(f" Prandtl: {r['prandtl']:.4f} | Conductivity: {r['k_gas']:.4f} W/m-K")
    print(f" Viscosity: {r['mu']:.4e} Pa-s | Density: {r['rho_throat']:.3f} kg/m³")
    print("="*60 + "\n")

# =====================================================
# === MAIN =============================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (SI Heat Transfer Additions) ===\n")
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N
    show_opt = input("Also show optimal O/F results? (y/n): ").strip().lower()

    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print(r)

        # Plotting logic remains the same...
        Isp_list = [compute_sizing(target_thrust_N, pc, of)['Isp'] for of in OF_RANGE]
        idx_isp = np.nanargmax(Isp_list)
        opt_Isp_OF = OF_RANGE[idx_isp]

        if show_opt == "y":
            print("=== OPTIMAL RESULTS SUMMARY ===")
            pretty_print(compute_sizing(target_thrust_N, pc, opt_Isp_OF))