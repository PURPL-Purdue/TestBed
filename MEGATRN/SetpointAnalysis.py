import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj_w_units import CEA_Obj 
import bisect
from matplotlib.widgets import Slider

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

OXIDIZER = "LOX"
FUEL = "RP1"

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
# === CONSTANTS & CONVERSIONS =========================
# =====================================================

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822
PSI_TO_MPA = 0.00689476  
R_UNIVERSAL = 8314.4621

# Manual conversion helpers for transport properties
MILLIPOISE_TO_PAS = 0.0001        
MCAL_CM_K_S_TO_W_MK = 0.4184      

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
    Using the 'cea_obj_w_units' wrapper.
    """
    cea = CEA_Obj(oxName=ox, fuelName=fuel, 
                  pressure_units='MPa', 
                  temperature_units='K',
                  density_units='kg/m^3',
                  specific_heat_units='J/kg-K')

    pc_mpa = pc_psi * PSI_TO_MPA

    # Temperatures (K)
    Tc = cea.get_Temperatures(Pc=pc_mpa, MR=OF_ratio, eps=1)[0]

    # MolWt and gamma
    MolWt, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_mpa, MR=OF_ratio)
    R_specific = R_UNIVERSAL / MolWt

    # Transport Properties (Chamber)
    trans_props = cea.get_Chamber_Transport(Pc=pc_mpa, MR=OF_ratio)
    cp = trans_props[0]                         # J/kg-K
    mu = trans_props[1] * MILLIPOISE_TO_PAS     # Pa-s
    k_gas = trans_props[2] * MCAL_CM_K_S_TO_W_MK # W/m-K
    prandtl = trans_props[3]                    

    # Density and Sonic Velocity at Throat
    rho_throat = cea.get_Densities(Pc=pc_mpa, MR=OF_ratio, eps=1)[1] 
    u_throat = cea.get_SonicVelocities(Pc=pc_mpa, MR=OF_ratio, eps=1)[1] # Mach 1 Velocity

    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)

    return {
        'cstar': cstar, 'Tc': Tc, 'gamma': gamma, 'MolWt': MolWt,
        'R_specific': R_specific, 'cp': cp, 'mu': mu, 'k_gas': k_gas,
        'prandtl': prandtl, 'rho_throat': rho_throat, 'u_throat': u_throat
    }

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
    
    cp, mu, k_gas = data['cp'], data['mu'], data['k_gas']
    prandtl, rho_throat, u_throat = data['prandtl'], data['rho_throat'], data['u_throat']

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

    cstar_eff = cstar * EFFICIENCY_FACTOR
    mdot = Pc_Pa * At / cstar_eff

    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Isp = (Cf * cstar) / G0
    Ve = Isp * G0 * EFFICIENCY_FACTOR

    D_chamber = 2.0 * Dt
    V_chamber = LSTAR * At
    L_chamber = V_chamber / (math.pi/4 * D_chamber**2)

    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    L_chamber_full = (L_chamber + L_converge) * 100.0

    result = {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber_full': L_chamber_full,
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
    print(f" isp = {r['Isp']:.3f} s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Full chamber L = {r['L_chamber_full']:.2f} cm")
    print(f" --- Heat Transfer Gas Props (Chamber) ---")
    print(f" Specific Heat (Cp): {r['cp']:.1f} J/kg-K")
    print(f" Prandtl: {r['prandtl']:.4f} | Conductivity: {r['k_gas']:.4f} W/m-K")
    print(f" Viscosity: {r['mu']:.4e} Pa-s | Density: {r['rho_throat']:.3f} kg/m³")
    print(f" Throat Velocity (Mach 1): {r['u_throat']:.1f} m/s")
    print(f" Cstar : {r['cstar']:.1f} m/s")
    print(f" Gamma : {r['gamma']:.4f}")
    print("="*60 + "\n")

# =====================================================
# === MAIN =============================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (Official SI Units Wrapper) ===\n")
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N
    
    show_opt = input("Also show optimal O/F results? (y/n): ").strip().lower()

    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print(r)

        if show_opt == "y":
            Isp_list = [compute_sizing(target_thrust_N, pc, of)['Isp'] for of in OF_RANGE]
            idx_isp = np.nanargmax(Isp_list)
            opt_Isp_OF = OF_RANGE[idx_isp]
            print("=== OPTIMAL RESULTS SUMMARY ===")
            pretty_print(compute_sizing(target_thrust_N, pc, opt_Isp_OF))