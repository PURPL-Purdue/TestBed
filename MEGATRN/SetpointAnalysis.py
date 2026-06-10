import numpy as np
import math
from rocketcea.cea_obj import CEA_Obj  # Using default native object for bug-free transport properties


# =====================================================
# === USER CONFIGURATION & CONSTANTS ==================
# =====================================================
OXIDIZER = "LOX"
FUEL = "ETHANOL"
TARGET_THRUST_LBF = 900.0
CHAMBER_PRESSURES = [200.0]
OF_TARGET = 1.0
OF_RANGE = np.linspace(0.5, 3.0, 50)
AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = 1
LSTAR = 1.0922

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822
R_UNIVERSAL = 8314.4621

# Accurate manual conversion factors from NASA native units
MILLIPOISE_TO_PAS = 0.0001        
MCAL_CM_K_S_TO_W_MK = 0.4184      
CAL_G_K_TO_J_KG_K = 4184.0        

# =====================================================
# === MATHEMATICAL HELPERS ============================
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
# === CORE GEOMETRY & SIZING ENGINE ===================
# =====================================================
def compute_sizing(F_newtons, pc_psi, of):
    # Step 1: Geometry Sizing (Using standard infinite chamber baseline)
    cea_init = CEA_Obj(oxName=OXIDIZER, fuelName=FUEL)
    _, gamma_init = cea_init.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=of)

    Pc_Pa = pc_psi * PSI_TO_PA
    Pa_Pa = AMBIENT_P_PSI * PSI_TO_PA
    Pa_over_Pc = Pa_Pa / Pc_Pa

    if Pa_over_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_over_Pc, gamma_init)
        eps = area_ratio_from_M(Me, gamma_init)
        Pe_Pc = pe_over_pc_from_M(Me, gamma_init)
    else:
        Me, eps, Pe_Pc = 1.0, 1.0, 1.0

    Cf = compute_Cf_ideal(gamma_init, Pe_Pc, Pa_over_Pc, eps)
    At = F_newtons / (Pc_Pa * Cf)
    Dt = math.sqrt(4.0 * At / math.pi)
    D_chamber = 2.0 * Dt
    A_chamber = (math.pi / 4.0) * (D_chamber ** 2)
    
    # Calculate exact contraction ratio
    contraction_ratio = A_chamber / At  

    # Step 2: Create the Finite Area Combustor (FAC) object
    cea_fac = CEA_Obj(oxName=OXIDIZER, fuelName=FUEL, fac_CR=contraction_ratio)

    # Step 3: Extract and isolate temperatures and densities safely
    ch_temp, th_temp, _ = cea_fac.get_Temperatures(Pc=pc_psi, MR=of, eps=eps)
    ch_rho_raw, th_rho_raw, _ = cea_fac.get_Densities(Pc=pc_psi, MR=of, eps=eps)
    ch_sonic_raw, th_sonic_raw, _ = cea_fac.get_SonicVelocities(Pc=pc_psi, MR=of, eps=eps)

    # Step 4: Isolate Mach Numbers without using the broken get_MachNumber unpacker
    # For an FAC engine, chamber Mach is solved via contraction ratio; throat is always choked (1.0)
    ch_mach = cea_fac.get_Chamber_MachNumber(Pc=pc_psi, MR=of, fac_CR=contraction_ratio)
    th_mach = 1.0

    # Step 5: Convert native metrics to pure SI scalars
    ch_vel_m_s = (ch_sonic_raw * ch_mach) * 0.3048
    ch_rho_m_s = ch_rho_raw * 16.0185
    
    th_vel_m_s = th_sonic_raw * th_mach * 0.3048 
    th_rho_m_s = th_rho_raw * 16.0185

    # Step 6: Get Singular Transport Properties
    ch_trans = cea_fac.get_Chamber_Transport(Pc=pc_psi, MR=of)
    ch_cp = ch_trans[0] * CAL_G_K_TO_J_KG_K
    ch_mu = ch_trans[1] * MILLIPOISE_TO_PAS
    ch_k  = ch_trans[2] * MCAL_CM_K_S_TO_W_MK
    ch_pr = ch_trans[3]

    th_trans = cea_fac.get_Throat_Transport(Pc=pc_psi, MR=of)
    th_cp = th_trans[0] * CAL_G_K_TO_J_KG_K
    th_mu = th_trans[1] * MILLIPOISE_TO_PAS
    th_k  = th_trans[2] * MCAL_CM_K_S_TO_W_MK
    th_pr = th_trans[3]

    _, ch_gamma = cea_fac.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=of, eps=eps)
    _, th_gamma = cea_fac.get_Throat_MolWt_gamma(Pc=pc_psi, MR=of, eps=eps)

    # Performance analytics calculations
    MolWt, gamma = cea_fac.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=of, eps=eps)
    R_specific = R_UNIVERSAL / MolWt
    cstar = cstar_from_T_gamma_R(ch_temp, gamma, R_specific)

    cstar_eff = cstar * EFFICIENCY_FACTOR
    mdot = Pc_Pa * At / cstar_eff
    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Isp = (Cf * cstar) / G0
    Ve = Isp * G0 * EFFICIENCY_FACTOR

    V_chamber = LSTAR * At
    L_chamber = V_chamber / A_chamber
    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    L_chamber_full = (L_chamber + L_converge) * 100.0

    return {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'cstar': cstar, 'Isp': Isp, 'Ve': Ve,
        'Cf': Cf, 'eps': eps, 'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber_full': L_chamber_full, 'contraction_ratio': contraction_ratio,
        
        # Chamber Data
        'ch_temp': ch_temp, 'ch_density': ch_rho_m_s, 'ch_vel': ch_vel_m_s,
        'ch_cp': ch_cp, 'ch_mu': ch_mu, 'ch_k': ch_k, 'ch_pr': ch_pr, 'ch_gamma': ch_gamma,
        # Throat Data
        'th_temp': th_temp, 'th_density': th_rho_m_s, 'th_vel': th_vel_m_s,
        'th_cp': th_cp, 'th_mu': th_mu, 'th_k': th_k, 'th_pr': th_pr, 'th_gamma': th_gamma,
    }
# =====================================================
# === DISPLAY FORMATTER ===============================
# =====================================================
def pretty_print_fac(r):
    print("="*60)
    print(f"FAC CONFIG: Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" Calculated Contraction Area Ratio (fac_CR): {r['contraction_ratio']:.4f}")
    print(f" Isp = {r['Isp']:.3f} s | Cstar = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Full chamber L = {r['L_chamber_full']:.2f} cm")
    print("-" * 60)
    print(f" --- 1D INPUT STRINGS: CHAMBER STATION ---")
    print(f"GasTemp = {r['ch_temp']:.2f}                # K")
    print(f"GasDensity = {r['ch_density']:.4f}              # kg/m^3")
    print(f"GasVelocity = {r['ch_vel']:.2f}            # m/s (True Subsonic Flow)")
    print(f"GasDynamicViscosity = {r['ch_visc']:.4e} # Pa-s") if 'ch_visc' in r else print(f"GasDynamicViscosity = {r['ch_mu']:.4e} # Pa-s")
    print(f"GasGamma = {r['ch_gamma']:.4f}           # Cp/Cv")
    print(f"GasSpecificHeat = {r['ch_cp']:.1f}        # J/kg-K")
    print(f"GasThermalConductivity = {r['ch_k']:.4f} # W/m-K")
    print(f"GasPr = {r['ch_pr']:.4f}                  # Dimensionless")
    print("-" * 60)
    print(f" --- 1D INPUT STRINGS: THROAT STATION ---")
    print(f"GasTemp = {r['th_temp']:.2f}                # K")
    print(f"GasDensity = {r['th_density']:.4f}              # kg/m^3")
    print(f"GasVelocity = {r['th_vel']:.2f}            # m/s (Choked Flow)")
    print(f"GasDynamicViscosity = {r['th_mu']:.4e} # Pa-s")
    print(f"GasGamma = {r['th_gamma']:.4f}           # Cp/Cv")
    print(f"GasSpecificHeat = {r['th_cp']:.1f}        # J/kg-K")
    print(f"GasThermalConductivity = {r['th_k']:.4f} # W/m-K")
    print(f"GasPr = {r['th_pr']:.4f}                  # Dimensionless")
    print("="*60 + "\n")


# =====================================================
# === MAIN =============================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (Finite Area Combustor Wrapper) ===\n")
    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N
    
    show_opt = input("Also show optimal O/F results? (y/n): ").strip().lower()

    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print_fac(r)

        if show_opt == "y":
            Isp_list = [compute_sizing(target_thrust_N, pc, of)['Isp'] for of in OF_RANGE]
            idx_isp = np.nanargmax(Isp_list)
            opt_Isp_OF = OF_RANGE[idx_isp]
            print("=== OPTIMAL RESULTS SUMMARY ===")
            pretty_print_fac(compute_sizing(target_thrust_N, pc, opt_Isp_OF))