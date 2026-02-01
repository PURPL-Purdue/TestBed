import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj import CEA_Obj

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0
LSTAR = 1.0922

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822
R_UNIVERSAL = 8314.4621
KG_TO_LB = 2.20462
GAL_PER_LITER = 0.264172
L_PER_M3 = 1000.0
M3_TO_GAL = GAL_PER_LITER * L_PER_M3

OXIDIZER = "LOX"
FUEL = "RP1"
DENSITY_RP1 = 810.0
DENSITY_LOX = 1140.0

# unit conversion constants used for CEA outputs
BTU_TO_J = 1055.05585
LBM_TO_KG = 0.45359237
LBM_PER_CUFT_TO_KG_PER_M3 = 16.01846337396
FT_TO_M = 0.3048
MILLIPOISE_TO_PA_S = 1e-4


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

def get_numbers_extended(OF_ratio, pc_psi, ox, fuel):

    # temperatures (degR), convert primary chamber temperature to K
    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0  # Rankine -> Kelvin

    # molecular weight and gamma (MolWt returned in lbm/lbmole)
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt, gamma = mw_gamma

    # R_specific
    MolWt_kg_per_kmol = MolWt * LBM_TO_KG
    R_specific = R_UNIVERSAL / MolWt_kg_per_kmol


    # cstar
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)

    cp_btu_per_lbm_degR = cea.get_Chamber_Cp(Pc=pc_psi, MR=OF_ratio, eps=1)
    cp = cp_btu_per_lbm_degR * BTU_TO_J / LBM_TO_KG * (9.0 / 5.0)

    return {
        'cstar': cstar,
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
        'cp': cp
    }

def CEA(F_newtons, of, pc_psi):
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

    mdot = Pc_Pa * At / cstar

    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Ve = Cf * cstar * EFFICIENCY_FACTOR
    Isp = Ve / G0
    


    #R_throat = 0.382 * Dt  # Throat radius-of-curvature approximation (not geometric throat radius)
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    #c
    
    result = {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'L_nozzle': L_nozzle*100.0,
        'cp': cp, 'mu': mu, 'rho_chamber': rho_chamber,
        'sonic_vel': sonic_vel, 'R_specific': R_specific
    }


    return result

def OFflowChecker(OF, mdot, FireTime):
        mdot_ox = mdot * OF / (1 + OF)
        mdot_fuel = mdot / (1 + OF)

        m_ox = mdot_ox * FireTime
        m_fuel = mdot_fuel * FireTime

        #feasability check based on amount of propellant available
        passfail = int(m_ox < 48.5 and m_fuel < 16.45)

        return passfail, mdot_ox, mdot_fuel

def plotter(result, FireTime, passfail, thrust):
    """
    Creates a single matrix row for trade study / plotting
    """
    row = np.array([
        passfail,                  # 0 feasibility
        result['pc_psi'],          # 1 chamber pressure
        result['of'],              # 2 OF ratio
        result['mdot'],            # 3 total mass flow
        thrust,                    # 4 thrust
        result['Tc'],              # 5 chamber temperature
        result['Dt']*100,          # 6 throat diameter (cm)
        result['Dt']*200           # 7 chamber diameter (2xthroat)
    ])


    return row

    return row






if __name__ == "__main__":
    ThrustStart = 200
    ThrustMax = 5000
    ThrustTemp = ThrustStart
    Pcstart = 100
    PcMax = 300
    OFstart = .5
    OFTemp = OFstart
    OFEnd = 4
    FireTime = 2

    print("HSK Operating Envelope Trade Study")
    print("----------------------------------")
    print("PassFail | Pc (psi) | OF Ratio | Total Mass Flow (kg/s) | Thrust (N) | Chamber Temp (K) | Throat Diameter (cm) | Chamber Diameter (cm)")
    
    cea = CEA_Obj(oxName=OXIDIZER, fuelName=FUEL)
    rows = []

    while Pcstart <= PcMax:
        ThrustStart = ThrustTemp
        while ThrustStart <= ThrustMax:
            OFstart = OFTemp
            while OFstart <= OFEnd:
                result = CEA(ThrustStart, OFstart, Pcstart)
                
                passfail, mdot_ox, mdot_fuel = OFflowChecker(OFstart, result['mdot'], FireTime)
                
                row = plotter(result, FireTime, passfail, ThrustStart)
                rows.append(row)
                print(f"{int(row[0])}        | {row[1]:.2f}    | {row[2]:.2f}   | {row[3]:.2f}        | {row[4]:.2f}  | {row[5]:.2f}      | {row[6]:.2f}     | {row[7]:.4f}")

                OFstart += 0.5
            ThrustStart += 20
        Pcstart += 5

    data_matrix = np.vstack(rows)
    np.savetxt("trade_study_matrix.csv", data_matrix, delimiter=",")

