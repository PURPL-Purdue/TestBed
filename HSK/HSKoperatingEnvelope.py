import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj_w_units import CEA_Obj
from scipy.integrate import solve_ivp
from thermal_solver import calculate_convection_coeff, heat_equation

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85

NOZZLE_HALF_ANGLE_DEG = 15.0
LSTAR = 1.0922

G0 = 9.80665

OXIDIZER = "LOX"
FUEL = "RP1"
DENSITY_RP1 = 810.0
DENSITY_LOX = 1140.0

R_UNIVERSAL = 8314.4621


# =================== PARAMETERS (THERMAL) ====================
#P_CHAMBER = 1.724e+6        # 250 psi in Pa
LENGTH = 0.2 * 0.0254       # m
THICKNESS = 0.01     # m

# ---------- Material properties ----------
material = "SS 316"
K_METAL = 16.3      # W/(m·K)  
C_METAL = 500       # J/(kg·K)   
RHO_METAL = 8000
ALPHA = K_METAL / (RHO_METAL * C_METAL)
T_MELT = 1370.0 + 273.15  # Celsius to Kelvin


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


def get_numbers_extended(OF_ratio, pc_psi):

    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0  # Rankine -> Kelvin

    MolWt, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)

    MolWt_kg_per_kmol = MolWt * 0.45359237  # lbm -> kg
    R_specific = R_UNIVERSAL / MolWt_kg_per_kmol

    return {
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
    }


def CEA(F_lbf, of, pc_psi):

    data = get_numbers_extended(of, pc_psi)

    gamma, Tc, R_specific = data['gamma'], data['Tc'], data['R_specific']

    PcPa = pc_psi * 6894.76
    Pc = pc_psi
    Pa = AMBIENT_P_PSI
    Pa_over_Pc = Pa / Pc

    # === Nozzle parameters ===
    if Pa_over_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_over_Pc, gamma)
        eps = area_ratio_from_M(Me, gamma)
        Pe_Pc = pe_over_pc_from_M(Me, gamma)
    else:
        Me, eps, Pe_Pc = 1.0, 1.0, 1.0

    Cf = compute_Cf_ideal(gamma, Pe_Pc, Pa_over_Pc, eps)

    # Throat area in in^2 (consistent with lbf and psi)
    At = F_lbf / (Pc * Cf)

    cstar = cea.get_Cstar(Pc=pc_psi, MR=of)  # now returns m/s (set in CEA object)

    # Removed in^2 → m^2 conversion (0.0006452)
    # Keep geometry in meters by converting once properly
    At_m2 = At * 0.00064516  # correct in^2 → m^2

    Dt = math.sqrt(4.0 * At_m2 / math.pi) # m

    # mdot from definition: mdot = Pc*At/c*
    # Convert Pc (psi) → Pa and At → m^2 for SI consistency
    mdot = (pc_psi * 6894.757 * At_m2) / cstar  # kg/s (removed 0.45359237 factor)

    Ae = At_m2 * eps
    De = math.sqrt(4.0 * Ae / math.pi)

    Ve = Cf * cstar * EFFICIENCY_FACTOR  # m/s
    Isp = Ve / G0                        # seconds
    
    V_c = mdot / (At_m2 * (pc_psi * 6894.757 / (8.314 * Tc))) # math.sqrt(gamma * R_specific * Tc)  # m/s
    ##Thermal Analysis
    # Convert Pc from psi to Pa for thermal analysis
    H_conv, T_gas = calculate_convection_coeff(OXIDIZER, FUEL, of, PcPa, Dt, V_c, eps)[0:2]
    # print(f"Convection Coefficient: {H_conv:.2f} W/m^2-K, Gas Temperature: {T_gas:.2f} K")
    T_initial = 25 + 273.15  # Celsius to Kelvin
    # ==================== SOLVE ====================
    n_nodes = 50
    r_nodes = np.linspace(Dt, Dt + THICKNESS, n_nodes)
    dr = r_nodes[1] - r_nodes[0]
    T0 = np.full(n_nodes, T_initial)
    t_final = 20 # seconds

    sol = solve_ivp(heat_equation, [0, t_final], T0, args=(r_nodes, dr, n_nodes, ALPHA, H_conv, T_gas, K_METAL), 
                    method='BDF', t_eval=np.linspace(0, t_final, 1000))

    # ==================== RESULTS ====================
    T_res = sol.y
    times = sol.t

    inner_wall_temp = T_res[1, :]

    # Find melting time
    melt_indices = np.where(inner_wall_temp >= T_MELT)[0]
    if len(melt_indices) > 0:
        t_melt = times[melt_indices[0]]
        # print(f"CRITICAL: Inner wall melts at {t_melt:.3f} seconds.")
    else:
        t_melt = -1
        # print("Wall did not melt within the time frame.")

    result = {
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
        'At': At_m2,
        'Dt': Dt,
        'Ae': Ae,
        'De': De,
        'R_specific': R_specific,
        'material': material,
        't_melt': t_melt
    }

    return result


def flowChecker(mdot_cea, maxThroatDia, Pc_psi, cstar):
    # Convert chamber pressure from psi to Pa
    Pc = Pc_psi * 6894.757

    # Throat area
    A_t = math.pi * maxThroatDia**2 / 4

    # Maximum allowable mass flow from rocket relation
    mdot_max = (Pc * A_t) / cstar

    # FAIL if requested mdot exceeds geometric limit
    passfail = int(mdot_cea < mdot_max)

    return passfail, mdot_max


def plotter(result, FireTime, passfail, thrust_lbf):

    row = np.array([
        passfail,
        result['pc_psi'],          # psi
        result['of'],
        result['mdot'],            # kg/s
        thrust_lbf,                # lbf
        result['Tc'],              # K
        result['Dt'] * 100,        # cm
        (10.7 / (result['Dt'] * 100)), #contraction Ratio
        result['eps'],
        result['t_melt'],           # seconds until melting (or -1 if no melt)
        result['Ve'],              # m/s
        #result['Dt'] * 200         # chamber diameter (cm)
    ])

    return row


# def thermalAnalysis(result, FireTime, passfail, thrust_lbf):

    

#     return row

if __name__ == "__main__":

    ThrustStart = 100
    ThrustMax = 1800
    ThrustTemp = ThrustStart
    Pcstart = 100
    PcMax = 250
    OFstart = 0.5
    OFTemp = OFstart
    OFEnd = 4
    FireTime = 3
    MaxThroatDia = 0.0535  #5.35 cm

    print("HSK Operating Envelope Trade Study")
    print("----------------------------------")
    print("PassFail | Pc (psi) | OF Ratio | Total Mass Flow (kg/s) | Thrust (lbf) | Chamber Temp (K) | Throat Diameter (cm) | EPS | Time to Melt (s) | Velocity")

    cea = CEA_Obj(
        oxName=OXIDIZER,
        fuelName=FUEL,
        pressure_units='psia',
        cstar_units='m/s'
    )

    rows = []

    while Pcstart <= PcMax:
        ThrustStart = ThrustTemp
        while ThrustStart <= ThrustMax:
            OFstart = OFTemp
            while OFstart <= OFEnd:

                result = CEA(ThrustStart, OFstart, Pcstart)

                passfail, maxMdot = flowChecker(result['mdot'], MaxThroatDia, Pcstart, result['cstar'])
                #passfail = OFflowChecker(OFstart, result['mdot'], MaxThroatDia)

                if passfail == 1:
                    row = plotter(result, FireTime, passfail, ThrustStart)
                    rows.append(row)

                    print(f"{int(row[0])} | {row[1]} | {row[2]:.2f} | {row[3]:.2f} | {row[4]:.2f} | {row[5]:.2f} | {row[6]:.2f} | {row[7]:.2f} | {row[8]:.2f} | {row[9]:.2f} | {row[10]:.2f}")

                OFstart += 0.5
            ThrustStart += 20
        Pcstart += 5

    data_matrix = np.vstack(rows)

    show_opt = input("Save? (y/n): ").strip().lower()

    if show_opt == "y":
        np.savetxt("LOXRPenvelope.csv", data_matrix, delimiter=",")
