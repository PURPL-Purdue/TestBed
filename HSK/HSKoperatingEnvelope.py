import numpy as np
import matplotlib.pyplot as plt
import math
from rocketcea.cea_obj_w_units import CEA_Obj   # changed to w_units

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

    cp = None  # not used downstream

    return {
        'Tc': Tc,
        'gamma': gamma,
        'MolWt': MolWt,
        'R_specific': R_specific,
        'cp': cp
    }


def CEA(F_lbf, of, pc_psi):

    data = get_numbers_extended(of, pc_psi)

    gamma, Tc, R_specific = data['gamma'], data['Tc'], data['R_specific']

    Pc = pc_psi                # removed PSI_TO_PA conversion (CEA expects psia)
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

    Dt = math.sqrt(4.0 * At_m2 / math.pi)

    # mdot from definition: mdot = Pc*At/c*
    # Convert Pc (psi) → Pa and At → m^2 for SI consistency
    mdot = (pc_psi * 6894.757 * At_m2) / cstar  # kg/s (removed 0.45359237 factor)

    Ae = At_m2 * eps
    De = math.sqrt(4.0 * Ae / math.pi)

    Ve = Cf * cstar * EFFICIENCY_FACTOR  # m/s
    Isp = Ve / G0                        # seconds

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
        'R_specific': R_specific
    }

    return result


def OFflowChecker(OF, mdot, FireTime):
    mdot_ox = mdot * OF / (1 + OF)
    mdot_fuel = mdot / (1 + OF)

    m_ox = mdot_ox * FireTime
    m_fuel = mdot_fuel * FireTime

    passfail = int(m_ox < 48.5 and m_fuel < 16.45)

    return passfail, mdot_ox, mdot_fuel


def plotter(result, FireTime, passfail, thrust_lbf):

    row = np.array([
        passfail,
        result['pc_psi'],          # psi
        result['of'],
        result['mdot'],            # kg/s
        thrust_lbf,                # lbf
        result['Tc'],              # K
        result['Dt'] * 100,        # cm
        result['Dt'] * 200         # chamber diameter (cm)
    ])

    return row


if __name__ == "__main__":

    ThrustStart = 1500
    ThrustMax = 1800
    ThrustTemp = ThrustStart
    Pcstart = 300
    PcMax = 325
    OFstart = 1.5
    OFTemp = OFstart
    OFEnd = 2.5
    FireTime = 3

    print("HSK Operating Envelope Trade Study")
    print("----------------------------------")
    print("PassFail | Pc (psi) | OF Ratio | Total Mass Flow (kg/s) | Thrust (lbf) | Chamber Temp (K) | Throat Diameter (cm) | Chamber Diameter (cm)")

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

                passfail, mdot_ox, mdot_fuel = OFflowChecker(OFstart, result['mdot'], FireTime)

                row = plotter(result, FireTime, passfail, ThrustStart)
                rows.append(row)

                print(f"{int(row[0])} | {row[1]:.2f} | {row[2]:.2f} | {row[3]:.2f} | {row[4]:.2f} | {row[5]:.2f} | {row[6]:.2f} | {row[7]:.2f}")

                OFstart += 0.5
            ThrustStart += 20
        Pcstart += 5

    data_matrix = np.vstack(rows)

    show_opt = input("Save? (y/n): ").strip().lower()

    if show_opt == "y":
        np.savetxt("trade_study_matrix.csv", data_matrix, delimiter=",")
