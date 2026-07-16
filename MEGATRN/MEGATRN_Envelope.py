import numpy as np
import math
from rocketcea.cea_obj_w_units import CEA_Obj

# Universal variables
MIN_CONTRACTION_RATIO = 2.0
MAX_CONTRACTION_RATIO = 4.0
LSTAR = 1.2  # meters

CHAMBER_DIA_CM = 10.7
CHAMBER_DIA_M = CHAMBER_DIA_CM / 100.0

AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85

OXIDIZER = "LOX"
FUEL = "RP1"
R_UNIVERSAL = 8314.4621


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
        'R_specific': R_specific,
    }


def CEA(F_lbf, of, pc_psi):
    data = get_numbers_extended(of, pc_psi)
    gamma, Tc = data['gamma'], data['Tc']

    Pc = pc_psi
    Pa = AMBIENT_P_PSI
    Pa_over_Pc = Pa / Pc

    if Pa_over_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_over_Pc, gamma)
        eps = area_ratio_from_M(Me, gamma)
        Pe_Pc = pe_over_pc_from_M(Me, gamma)
    else:
        Me, eps, Pe_Pc = 1.0, 1.0, 1.0

    Cf = compute_Cf_ideal(gamma, Pe_Pc, Pa_over_Pc, eps)
    At = F_lbf / (Pc * Cf)
    cstar = cea.get_Cstar(Pc=pc_psi, MR=of)

    At_m2 = At * 0.00064516
    Dt = math.sqrt(4.0 * At_m2 / math.pi)
    mdot = (pc_psi * 6894.757 * At_m2) / cstar

    result = {
        'pc_psi': pc_psi,
        'of': of,
        'Tc': Tc,
        'cstar': cstar,
        'mdot': mdot,
        'Dt': Dt,
        'eps': eps,
    }

    return result


def plotter(result, thrust_lbf, chamber_length):
    row = np.array([
        1,
        result['pc_psi'],
        result['of'],
        result['mdot'],
        thrust_lbf,
        result['Tc'],
        result['Dt'] * 100,
        ((CHAMBER_DIA_CM / (result['Dt'] * 100))**2),  # contraction AREA ratio
        result['eps'],
        chamber_length * 100  # Save chamber length in cm
    ])
    return row


if __name__ == "__main__":
    ThrustStart = 100
    ThrustMax = 1800
    ThrustTemp = ThrustStart
    Pcstart = 180
    PcMax = 205
    OFstart = 0.5
    OFTemp = OFstart
    OFEnd = 4

    print("Cali Operating Envelope Trade Study")
    print("----------------------------------")
    print("PassFail | Pc (psi) | OF Ratio | Total Mass Flow (kg/s) | Thrust (lbf) | Chamber Temp (K) | Throat Diameter (cm) | Contraction Area Ratio | EPS | Chamber Length (cm)")

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
                
                # Check 1: Contraction Area Ratio Limits
                contraction_area_ratio = (CHAMBER_DIA_CM / (result['Dt'] * 100)) ** 2
                is_ratio_valid = MIN_CONTRACTION_RATIO <= contraction_area_ratio <= MAX_CONTRACTION_RATIO

                if is_ratio_valid:
                    # Check 2: Lstar length constraint (L_chamber >= 1.25 * Chamber_Diameter)
                    chamber_length = LSTAR / contraction_area_ratio
                    is_length_valid = chamber_length >= (1.25 * CHAMBER_DIA_M)

                    if is_length_valid:
                        row = plotter(result, ThrustStart, chamber_length)
                        rows.append(row)
                        print(f"{int(row[0])} | {row[1]} | {row[2]:.2f} | {row[3]:.2f} | {row[4]:.2f} | {row[5]:.2f} | {row[6]:.2f} | {row[7]:.2f} | {row[8]:.2f} | {row[9]:.2f}")

                OFstart += 0.5
            ThrustStart += 20
        Pcstart += 5

    if rows:
        data_matrix = np.vstack(rows)
        show_opt = input("Save? (y/n): ").strip().lower()
        if show_opt == "y":
            np.savetxt("LOXRPenvelope.csv", data_matrix, delimiter=",")
    else:
        print("No configurations met the contraction ratio and L* chamber length constraints.")