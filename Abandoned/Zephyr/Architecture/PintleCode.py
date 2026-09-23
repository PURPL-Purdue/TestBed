import math
from rocketcea.cea_obj import CEA_Obj

# =====================================================
# ================== CONSTANTS ========================
# =====================================================

G0 = 9.80665                 # m/s^2
PSI_TO_PA = 6894.75729
PA_TO_PSI = 1.0 / PSI_TO_PA
LBF_TO_N = 4.4482216152605

EFFICIENCY_FACTOR = 0.95     # c* efficiency
AMBIENT_P_PSI = 14.7
AMBIENT_P_PA = AMBIENT_P_PSI * PSI_TO_PA

OXIDIZER = "LOX"
FUEL = "RP-1"

# =====================================================
# =============== ISENTROPIC RELATIONS ================
# =====================================================

def pe_over_pc_from_M(M, gamma):
    return (1 + 0.5*(gamma-1)*M*M)**(-gamma/(gamma-1))

def area_ratio_from_M(M, gamma):
    term = (2/(gamma+1))*(1 + 0.5*(gamma-1)*M*M)
    exp = (gamma+1)/(2*(gamma-1))
    return (1/M)*(term**exp)

def Mach_from_pe_pc(pe_pc, gamma):
    M = 3.0
    for _ in range(40):
        f = pe_over_pc_from_M(M, gamma) - pe_pc
        df = (
            pe_over_pc_from_M(M+1e-6, gamma)
            - pe_over_pc_from_M(M-1e-6, gamma)
        ) / 2e-6
        M -= f/df
    return max(M, 0.01)

# =====================================================
# =============== THRUST COEFFICIENT ==================
# =====================================================

def compute_Cf(gamma, Pe_Pc, Pa_Pc, Ae_At):
    term1 = math.sqrt(
        (2*gamma**2/(gamma-1)) *
        (2/(gamma+1))**((gamma+1)/(gamma-1)) *
        (1 - Pe_Pc**((gamma-1)/gamma))
    )
    term2 = (Pe_Pc - Pa_Pc) * Ae_At
    return term1 + term2

# =====================================================
# ============ mdot + At → Pc SOLVER ==================
# =====================================================

def solve_pc_from_mdot(
    mdot,
    At,
    of,
    pc_guess_psi=300.0,
    tol=1e-5,
    max_iter=40
):
    """
    Fixed inputs:
      - mdot (kg/s)
      - throat area At (m^2)
      - mixture ratio (O/F)

    Solves:
      - chamber pressure Pc
      - thrust (output)
    """

    cea = CEA_Obj(oxName=OXIDIZER, fuelName=FUEL)
    pc_psi = pc_guess_psi

    for _ in range(max_iter):

        # === CEA call ===
        cstar_ft_s = cea.get_Cstar(
            Pc=pc_psi,
            MR=of
        )
        gamma = cea.get_Gamma(
            Pc=pc_psi,
            MR=of,
            eps=40.0
        )

        cstar = cstar_ft_s * 0.3048  # ft/s → m/s
        cstar_eff = cstar * EFFICIENCY_FACTOR

        Pc = pc_psi * PSI_TO_PA

        # === mdot implied by Pc ===
        mdot_calc = Pc * At / cstar_eff
        err = mdot_calc - mdot

        if abs(err/mdot) < tol:
            break

        Pc *= mdot / mdot_calc
        pc_psi = Pc * PA_TO_PSI

    # =================================================
    # === NOZZLE PERFORMANCE (OUTPUTS) ================
    # =================================================

    Pa_Pc = AMBIENT_P_PA / Pc

    if Pa_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_Pc, gamma)
        eps = area_ratio_from_M(Me, gamma)
        Pe_Pc = pe_over_pc_from_M(Me, gamma)
    else:
        eps = 1.0
        Pe_Pc = 1.0

    Cf = compute_Cf(gamma, Pe_Pc, Pa_Pc, eps)

    thrust = Cf * Pc * At
    Isp = (Cf * cstar) / G0

    return {
        "Pc_psi": pc_psi,
        "Pc_Pa": Pc,
        "mdot": mdot,
        "At": At,
        "Dt": math.sqrt(4*At/math.pi),
        "gamma": gamma,
        "cstar": cstar,
        "Cf": Cf,
        "eps": eps,
        "thrust_N": thrust,
        "thrust_lbf": thrust / LBF_TO_N,
        "Isp": Isp
    }

# =====================================================
# ==================== MAIN ===========================
# =====================================================

if __name__ == "__main__":

    # -------- FIXED INPUTS --------
    MDOT = 12.0                 # kg/s
    OF = 2.6
    THROAT_DIAMETER_MM = 70.0   # FIXED HARDWARE
    At = math.pi * (THROAT_DIAMETER_MM/1000)**2 / 4

    result = solve_pc_from_mdot(
        mdot=MDOT,
        At=At,
        of=OF,
        pc_guess_psi=400.0
    )

    print("\n=== Fixed Throat Area, mdot → Pc ===\n")
    print(f"Mass flow        : {result['mdot']:.3f} kg/s")
    print(f"Throat diameter  : {result['Dt']*1000:.2f} mm")
    print(f"Chamber pressure : {result['Pc_psi']:.1f} psi")
    print(f"Gamma            : {result['gamma']:.3f}")
    print(f"C*               : {result['cstar']:.1f} m/s")
    print(f"Cf               : {result['Cf']:.3f}")
    print(f"Expansion ratio  : {result['eps']:.1f}")
    print(f"Thrust           : {result['thrust_lbf']:.1f} lbf")
    print(f"Isp              : {result['Isp']:.1f} s")
