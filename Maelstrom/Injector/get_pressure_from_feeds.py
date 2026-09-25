import numpy as np
from pathlib import Path
from ruamel.yaml import YAML
from rocketcea.cea_obj import CEA_Obj

# Unit conversions
psi_to_pa = 6894.76
ft_to_m = 0.3048
m_to_in = 39.3701

# Oxygen gas constant
R_UNIV = 8.314462618
M_O2 = 0.031998
R_O2 = R_UNIV / M_O2  # J/(kg*K)

# ---------------------------------------------------------
# Helper: find YAML
# ---------------------------------------------------------
def find_yaml(filename="Maelstrom.yaml", start_dir=None):
    start = Path(start_dir or Path.cwd())
    for path in start.rglob(filename):
        return path
    raise FileNotFoundError(f"{filename} not found")


# ---------------------------------------------------------
# Internal GOX injector pressure solver
# ---------------------------------------------------------
def _required_gox_feed_pressure(m_dot, A, Cd, T, gamma, p_chamber,
                                assume_choked=True, max_iter=50, tol=1e-5):

    # Try choked solution first
    if assume_choked:
        Cg = np.sqrt(gamma / (R_O2 * T)) * (2/(gamma + 1)) ** ((gamma + 1)/(2*(gamma - 1)))
        p_feed = m_dot / (Cd * A * Cg)
        p_cr = (2/(gamma + 1)) ** (gamma/(gamma - 1))

        if p_chamber / p_feed < p_cr:
            return p_feed  # valid choked solution

    # Subsonic fallback (Newton iteration)
    p = max(p_chamber * 1.01, p_chamber + 5e3)
    for _ in range(max_iter):
        dp = max(p - p_chamber, 1.0)
        rho = p / (R_O2 * T)

        f = Cd * A * np.sqrt(2 * rho * dp) - m_dot
        if abs(f) < tol * m_dot:
            break

        # Numerical derivative
        p2 = p + max(p * 1e-4, 100.0)
        dp2 = max(p2 - p_chamber, 1.0)
        rho2 = p2 / (R_O2 * T)
        f2 = Cd * A * np.sqrt(2 * rho2 * dp2) - m_dot

        df_dp = (f2 - f) / (p2 - p)
        if df_dp == 0:
            break

        p -= f / df_dp
        p = max(p_chamber * 1.001, p)

    return p


# ---------------------------------------------------------
# MAIN INVERSE SOLVER FUNCTION
# ---------------------------------------------------------
def solve_chamber_conditions_from_feed_pressures(
    p_fu_feed_psi,
    p_ox_feed_psi,
    OF_range=(0.5, 2.5),
    Pc_range=(100, 500),
    n_OF=200,
    n_Pc=200
):
    """
    Returns the Pc [psi] and OF that produce the required
    fuel and oxidizer feed pressures given injector geometry.
    """

    # -----------------------------------------------------
    # Load YAML
    # -----------------------------------------------------
    yaml_path = find_yaml()
    yaml = YAML()
    yaml.preserve_quotes = True
    with open(yaml_path, "r") as f:
        data = yaml.load(f)

    # -----------------------------------------------------
    # Extract geometry from YAML
    # -----------------------------------------------------
    throat_diameter_m = data["throat_diameter"] * 0.0254
    A_t = np.pi * (throat_diameter_m/2)**2

    # Fuel injector total area
    d_fu = data["rp_orifice_dia"] * 0.0254
    A_fu = np.pi * (d_fu/2)**2 * data["rp_orifice_number"]

    # Ox injector total area
    d_ox = data["gox_orifice_dia"] * 0.0254
    A_ox = np.pi * (d_ox/2)**2 * data["gox_orifice_number"]

    Cd_fu = data["rp_discharge_coeff"]
    Cd_ox = data["gox_discharge_coeff"]

    rho_fu = 810.0  # you already assume this constant

    T_ox = 293.15
    gamma_ox = 1.4

    cstar_eff = data["cstar_eff"]

    # Convert target pressures
    p_fu_target_pa = p_fu_feed_psi * psi_to_pa
    p_ox_target_pa = p_ox_feed_psi * psi_to_pa

    # -----------------------------------------------------
    # Setup CEA
    # -----------------------------------------------------
    engine = CEA_Obj(oxName="GOX", fuelName="JetA")

    # -----------------------------------------------------
    # Search space
    # -----------------------------------------------------
    OFs = np.linspace(OF_range[0], OF_range[1], n_OF)
    Pcs = np.linspace(Pc_range[0], Pc_range[1], n_Pc)

    best_err = 1e99
    best_Pc = None
    best_OF = None

    # -----------------------------------------------------
    # Brute-force scan
    # -----------------------------------------------------
    for OF in OFs:
        for Pc in Pcs:

            p_c_pa = Pc * psi_to_pa

            # --- CEA c* ---
            try:
                cstar = engine.get_Cstar(Pc=Pc, MR=OF)
            except Exception:
                continue

            cstar_real = cstar * cstar_eff
            if cstar_real <= 0 or not np.isfinite(cstar_real):
                continue

            # m_dot = Pc * A_t / c*   (RocketCEA uses Pc in psi but c* returns ft/s)
            m_dot_total = p_c_pa * A_t / (cstar_real * ft_to_m)
            if m_dot_total <= 0:
                continue

            m_dot_fu = m_dot_total / (1 + OF)
            m_dot_ox = m_dot_total - m_dot_fu

            # --- Fuel feed pressure ---
            dp_fu = (m_dot_fu / (Cd_fu * A_fu))**2 / (2 * rho_fu)
            p_fu_req = p_c_pa + dp_fu

            # --- Ox feed pressure ---
            p_ox_req = _required_gox_feed_pressure(
                m_dot=m_dot_ox,
                A=A_ox,
                Cd=Cd_ox,
                T=T_ox,
                gamma=gamma_ox,
                p_chamber=p_c_pa
            )

            if p_ox_req <= p_c_pa:
                continue

            # Error metric
            err = ((p_fu_req - p_fu_target_pa)/p_fu_target_pa)**2 + \
                  ((p_ox_req - p_ox_target_pa)/p_ox_target_pa)**2

            if err < best_err:
                best_err = err
                best_Pc = Pc
                best_OF = OF

    return best_Pc, best_OF