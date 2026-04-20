import numpy as np
from rocketcea.cea_obj import CEA_Obj
from scipy.interpolate import LinearNDInterpolator

# Unit conversions
psi_to_pa = 6894.76
ft_to_m = 0.3048

# Gas constants
R_UNIV = 8.314462618
M_O2 = 0.031998
M_H2 = 0.002016
R_O2 = R_UNIV / M_O2   # J/(kg·K)
R_H2 = R_UNIV / M_H2   # J/(kg·K)
gamma_O2 = 1.4
gamma_H2 = 1.4


def finite_minmax(a):
    a = np.asarray(a, dtype=float)
    m = np.isfinite(a)
    if not np.any(m):
        raise ValueError("Array has no finite values.")
    return float(np.min(a[m])), float(np.max(a[m]))


def required_feed_pressure_orifice(
        m_dot, CdA, T, gamma, p_chamber, R, max_iter=80, rtol=1e-8):
    """
    Find upstream feed pressure (PSI) required to deliver m_dot (kg/s)
    through an orifice of effective area CdA (m²) into a chamber at
    p_chamber (PSI), using ideal-gas density.

    Args:
        m_dot:      required mass flow rate [kg/s]
        CdA:        discharge coefficient * orifice area [m²]
        T:          propellant temperature [K]
        gamma:      specific heat ratio [-]
        p_chamber:  chamber pressure [PSI]
        R:          specific gas constant [J/(kg·K)]
        max_iter:   bisection iteration limit
        rtol:       relative tolerance on m_dot

    Returns:
        upstream feed pressure [PSI]
    """
    if m_dot <= 0:
        return float(p_chamber)

    pr_crit = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))

    def get_density(p_psi):
        """Ideal-gas density [kg/m³] at given pressure and fixed temperature."""
        p_pa = p_psi * psi_to_pa
        return p_pa / (R * T)

    def mdot_of_p0(p0_psi):
        """Mass flow through orifice given upstream pressure in PSI."""
        if p0_psi <= p_chamber:
            return 0.0

        p0_pa = p0_psi * psi_to_pa
        p_chamber_pa = p_chamber * psi_to_pa
        rho = get_density(p0_psi)
        pr = p_chamber_pa / p0_pa  # dimensionless

        if pr <= pr_crit:
            # Choked (sonic)
            return CdA * np.sqrt(
                gamma * rho * p0_pa
                * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0))
            )
        else:
            # Subsonic compressible
            term = pr ** (2.0 / gamma) - pr ** ((gamma + 1.0) / gamma)
            if term <= 0.0:
                return 0.0
            return CdA * np.sqrt(
                2.0 * rho * p0_pa * (gamma / (gamma - 1.0)) * term
            )

    # --- Bracket the solution ---
    p_lo_psi = p_chamber * 1.000001
    f_lo = mdot_of_p0(p_lo_psi) - m_dot  # expected < 0

    p_hi_psi = max(p_chamber * 1.1, p_chamber * 2.0)
    f_hi = mdot_of_p0(p_hi_psi) - m_dot

    for _ in range(60):
        if f_hi >= 0:
            break
        p_hi_psi *= 1.5
        f_hi = mdot_of_p0(p_hi_psi) - m_dot
    else:
        raise RuntimeError(
            f"Could not bracket feed pressure for m_dot={m_dot:.4f} kg/s. "
            "Check CdA, temperature, and chamber pressure."
        )

    # --- Bisection ---
    for _ in range(max_iter):
        p_mid_psi = 0.5 * (p_lo_psi + p_hi_psi)
        f_mid = mdot_of_p0(p_mid_psi) - m_dot

        if abs(f_mid) <= rtol * m_dot:
            return float(p_mid_psi)

        if f_mid > 0:
            p_hi_psi = p_mid_psi
        else:
            p_lo_psi = p_mid_psi

    return float(0.5 * (p_lo_psi + p_hi_psi))


def build_inverse_interpolators(fu_p_map, ox_p_map, pc_axis, of_axis, fill_value=np.nan):
    """
    Builds inverse interpolators:
      (p_fu, p_ox) -> Pc
      (p_fu, p_ox) -> OF
    """
    PC, OF = np.meshgrid(pc_axis, of_axis, indexing="ij")

    pf = fu_p_map.ravel()
    po = ox_p_map.ravel()
    pc_vals = PC.ravel()
    of_vals = OF.ravel()

    mask = np.isfinite(pf) & np.isfinite(po) & np.isfinite(pc_vals) & np.isfinite(of_vals)
    pts = np.column_stack([pf[mask], po[mask]])

    pc_interp = LinearNDInterpolator(pts, pc_vals[mask], fill_value=fill_value)
    of_interp = LinearNDInterpolator(pts, of_vals[mask], fill_value=fill_value)
    return pc_interp, of_interp


def pressure_map(minOF, maxOF, maxPc, resolution, cstar_eff,
                 throat_area, fuel_CdA, ox_CdA):

    engine = CEA_Obj(oxName="GOX", fuelName="GH2")

    pc_scale = np.linspace(1.0, maxPc, resolution)   # PSI
    OF_scale = np.linspace(minOF, maxOF, resolution)  # -

    fu_p_map = np.full((resolution, resolution), np.nan, dtype=float)
    ox_p_map = np.full((resolution, resolution), np.nan, dtype=float)

    for i, pc in enumerate(pc_scale):
        for j, OF in enumerate(OF_scale):
            cstar = engine.get_Cstar(Pc=pc, MR=OF)   # ft/s
            p_c_pa = pc * psi_to_pa
            cstar_real = cstar * cstar_eff * ft_to_m  # m/s

            m_dot_total = p_c_pa * throat_area / cstar_real  # kg/s
            m_dot_fu = m_dot_total / (1.0 + OF)
            m_dot_ox = m_dot_total - m_dot_fu

            fu_p_map[i, j] = required_feed_pressure_orifice(
                m_dot_fu, fuel_CdA, T=293.15, gamma=gamma_H2,
                p_chamber=pc, R=R_H2
            )
            ox_p_map[i, j] = required_feed_pressure_orifice(
                m_dot_ox, ox_CdA, T=293.15, gamma=gamma_O2,
                p_chamber=pc, R=R_O2
            )

        print(f"Chamber pressure: {pc:.1f} psi")

    # Build inverse interpolators
    pc_interp, of_interp = build_inverse_interpolators(fu_p_map, ox_p_map, pc_scale, OF_scale)

    # Inverse map grid
    pf_min, pf_max = finite_minmax(fu_p_map)
    po_min, po_max = finite_minmax(ox_p_map)
    pf_axis = np.linspace(pf_min, pf_max, resolution)
    po_axis = np.linspace(po_min, po_max, resolution)
    PF, PO = np.meshgrid(pf_axis, po_axis, indexing="ij")

    PC_pfpo = pc_interp(PF, PO)
    OF_pfpo = of_interp(PF, PO)

    # Mask edge artifacts on all four boundaries
    bad = (
        (PC_pfpo >= 0.999 * maxPc) |
        (PC_pfpo <= 1.001 * 1.0) |
        (OF_pfpo < 1.001 * minOF) |
        (OF_pfpo > 0.999 * maxOF)
    )
    PC_pfpo[bad] = np.nan
    OF_pfpo[bad] = np.nan

    return pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo