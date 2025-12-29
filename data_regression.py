import os
from pathlib import Path

import numpy as np
import pandas as pd
from pyfluids import Fluid, FluidsList, Input
from rocketcea.cea_obj import CEA_Obj

from Injector.get_mdot_from_area import get_mdot_from_area


def data_regression(
    file_name: str,
    fuel_choice: str,
    fuel_CdA: float,
    ox_CdA: float,
    A_star: float,
    *,
    csv_dir: str | Path = r"C:\Users\Precision\Downloads",
    ox_name: str = "GOX",
    T_C: float = 20.0,
    accept_frac: float = 0.10,   # median ±10% band
    ATM_PSI: float = 14.6959,
    psi_to_pa: float = 6894.757,
    ft_to_m: float = 0.3048,
    # column names
    state_ox: str = "SN-O2-02 State",
    state_fu: str = "SN-FU-01 State",
    fuel_upstream_col: str = "PT-FU-01 Pressure",
    ox_upstream_col: str = "PT-O2-04 Pressure",
    downstream_col: str = "PT-GG-01 Pressure",   # used for mdot + CEA Pc
    downstream_col2: str = "PT-GG-02 Pressure",  # used for c*_act
    # fuel density
    rho_fu_const: float = 810.0,
    clear_terminal: bool = True,
    print_summary: bool = True,
):
    """
    Returns:
      summary: dict of mean values + counts + overall pressure deviation
      df_ss: steady-state dataframe (filtered)
    """

    # -------------------------
    # Load
    # -------------------------
    csv_path = Path(csv_dir) / file_name
    df = pd.read_csv(csv_path)

    # -------------------------
    # Filter: both states TRUE
    # -------------------------
    def _truthy(series: pd.Series) -> pd.Series:
        return (
            series.astype(str)
            .str.strip()
            .str.upper()
            .isin(["TRUE", "1", "T", "YES", "Y"])
        )

    df[state_ox] = _truthy(df[state_ox])
    df[state_fu] = _truthy(df[state_fu])

    df_active = df.loc[df[state_ox] & df[state_fu]].copy()

    # -------------------------
    # "Steady" filter: median ±accept_frac for key pressures
    # -------------------------
    steady_pressure_cols = [ox_upstream_col, fuel_upstream_col, downstream_col, downstream_col2]

    for c in steady_pressure_cols:
        df_active[c] = pd.to_numeric(df_active[c], errors="coerce")

    med = df_active[steady_pressure_cols].median(skipna=True)
    lo = (1.0 - accept_frac) * med
    hi = (1.0 + accept_frac) * med

    steady_mask = np.ones(len(df_active), dtype=bool)
    for c in steady_pressure_cols:
        steady_mask &= df_active[c].between(lo[c], hi[c], inclusive="both").to_numpy()

    df_ss = df_active.loc[steady_mask].copy()

    # -------------------------
    # Clear terminal + header
    # -------------------------
    if clear_terminal:
        os.system("cls" if os.name == "nt" else "clear")

    if print_summary:
        print(f"\n==={file_name}===")
        print(f"Active samples (both states TRUE): {len(df_active)}")
        print(f"Median±{int(accept_frac*100)}% steady-state samples: {len(df_ss)}")

    if len(df_ss) == 0:
        summary = {
            "file_name": file_name,
            "fuel_choice": fuel_choice,
            "n_active": int(len(df_active)),
            "n_steady": 0,
            "note": "No steady-state samples passed the median band filter.",
        }
        return summary, df_ss

    # -------------------------
    # Pressures (psig -> psia -> Pa)
    # -------------------------
    p_ox_pa    = (df_ss[ox_upstream_col].to_numpy()   + ATM_PSI) * psi_to_pa
    p_fu_pa    = (df_ss[fuel_upstream_col].to_numpy() + ATM_PSI) * psi_to_pa
    p_down_pa  = (df_ss[downstream_col].to_numpy()    + ATM_PSI) * psi_to_pa
    p_down2_pa = (df_ss[downstream_col2].to_numpy()   + ATM_PSI) * psi_to_pa

    # -------------------------
    # Fuel density (constant placeholder)
    # -------------------------
    df_ss["rho_fu"] = float(rho_fu_const)

    # -------------------------
    # Ox density + gamma (from pyfluids)
    # -------------------------
    rho_ox = np.empty(len(df_ss), dtype=float)
    gamma_ox = np.empty(len(df_ss), dtype=float)

    for i in range(len(df_ss)):
        f = Fluid(FluidsList.Oxygen).with_state(
            Input.pressure(float(p_ox_pa[i])),
            Input.temperature(T_C),
        )
        rho_ox[i] = f.density
        gamma_ox[i] = (f.sound_speed ** 2) * f.density / f.pressure

    df_ss["rho_ox"] = rho_ox
    df_ss["gamma_ox"] = gamma_ox

    # -------------------------
    # RocketCEA engine + flows + c*
    # -------------------------
    engine = CEA_Obj(oxName=ox_name, fuelName=fuel_choice)

    m_dot_fu = np.empty(len(df_ss), dtype=float)
    m_dot_ox = np.empty(len(df_ss), dtype=float)
    OF_ratio = np.empty(len(df_ss), dtype=float)
    cstar_the = np.empty(len(df_ss), dtype=float)

    for i in range(len(df_ss)):
        m_dot_fu[i] = get_mdot_from_area(
            "liquid", fuel_CdA, p_fu_pa[i], p_down_pa[i], df_ss["rho_fu"].iat[i]
        )
        m_dot_ox[i] = get_mdot_from_area(
            "gas", ox_CdA, p_ox_pa[i], p_down_pa[i], df_ss["rho_ox"].iat[i], gamma=float(gamma_ox[i])
        )
        OF_ratio[i] = m_dot_ox[i] / m_dot_fu[i]
        cstar_the[i] = engine.get_Cstar(Pc=float(p_down_pa[i]), MR=float(OF_ratio[i]))

    df_ss["m_dot_fu"] = m_dot_fu
    df_ss["m_dot_ox"] = m_dot_ox
    df_ss["OF_ratio"] = OF_ratio

    df_ss["cstar_the"] = cstar_the * ft_to_m                 # rocketcea often returns ft/s
    df_ss["cstar_act"] = p_down2_pa * A_star / (m_dot_fu + m_dot_ox)
    df_ss["cstar_eff"] = df_ss["cstar_act"] / df_ss["cstar_the"]

    # -------------------------
    # Pressure deviation metric (psi): mean abs dev from median, averaged across channels
    # -------------------------
    dev_psi = {}
    for c in steady_pressure_cols:
        x = df_ss[c].to_numpy(dtype=float)  # still psig
        x_med = np.nanmedian(x)
        dev_psi[c] = float(np.nanmean(np.abs(x - x_med)))

    overall_dev_psi = float(np.mean(list(dev_psi.values()))) if dev_psi else float("nan")

    # -------------------------
    # Summary stats
    # -------------------------
    summary = {
        "file_name": file_name,
        "fuel_choice": fuel_choice,
        "n_active": int(len(df_active)),
        "n_steady": int(len(df_ss)),
        "m_dot_fu_mean": float(df_ss["m_dot_fu"].mean()),
        "m_dot_ox_mean": float(df_ss["m_dot_ox"].mean()),
        "OF_mean": float(df_ss["OF_ratio"].mean()),
        "p_fu_mean_psig": float(df_ss[fuel_upstream_col].mean()),
        "p_ox_mean_psig": float(df_ss[ox_upstream_col].mean()),
        "pc_mean_psig": float(df_ss[downstream_col2].mean()),
        "overall_pressure_dev_psig": float(overall_dev_psi),
        "rho_fu_mean": float(df_ss["rho_fu"].mean()),
        "rho_ox_mean": float(df_ss["rho_ox"].mean()),
        "cstar_expected_mean_mps": float(df_ss["cstar_the"].mean()),
        "cstar_actual_mean_mps": float(df_ss["cstar_act"].mean()),
        "cstar_eff_mean": float(df_ss["cstar_eff"].mean()),
    }

    # -------------------------
    # Print summary (same vibe as your original)
    # -------------------------
    if print_summary:
        print("\n=== FLOWS ===")
        print(f"Fuel mdot           : {summary['m_dot_fu_mean']:.3g} kg/s")
        print(f"Ox mdot             : {summary['m_dot_ox_mean']:.3g} kg/s")
        print(f"OF Ratio            : {summary['OF_mean']:.3g}\n")

        print("=== PRESSURES ===")
        print(f"Fuel feed pressure  : {summary['p_fu_mean_psig']:.3g} psi")
        print(f"Ox feed pressure    : {summary['p_ox_mean_psig']:.3g} psi")
        print(f"Chamber pressure    : {summary['pc_mean_psig']:.3g} psi")
        print(f"Pressure deviation  : {summary['overall_pressure_dev_psig']:.1f} psi\n")

        print("=== DENSITIES ===")
        print(f"Fuel density        : {summary['rho_fu_mean']:.2f} kg/m^3")
        print(f"Ox density          : {summary['rho_ox_mean']:.2f} kg/m^3\n")

        print("=== C* PERFORMANCE ===")
        print(f"Expected C*         : {summary['cstar_expected_mean_mps']:.0f} m/s")
        print(f"Actual C*           : {summary['cstar_actual_mean_mps']:.0f} m/s")
        print(f"Cstar Efficiency    : {summary['cstar_eff_mean']:.3g}")

    return summary, df_ss