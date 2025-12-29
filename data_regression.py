import pandas as pd
import numpy as np
from pathlib import Path
from pyfluids import Fluid, FluidsList, Input
from Injector.get_mdot_from_area import get_mdot_from_area
from rocketcea.cea_obj import CEA_Obj
import os

# =================================================
# USER INPUTS
# =================================================

file_name = "GG_test_low_20250526_205511.csv"

# file_name = "GG_test_low_20250527_024632.csv"
# file_name = "GG_test_low_20250527_033233.csv"

csv_path = Path(r"C:\Users\Precision\Downloads") / file_name

pressure_cols = [
    "PT-GG-01 Pressure",
    "PT-GG-02 Pressure",
    "PT-O2-04 Pressure",
    "PT-FU-01 Pressure",
]

state_ox = "SN-O2-02 State"
state_fu = "SN-FU-01 State"

ref_pressure = "PT-O2-04 Pressure"

fuel_upstream_col = "PT-FU-01 Pressure"
ox_upstream_col   = "PT-O2-04 Pressure"
downstream_col    = "PT-GG-01 Pressure"   # used for mdot + CEA Pc (edit if needed)
downstream_col2   = "PT-GG-02 Pressure"   # used for c*_act (edit if needed)

psi_to_pa = 6894.757293168
pa_to_psi = 1.0 / psi_to_pa
ft_to_m = 0.3048
ATM_PSI = 14.6959  # set to 0.0 if sensors are absolute (psia)

engine = CEA_Obj(oxName="GOX", fuelName="JetA")

fuel_CdA = 7.21e-06
ox_CdA   = 1.87e-05
A_star   = 0.000160041
gamma_ox = 1.4

# "Steady-ish" detector (no fixed pressure cutoff)
window = 10
deriv_thresh = 2.0
std_thresh   = 3.0

# Acceptance interval around median (two-sided)
accept_frac = 0.1  # ±10%

T_C = 20.0  # pyfluids temp input in Celsius (per you)

# =================================================
# LOAD CSV
# =================================================
df = pd.read_csv(csv_path)

# =================================================
# FILTER: BOTH STATES TRUE
# =================================================
df[state_ox] = df[state_ox].astype(str).str.strip().str.upper().isin(["TRUE", "1", "T", "YES", "Y"])
df[state_fu] = df[state_fu].astype(str).str.strip().str.upper().isin(["TRUE", "1", "T", "YES", "Y"])

df_active = df.loc[df[state_ox] & df[state_fu]].copy()

# =================================================
# "STEADY" FILTER: median ±10% band (two-sided) across pressures you care about
# =================================================
accept_frac = 0.10  # 10% two-sided acceptance about the median

# choose which pressures define "steady" (use upstream/downstream that drive mdot/c*)
steady_pressure_cols = [ox_upstream_col, fuel_upstream_col, downstream_col, downstream_col2]

# ensure numeric
for c in steady_pressure_cols:
    df_active[c] = pd.to_numeric(df_active[c], errors="coerce")

# median-based band (robust to rugged data, no derivatives)
med = df_active[steady_pressure_cols].median(skipna=True)
lo = (1.0 - accept_frac) * med
hi = (1.0 + accept_frac) * med

steady_mask = np.ones(len(df_active), dtype=bool)
for c in steady_pressure_cols:
    steady_mask &= df_active[c].between(lo[c], hi[c], inclusive="both").to_numpy()

df_ss = df_active.loc[steady_mask].copy()

os.system('cls' if os.name == 'nt' else 'clear')

print("\n===" + file_name + "===")
print(f"Active samples (both states TRUE): {len(df_active)}")
print(f"Median±{int(accept_frac*100)}% steady-state samples: {len(df_ss)}")

# =================================================
# PRESSURES -> Pa (psig -> psia) using df_ss
# =================================================
p_ox_pa    = (df_ss[ox_upstream_col].to_numpy()    + ATM_PSI) * psi_to_pa
p_fu_up_pa = (df_ss[fuel_upstream_col].to_numpy()  + ATM_PSI) * psi_to_pa
p_down_pa  = (df_ss[downstream_col].to_numpy()     + ATM_PSI) * psi_to_pa
p_down2_pa = (df_ss[downstream_col2].to_numpy()    + ATM_PSI) * psi_to_pa

# =================================================
# rho / mdots / OF / cstar_eff (STEADY-ONLY)
# =================================================
df_ss["rho_fu"] = 810.0  # placeholder constant

rho_ox   = np.empty(len(df_ss), dtype=float)
gamma_ox = np.empty(len(df_ss), dtype=float)

for i in range(len(df_ss)):
    f = Fluid(FluidsList.Oxygen).with_state(
        Input.pressure(float(p_ox_pa[i])),  # Pa
        Input.temperature(T_C),              # °C
    )

    rho_ox[i] = f.density
    gamma_ox[i] = (f.sound_speed ** 2) * f.density / f.pressure

df_ss["rho_ox"] = rho_ox
df_ss["gamma_ox"] = gamma_ox

m_dot_fu = np.empty(len(df_ss), dtype=float)
m_dot_ox = np.empty(len(df_ss), dtype=float)
OF_ratio = np.empty(len(df_ss), dtype=float)
cstar_the = np.empty(len(df_ss), dtype=float)

for i in range(len(df_ss)):
    m_dot_fu[i] = get_mdot_from_area("liquid", fuel_CdA, p_fu_up_pa[i], p_down_pa[i], df_ss["rho_fu"].iat[i])
    m_dot_ox[i] = get_mdot_from_area("gas",   ox_CdA,   p_ox_pa[i],    p_down_pa[i], df_ss["rho_ox"].iat[i], gamma=gamma_ox[i])
    OF_ratio[i] = m_dot_ox[i] / m_dot_fu[i]
    cstar_the[i] = engine.get_Cstar(Pc=p_down_pa[i], MR=OF_ratio[i])

df_ss["m_dot_fu"] = m_dot_fu
df_ss["m_dot_ox"] = m_dot_ox
df_ss["OF_ratio"] = OF_ratio

df_ss["cstar_the"] = cstar_the * ft_to_m              # (m/s) if rocketcea returns (ft/s)
df_ss["cstar_act"] = p_down2_pa * A_star / (m_dot_fu + m_dot_ox)
df_ss["cstar_eff"] = df_ss["cstar_act"] / df_ss["cstar_the"]

# =================================================
# Averages you care about
# =================================================
OF_mean = float(df_ss["OF_ratio"].mean())
cstar_eff_mean = float(df_ss["cstar_eff"].mean())

# =================================================
# "average deviation of the pressure (psi)" for the steady set
# here: mean absolute deviation from the steady-set median, per pressure channel
# =================================================
dev_psi = {}
for c in steady_pressure_cols:
    x = df_ss[c].to_numpy(dtype=float)  # still psig in the CSV
    x_med = np.nanmedian(x)
    dev_psi[c] = float(np.nanmean(np.abs(x - x_med)))

overall_dev_psi = float(np.mean(list(dev_psi.values()))) if dev_psi else np.nan
# =================================================
# Print summary
# =================================================
print("\n=== FLOWS ===")
print(f"Fuel mdot           : {df_ss['m_dot_fu'].mean():.3g} kg/s")
print(f"Ox mdot             : {df_ss['m_dot_ox'].mean():.3g} kg/s")
print(f"OF Ratio            : {OF_mean:.3g}\n")

print("=== PRESSURES ===")
print(f"Fuel feed pressure  : {df_ss[fuel_upstream_col].mean():.3g} psi")
print(f"Ox feed pressure    : {df_ss[ox_upstream_col].mean():.3g} psi")
print(f"Chamber pressure    : {df_ss[downstream_col2].mean():.3g} psi\n")

print("=== DENSITIES ===")
print(f"Fuel density        : {df_ss['rho_fu'].mean():.2f} kg/m^3")
print(f"Ox density          : {df_ss['rho_ox'].mean():.2f} kg/m^3\n")

print("=== C* PERFORMANCE ===")
print(f"Expected C*         : {df_ss['cstar_the'].mean():.0f} m/s")
print(f"Actual C*           : {df_ss['cstar_act'].mean():.0f} m/s")
print(f"Cstar Efficiency    : {cstar_eff_mean:.3g}")
