import pandas as pd
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Constants
orifice_dia = 0.144  # inches
rho = 997  # kg/m³ (water at 25°C)
in_to_m = 0.0254
psi_to_Pa = 6894.76

orifice_dia_m = orifice_dia * in_to_m
expected_cd = 0.905
throat_area = np.pi * (orifice_dia_m / 2) ** 2
expected_cda = expected_cd * throat_area
print(f"Expected CdA (m²): {expected_cda:.6e}")

timestep_1 = 1.8
timestep_2 = 0.9

# Test-specific mass collected (g) and fixed durations (s)
mass_data = {
    1: (162, timestep_1), 3: (467, timestep_1), 4: (490, timestep_1), 5: (480, timestep_1),
    6: (328, timestep_1), 7: (331, timestep_1), 8: (282, timestep_1),
    9: (278, timestep_2), 10: (277, timestep_2), 11: (279, timestep_2), 12: (399, timestep_2),
    13: (202, timestep_2), 14: (356, timestep_2), 15: (355, timestep_2), 16: (355, timestep_2),
    17: (350, timestep_2), 18: (434, timestep_2)
}

data_dir = "/Users/dominiksloup/Desktop/Test raw"
files = sorted(glob(os.path.join(data_dir, "VENTURI_test_low_*.csv")))

results = []

for i, file in enumerate(files, start=1):
    if i not in mass_data:
        continue

    try:
        df = pd.read_csv(file)
        df.columns = df.columns.str.strip()

        # Filter valve states: True, False, True
        mask = (df["SN-N2-04 State"] == True) & \
               (df["SN-N2-07 State"] == False) & \
               (df["SN-FU-01 State"] == True)
        filtered = df[mask]

        if len(filtered) < 4:
            print(f"Test {i}: not enough valid rows in valve state window")
            continue

        filtered = filtered.iloc[1:-1]

        delta_p_psi = filtered["PT-VT-01 Pressure"] - filtered["PT-VT-02 Pressure"]
        delta_p_Pa = delta_p_psi.mean() * psi_to_Pa

        m_grams, duration = mass_data[i]
        m_dot = (m_grams / 1000) / duration  # kg/s

        results.append({
            "Test": i,
            "Avg ΔP (Pa)": delta_p_Pa,
            "Mass Flow Rate (kg/s)": m_dot
        })

    except Exception as e:
        print(f"Error in test {i}: {e}")
        continue
    
result_df = pd.DataFrame(results)

sqrt_dp = np.sqrt(result_df["Avg ΔP (Pa)"])

dp_theoretical = np.linspace(result_df["Avg ΔP (Pa)"].min(),
                             result_df["Avg ΔP (Pa)"].max(), 50)

slope, intercept, r_value, p_value, std_err = linregress(sqrt_dp, result_df["Mass Flow Rate (kg/s)"])

regressed_cda = slope / np.sqrt(2 * rho)
regressed_cd = regressed_cda / throat_area

sqrt_dp_theoretical = np.sqrt(dp_theoretical)
mdot_theoretical = expected_cda * np.sqrt(2 * rho * dp_theoretical)

print(f"Actual CdA (m²): {regressed_cda:.6e}")

print(f"Expected Cd: {expected_cd:.3f}")
print(f"Actual Cd: {regressed_cd:.3f}")


plt.figure(figsize=(10, 6))

plt.scatter(sqrt_dp, result_df["Mass Flow Rate (kg/s)"], label="Measured data", color='blue')
plt.plot(sqrt_dp_theoretical, mdot_theoretical, label="Theoretical model", color='red')

mdot_fit = slope * sqrt_dp_theoretical + intercept
plt.plot(sqrt_dp_theoretical, mdot_fit, label="Fitted line", color='green', linestyle='--')

plt.xlabel("Sqrt(Avg ΔP) (Pa$^{0.5}$)")
plt.ylabel("Mass Flow Rate (kg/s)")
plt.title("Mass Flow Rate vs. Sqrt of Pressure Drop")
plt.legend()
plt.grid(True)

# Equation text with R²
equation_text = f"$y = {slope:.5f} \\, \\sqrt{{\\Delta P}} + {intercept:.3f}$\n$R^2 = {r_value**2:.3f}$"

# Place text near the upper left or a good location
x_text = sqrt_dp_theoretical.min() + 0.1 * (sqrt_dp_theoretical.max() - sqrt_dp_theoretical.min())
y_text = mdot_fit.max() * 0.8

plt.text(x_text, y_text, equation_text, fontsize=12, color='green',
         bbox=dict(facecolor='white', alpha=0.7))

plt.tight_layout()
plt.show()
