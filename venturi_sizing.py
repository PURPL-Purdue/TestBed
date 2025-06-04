# Maelstrom Venturi Sizing Code
# Authors: Dominik Sloup
# First Created: 06/01/2025
# Last Updated: 06/01/2025
# Calculations done in SI units

import numpy as np
import matplotlib.pyplot as plt
import CEA_Wrap as CEA
from matplotlib.widgets import TextBox
from pyfluids import Fluid, FluidsList, Input

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

# rho = pyfluids.Fluid(FluidsList) # Water density at ambient pressure and temp
C_d_vent = 0.905 # Discharge coefficient guess (based on PSP Liquids' testing)
# https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/936903267/Successful+Water+Flow+1
C_d_inj = 0.8 # Injector discharge coefficient guess
p_1 = 680 # Pressure upstream of cavitating venturi (psi)

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

p_1_pa = p_1 * 6894.76 # Pressure upstream of cavitating venturi (Pa)
m_dot_kg = 0.961 # Desired mass flow (kg/s)
A_f_m2 = 3.2e-05 # Combined fuel injection area (m^2)

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

water = Fluid(FluidsList.Water).with_state(Input.pressure(p_1_pa), Input.temperature(20))
rho = water.density # Water density (kg/m^3)

p_sat_pa = Fluid(FluidsList.Water).with_state(Input.quality(0.0), Input.temperature(20)).pressure # Water saturation pressure at 20 degrees celsius
A_m2 = m_dot_kg / (C_d_vent * np.sqrt(2 * rho * (p_1_pa - p_sat_pa))) # Throat area (m^2)
d_m = 2 * np.sqrt(A_m2 / np.pi) # Venturi throat diameter (m)

K = ((C_d_vent * A_m2) / (C_d_inj * A_f_m2)) ** 2
p_line = (K * p_1 + 625)/(1 + K)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

d = d_m * 39.3701 # Venturi throat diameter (in)

print(f"Anticipated C_d: {C_d_vent}")
print(f"Upstream venturi pressure: {p_1} psi\n")
print(f"Downstream venturi pressure: {p_line} psi\n")
print(f"Venturi throat diameter: {round(d,3)} in")