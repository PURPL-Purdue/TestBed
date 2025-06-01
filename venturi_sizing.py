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
C_d_inj = 0.7 # Injector discharge coefficient guess
m_dot = 0.58875 # Desired fuel mass flow (lb/s)
p_1 = 650 # Pressure upstream of cavitating venturi (psi)
p_c = 500 # Nominal chamber pressure (psi)
A_f = 6 * 3 * np.pi * (0.034 / 2) ** 2 # Combined fuel injection area (in^2)
# TODO: Need to calculate the actual line pressure given the mass flow in from the sonic nozzle and out the injectro assuming steady state

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

p_1_pa = p_1 * 6894.76 # Pressure upstream of cavitating venturi (Pa)
p_c_pa = p_c * 6894.76
m_dot_kg = m_dot * 0.453592 # Desired mass flow (kg/s)
A_f_m2 = A_f * 0.00064516 # Combined fuel injection area (m^2)

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

water = Fluid(FluidsList.Water).with_state(Input.pressure(p_1_pa), Input.temperature(20))
rho = water.density # Water density (kg/m^3)

p_sat_pa = Fluid(FluidsList.Water).with_state(Input.quality(0.0), Input.temperature(20)).pressure # Water saturation pressure at 20 degrees celsius
A_m2 = m_dot_kg / (C_d_vent * np.sqrt(2 * rho * (p_1_pa - p_sat_pa))) # Throat area (m^2)
d_m = 2 * np.sqrt(A_m2 / np.pi) # Venturi throat diameter (m)

K = ((C_d_vent * A_m2) / (C_d_inj * A_f_m2)) ** 2

p_line = (K * p_1 + p_c)/(1 + K)

# mass flow without the venturi:

m_dot_orig = 0.8 * A_f_m2 * np.sqrt(2 * rho * p_1_pa - p_c_pa)

print(m_dot_orig / 0.453592)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

d = d_m * 39.3701 # Venturi throat diameter (in)

print(f"Anticipated C_d: {C_d_vent}")
print(f"Upstream venturi pressure: {p_1} psi\n")
print(f"Line pressure: {p_line} psi\n")
print(f"Chamber pressure: {p_c} psi\n")

print(f"Venturi throat diameter: {d} in")