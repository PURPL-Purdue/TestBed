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

C_d = 0.905 # Discharge coefficient guess (based on PSP Liquids' testing) https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/936903267/Successful+Water+Flow+1
p_1 = 650 # Pressure downstream of cavitating venturi (psi)
x = 0.9 # Anticipated pressure recovery across the venturi
m_dot_kg = 0.961 # Desired mass flow (kg/s)
rho = 810 # RP-1 density (from internet :) )
p_0 = p_1 / x # Pressure upstream of the cavitating venturi (psi)
p_0_pa = p_0 * 6894.76 # Pressure upstream of cavitating venturi (Pa)
d = 0.426 # Line inner diameter (in)
A = np.pi * (d / 2) ** 2 # Cross sectional are of a -8 tube (in^2)
A_m2 = A / 1550 # Cross sectional are of a -8 tube (m^2)

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

# rho = Fluid(FluidsList.Water).with_state(Input.pressure(p_0_pa), Input.temperature(20)).density
p_sat_pa = Fluid(FluidsList.Water).with_state(Input.quality(0.0), Input.temperature(20)).pressure # Water saturation pressure at 20 degrees celsius
# This should be changed because RP-1's saturation pressure is different

# We assume cavitation to calculate the venturi throat area
A_star_m2 = m_dot_kg / (C_d * np.sqrt(2 * rho * (p_0_pa - p_sat_pa))) # Throat area (m^2)

# Modified bernoulli equation used to calculate throat area
# A_star_m2 = 1 / np.sqrt( (1 / A_m2**2) - (2 * rho * (p_sat_pa - p_0_pa)) / m_dot_kg**2 )

d_star_m = 2 * np.sqrt(A_star_m2 / np.pi) # Venturi throat diameter (m)
d_star_in = d_star_m * 39.3701 # Venturi throat diameter (in)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

print(f"\nAnticipated C_d:                 {C_d}")
print(f"Desired mass flow:               {m_dot_kg:.3f} kg/s")

print(f"\nUpstream venturi pressure:       {p_0:.0f} psi")
print(f"Expected throat pressure:        {(p_sat_pa / 6894.76):.1f} psi")
print(f"Downstream venturi pressure:     {p_1:.0f} psi")

print(f"Line inner diameter:             {0.426:.3f} in")  # assuming -8 tube ID
print(f"Venturi throat diameter:         {(d_star_m * 39.3701):.3f} in\n")
