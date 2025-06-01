# Maelstrom Venturi Sizing Code
# Authors: Dominik Sloup
# First Created: 06/01/2025
# Last Updated: 06/01/2025
# Calculations done in SI units

import math
import numpy as np
import matplotlib.pyplot as plt
import CEA_Wrap as CEA
from matplotlib.widgets import TextBox
from pyfluids import Fluid, FluidsList, Input

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

# rho = pyfluids.Fluid(FluidsList) # Water density at ambient pressure and temp
C_d = 0.905 # Discharge coefficient guess (based on PSP Liquids' testing) https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/936903267/Successful+Water+Flow+1
g = 9.81 # Gravitational constant

water = Fluid(FluidsList.Water).with_state(Input.pressure(101325), Input.temperature(20))
rho = water.density # Water density (kg/m^3)

P_sat = Fluid(FluidsList.Water).with_state(Input.quality(0.0), Input.temperature(20)).pressure # Water saturation pressure at 20 degrees celsius




