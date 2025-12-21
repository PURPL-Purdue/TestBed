# Maelstrom Venturi Sizing Code
# Authors: Dominik Sloup
# First Created: 12/20/2025
# Last Updated: 12/20/2025
# Calculations done in SI units

import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
from pathlib import Path
from ruamel.yaml import YAML
import yaml

def find_yaml(filename="Maelstrom.yaml", start_dir=None):
    start_dir = Path(start_dir or Path.cwd())
    for path in start_dir.rglob(filename):
        return path
    raise FileNotFoundError(f"{filename} not found starting from {start_dir}")

yaml_path = find_yaml()

    # Load YAML with formatting preserved
yaml = YAML()
yaml.preserve_quotes = True

with open(yaml_path, "r") as f:
    data = yaml.load(f)

# ------------------------------------------------------------
# UNIT CONSTANTS
# ------------------------------------------------------------
psi_to_pa = 6894.76
in_to_m = 0.0254
lb_to_kg = 0.453592

# ------------------------------------------------------------
# FLUID + GEOMETRY PARAMETERS
# ------------------------------------------------------------

d_u = data["fuel_tube_inner_dia"] # Oxidizer line diameter (in)
p_f = (data["fuel_injector_pressure"] + data["regen_pressure_drop"]) * psi_to_pa # Feed pressure (Pa)
mdot = data["fuel_design_mdot"] * lb_to_kg # Fuel mass flow rate (kg/s)

d_th_guess = 0.122 # Guess for venturi throat area (in) to estimate beta ratio
C_d = 0.931     # Venturi discharge coefficient (Tested on prototype, to be corrected)
C_p = 0.8 # Pressure recovery factor (literature source)
beta = d_th_guess / d_u

rho_fu = Fluid(FluidsList.Ethanol).with_state(Input.pressure(p_f), Input.temperature(20)).density

P_cr = (C_p * C_d ** 2) / (1 - beta ** 4) # Critical pressure ratio (above this value, venturi uncavitates)
p_u = p_f / (P_cr - 0.1) # Upstream venturi pressure (Pa)

# Set throat pressure to be equal to saturation pressure of ethanol (Pa)
p_th = Fluid(FluidsList.Ethanol).with_state(Input.quality(0.0),Input.temperature(20)).pressure

d_u = d_u * in_to_m # Convert the line diameter to meters

K = C_d * (np.pi / 4) * np.sqrt(2 * rho_fu * (p_u - p_th))
d_th = np.pow(mdot ** 2 / (K ** 2 + mdot ** 2 / d_u ** 4), 0.25) 

d_th = d_th / in_to_m # Convert throat diameter to inches

print(f"Venturi upstream pressure: {p_u / psi_to_pa:.0f} psi")
print(f"Venturi dpwnstream pressure: {p_f / psi_to_pa:.0f} psi")
print(f"Critical pressure ratio: {P_cr:.2f}")
print(f"Actual pressure ratio: {p_f / p_u:.2f}")
print(f"Ethanol saturation pressure: {p_th / psi_to_pa:.2f} psi")
print(f"Cavitating venturi throat diameter: {d_th:.3f} in")

data["fuel_feed_pressure"] = int(p_u / psi_to_pa)
data["critical_pressure_ratio"] = float(round(P_cr,2))
data["minumum_pressure_ratio"] = float(round(p_f / p_u,2))
data["venturi_throat_diameter"] = float(round(d_th,3))

with open(yaml_path, "w") as f:
    yaml.dump(data, f)