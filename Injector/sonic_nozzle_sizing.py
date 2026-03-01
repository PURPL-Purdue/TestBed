# Maelstrom Sonic Nozzle Sizing Code
# Authors: Chanak Gautam
# First Created: 2/17/2026
# Last Updated: 2/28/2026
# Calculations done in SI unit

import numpy as np
from pathlib import Path
from ruamel.yaml import YAML
import yaml
import math

def find_yaml(filename="Maelstrom.yaml", start_dir=None):
    start_dir = Path(start_dir or Path.cwd())
    for path in start_dir.rglob(filename):
        return path
    raise FileNotFoundError(f"{filename} not found starting from {start_dir}")

def p1_from_mach(p1):
    upstream_mach = (mdot / (p1 * A_2)) * math.sqrt(R * t_0 / gamma)
    return p_star * (
        ((gamma + 1) / 2) /
        (1 + (gamma - 1)/2 * upstream_mach**2)
    )**(gamma / (gamma - 1))

def fp(p1):
    return p1 - p1_from_mach(p1)
    
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

mdot = 0.8 * data["gox_design_mdot"] * lb_to_kg # 80% GOx mass flow rate (kg/s)
p_2 = data["gox_feed_pressure"] * psi_to_pa # Dowstream pressure relative to nozzle (Pa)
t_0 = 300 # Ambient temperature (K)
A_2 = (data["gox_tube_inner_dia"] * in_to_m / 2)**2 * np.pi # pipe cross-sectional area (in^2)
gamma = 1.4 # specific heat ratio GOx
R = 8.314 / 0.032 # specific gas constant GOx (J/kg*K)

downstream_mach = (mdot / (p_2 * A_2)) * math.sqrt(R * t_0 / gamma) # downstream mach number

# stagnation pressure from static pressure 
#p_0 = p_2 * (1 + (gamma - 1)/2 * downstream_mach**2)**(gamma/(gamma - 1))

# area ratio from isentropic relation
area_ratio = ((1 / downstream_mach) * ((2 / (gamma + 1)) * (1 + (gamma - 1)/2 * downstream_mach**2))
              **((gamma + 1) / (2 * (gamma - 1))))

# critical area
A_star = A_2 / area_ratio

# calculating critical pressure (Mach number @ throat = 1)
p_star = p_2 * ((1 + (gamma - 1)/2 * downstream_mach**2) * (2 / (gamma + 1)))**(gamma / (gamma - 1))

# ------------------------------------------------------------
# USING SECANT METHOD TO FIND UPSTREAM PRESSURE
# ------------------------------------------------------------

# two starting guesses for upstream pressure
p_prev = 0. * p_2
p_curr = 1.1 * p_2

# calculating the zeroing function
f_prev = fp(p_prev)
f_curr = fp(p_curr)

iterations = 0

for i in range(80):
    diff = (f_curr - f_prev)
    p_next = p_curr - f_curr * (p_curr - p_prev) / diff

    # checking if converged
    if abs(p_next - p_curr) / abs(p_next) < 0.00000001:
        p_1 = p_next
        break
    
    iterations = iterations + 1
    print(p_next, p_curr)

    # setting previous set of pressure guesses equivalent to new calculated values
    p_prev, f_prev = p_curr, f_curr
    p_curr = p_next
    f_curr = fp(p_curr)

print(f"=======SONIC NOZLE SIZING=======")
print(f"Downstream Mach Number: {downstream_mach:.2f}")
print(f"Sonic Nozzle Critical Area: {A_star * 1550:.3f} in^2")
print(f"Sonic Nozzle Critical Diameter: {math.sqrt(4 * A_star / np.pi) / in_to_m:.3f} in")
print(f"Sonic Nozzle Critical Pressure: {p_star / psi_to_pa:.2f} psi")
print(f"Required Upstream Static Pressure: {p_1 / psi_to_pa:.2f} psi (Iterations: {iterations})")
print(f"Upstream Mach Number: {(mdot / (p_1 * A_2)) * math.sqrt(R * t_0 / gamma):.2f}")