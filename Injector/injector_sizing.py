# Maelstrom Injector Sizing Code
# Authors: Dominik Sloup
# First Created: 06/06/2025
# Last Updated: 06/06/2025
# Calculations done in SI units

import os 
import numpy as np
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
import CEA_Wrap as CEA
from matplotlib.widgets import TextBox
from pyfluids import Fluid, FluidsList, Input
from pathlib import Path
from ruamel.yaml import YAML
import yaml

def clc():
    os.system('cls' if os.name == 'nt' else 'clear')

def find_yaml(filename="params.yaml", start_dir=None):
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
# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

gamma = 1.4 # GOx specific heat ratio
C_d_ox = 0.8 # Oxidizer orifice anticipated C_d (N/A)
C_d_fu = 0.8 # Fuel orifice anticipated C_d (N/A)
F = 500 # Desired thrust of the engine (lbf)
p_c = 500 # Optimal chamber pressure (psi)
OF = 1 # Nominal OF Ratio
K_fu = 0.5 # Desired Fuel injector stiffness (N/A)
K_ox = 1.2 # Desired Ox injector stiffness (N/A)
g = 9.81 # Gravitational constant (m/s^2)
g0 = 32.174 # Gravitational constant (ft/s^2)
d_line_fu = 0.4 # Fuel Line diameter (in)
d_line_ox = 0.65 # Oxidizer line diameter (in)
fu_orifice_num = 12 # Number of fuel orifices
ox_orifice_num = 24 # Number of oxidizer orifices
Isp_eff = 0.9 # Mixing efficiency factor (N/A)

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

lbf_to_N = 4.44822
psi_to_pa = 6894.76
ft_to_m = 0.3048
lb_to_kg = 0.453592
m_to_in = 39.3701

p_c_pa = p_c * psi_to_pa # Chamber pressure (Pa)
F_N = F * lbf_to_N # Desired thrust (N)

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

clc()
# Specific impulse of engine at ambient pressure

engine = CEA_Obj(oxName='GOX', fuelName='JetA')

cstar = engine.get_Cstar(Pc = p_c, MR = OF)
Cf, IspVac, IspSL = engine.get_PambCf(Pamb=14.7, Pc=p_c, MR=OF, eps=4.43)
Isp = (cstar * Cf / g0)

Isp_real = Isp * Isp_eff

# print(engine.get_full_cea_output(Pc = p_c, MR = OF, eps = 4.43))

cstar_m = cstar * ft_to_m

m_dot_lb = F / Isp_real # Total combined mass flow (lb/s)
m_dot = m_dot_lb * lb_to_kg # Total combined mass flow (kg/s)

m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

dp_fu = K_fu * p_c # Desired pressure drop across fuel injector (psi)
dp_ox = K_ox * p_c # Desired pressure drop across ox injector (psi)

p_m_fu = p_c + dp_fu # Manifold pressure (psi)
p_m_ox = p_c + dp_ox # Ox manifold pressure (psi)

p_m_ox_pa = p_m_ox * 6894.76 # Manifold pressure (Pa)

dp_fu_pa = dp_fu * 6894.76 # Pressure drop (Pa)

oxygen = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_m_ox_pa), Input.temperature(20))
rho_ox = oxygen.density # Oxygen density at manifold pressure and ambient temperature (kg/m^3)

rho_fu = 810 # RP-1 density at ambient temp and pressure (may want to update this to match manifold press), (kg/m^3)

# Compressible flow equation
A_fu = m_dot_fu / (C_d_fu * np.sqrt(2 * rho_fu * dp_fu_pa)) # Total fuel injection area (m^2)
d_fu = 2 * np.sqrt((A_fu / fu_orifice_num) / np.pi) # Fuel orifice diameter (m)
d_fu_in = d_fu * 39.3701 # Fuel orifice diameter (in)

# Incompressible flow equation (Choked condition)

A_ox = m_dot_ox / (C_d_ox * np.sqrt(rho_ox * p_m_ox_pa * (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))))
A_t = m_dot * cstar_m / p_c_pa

#  Incompressible flow equation (Unchoked condition)
# A_ox = m_dot_ox / (C_d_ox * np.sqrt(2 * rho_ox * p_m_ox_pa * (gamma / (gamma - 1))
# * ((p_c_pa / p_m_ox_pa) ** (2 / gamma) - (p_c_pa / p_m_ox_pa) ** ((gamma + 1) / gamma)))) # Total ox injection area (m^2)
d_ox = 2 * np.sqrt((A_ox / ox_orifice_num) / np.pi) # Ox orifice diameter (m)
d_ox_min_in = 2 * np.sqrt(A_ox / np.pi) * 39.3701 # Minimum ox inlet diameter (in)

d_t = 2 * np.sqrt(A_t / np.pi)

d_t_in = d_t * m_to_in
d_ox_in = d_ox * m_to_in # Ox orifice diameter (in)

d_line_fu_m = d_line_fu * 0.0254
d_line_ox_m = d_line_ox * 0.0254

A_line_fu = np.pi * (d_line_fu_m/2) ** 2 
A_line_ox = np.pi * (d_line_ox_m/2) ** 2 

# Line velocity calculations:
v_fu = m_dot_fu / (A_line_fu * rho_fu)
v_ox = m_dot_ox / (A_line_ox * rho_ox)

v_ox_ft = v_ox / ft_to_m
v_fu_ft = v_fu / ft_to_m

rho_ox_amb = Fluid(FluidsList.Oxygen).with_state(Input.pressure(101325), Input.temperature(20)).density
rho_n2_amb = Fluid(FluidsList.Nitrogen).with_state(Input.pressure(101325), Input.temperature(20)).density

V_dot_ox = m_dot_ox / rho_ox_amb # Oxygen volumetric flow at STP (m^3/s)
V_dot_n2 = m_dot_fu / rho_n2_amb # Nitrogen tank press volumetric flow at STP (m^3/s)

SCFM_ox =  V_dot_ox / np.pow(ft_to_m, 3) * 60
SCFM_n2 =  V_dot_n2 / np.pow(ft_to_m, 3) * 60

data["of_ratio"] = OF
data["gox_design_mdot"] = float(np.round(m_dot_ox / lb_to_kg, 3))
data["rp_design_mdot"] = float(np.round(m_dot_fu / lb_to_kg, 3))
data["gox_feed_pressure"] = p_m_ox
data["rp_feed_pressure"] = p_m_fu
data["gox_stiffness"] = K_ox
data["rp_stiffness"] = K_fu
data["gox_line_velocity"] = float(np.round(v_ox_ft, 2))
data["rp_line_velocity"] = float(np.round(v_fu_ft, 2))
data["gox_tube_inner_dia"] = d_line_ox
data["rp_tube_inner_dia"] = d_line_fu
data["rp_tube_inner_dia"] = d_line_fu
data["gox_SCFM"] = float(np.round(SCFM_ox, 2))
data["n2_SCFM"] = float(np.round(SCFM_n2, 2))

with open(yaml_path, "w") as f:
    yaml.dump(data, f)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

print(f"Specific impulse ({Isp_eff} efficiency): {math.floor(Isp_real)} s")
print(f"Total mass flow: {m_dot:.3f} kg/s")
print(f"Fuel mass flow: {m_dot_fu:.3f} kg/s")
print(f"Ox mass flow: {m_dot_ox:.3f} kg/s")

print(f"\nThroat diameter: {d_t_in:.2f} in")
print(f"\nTotal fuel injection area: {A_fu:.6f} m^2")
print(f"Fuel orifice diameter: {d_fu_in:.3f} in")
print(f"Fuel feed pressure: {p_m_fu:.2f} psi")
print(f"\nTotal oxidizer injection area: {A_ox:.6f} m^2")
print(f"Oxidizer orifice diameter: {d_ox_in:.3f} in")
print(f"Oxidizer minimum inlet diameter: {d_ox_min_in:.3f} in")
print(f"Oxidizer feed pressure: {p_m_ox:.2f} psi")
print("")
print(f"Fuel line velocity: {v_fu_ft:.3f} ft/s")
print(f"Ox line velocity: {v_ox_ft:.3f} ft/s")