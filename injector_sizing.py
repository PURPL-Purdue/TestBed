# Maelstrom Injector Sizing Code
# Authors: Dominik Sloup
# First Created: 06/04/2025
# Last Updated: 06/04/2025
# Calculations done in SI units

import os 
import numpy as np
import math
import matplotlib.pyplot as plt
import CEA_Wrap as CEA
from matplotlib.widgets import TextBox
from pyfluids import Fluid, FluidsList, Input

def clc():
    os.system('cls' if os.name == 'nt' else 'clear')

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

fuel = CEA.Fuel("Jet-A(L)")
oxidizer = CEA.Oxidizer("O2")
gamma = 1.4 # GOx specific heat ratio
C_d_ox = 0.8 # Oxidizer orifice anticipated C_d (N/A)
C_d_fu = 0.8 # Fuel orifice anticipated C_d (N/A)
F = 500 # Desired thrust of the engine (lbf)
p_c = 500 # Optimal chamber pressure (psi)
OF = 1.2 # Nominal OF Ratio
K_fu = 0.35 # Desired Fuel injector stiffness (N/A)
K_ox = 1 # Desired Ox injector stiffness (N/A)
g = 9.81 # Gravitational constant (m/s^2)
ox_orifice_num = 12 # Number of oxidizer orifices
fuel_orifice_num = 2 * ox_orifice_num # Number of fuel orifices

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

p_c_pa = p_c * 6894.76
F_N = F * 4.44822 # Desired thrust (N)
# m_dot_kg = m_dot * 0.453592 Desired mass flow (kg/s)
# A_f_m2 = A_f * 0.00064516 # Combined fuel injection area (m^2)

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

clc()
# Specific impulse of engine at ambient pressure
Isp = CEA.RocketProblem(materials=[fuel,oxidizer], pressure=p_c, o_f=OF, analysis_type="frozen").run_cea().isp

m_dot = F_N / (Isp * g) # Total combined mass flow (kg/s)
m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

dp_fu = K_fu * p_c # Desired pressure drop across injector (psi)
dp_ox = K_ox * p_c

p_m_fu = p_c + dp_fu # Manifold pressure (psi)
p_m_ox = p_c + dp_ox # Ox manifold pressure (psi)

p_m_ox_pa = p_m_ox * 6894.76 # Manifold pressure (Pa)

dp_fu_pa = dp_fu * 6894.76 # Pressure drop (Pa)

oxygen = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_m_ox_pa), Input.temperature(20))
rho_ox = oxygen.density # Oxygen density at manifold pressure and ambient temperature (kg/m^3)

rho_fu = 810 # RP-1 density at ambient temp and pressure (may want to update this to match manifold press), (kg/m^3)

# Compressible flow equation
A_fu = m_dot_fu / (C_d_fu * np.sqrt(2 * rho_fu * dp_fu_pa)) # Total fuel injection area (m^2)
d_fu = 2 * np.sqrt((A_fu / fuel_orifice_num) / np.pi) # Fuel orifice diameter (m)
d_fu_in = d_fu * 39.3701 # Fuel orifice diameter (in)

# Incompressible flow equation (Choked condition)
A_ox = m_dot_ox / (C_d_ox * np.sqrt(2 * rho_ox * p_m_ox_pa * (gamma / (gamma - 1))
* ((p_c_pa / p_m_ox_pa) ** (2 / gamma) - (p_c_pa / p_m_ox_pa) ** ((gamma + 1) / gamma)))) # Total ox injection area (m^2)
d_ox = 2 * np.sqrt((A_ox / ox_orifice_num) / np.pi) # Ox orifice diameter (m)
d_ox_in = d_ox * 39.3701 # Ox orifice diameter (in)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

print(f"Specific impulse: {math.floor(Isp)} s")
print(f"Total mass flow: {m_dot:.3f} kg/s")
print(f"Fuel mass flow: {m_dot_fu:.3f} kg/s")
print(f"Ox mass flow: {m_dot_ox:.3f} kg/s")

print(f"\nTotal fuel injection area: {A_fu:.6f} m^2")
print(f"Fuel orifice diameter: {d_fu_in:.3f} in")
print(f"Fuel feed pressure: {p_m_fu:.2f} psi")
print(f"\nTotal oxidizer injection area: {A_ox:.6f} m^2")
print(f"Oxidizer orifice diameter: {d_ox_in:.3f} in")
print(f"Oxidizer feed pressure: {p_m_ox:.2f} psi")
print("")