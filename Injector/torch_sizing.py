# Torch Sizing Code
# Authors: Dominik Sloup
# First Created: 06/06/2025
# Last Updated: 06/06/2025
# Calculations done in SI units

import numpy as np
import math
from rocketcea.cea_obj import CEA_Obj
from pyfluids import Fluid, FluidsList, Input

# ──────────────────────────────────────────────────────────────
#  INPUTS
# ──────────────────────────────────────────────────────────────

F = 2.24809 # Desired thrust of the torch (lbf)
p_c = 130 # Optimal chamber pressure (psi)
OF = 1.4 # Nominal OF Ratio

# engine = CEA_Obj(oxName='GOX', fuelName='JetA')
engine = CEA_Obj(oxName='GOX', fuelName='CH4(g)')
#engine = CEA_Obj(oxName='GOX', fuelName='CH4(g)')

p_ox = 400 # Oxidizer feed pressure (psi)
ox_orifice_num = 1 # Number of oxidizer orifices
C_d_ox = 0.8 # Oxidizer orifice anticipated C_d (N/A)

p_fu = 400 # Fuel feed pressure (psi)
fu_orifice_num = 1 # Number of fuel orifices
C_d_fu = 0.8 # Fuel orifice anticipated C_d (N/A)

cstar_eff = 0.7 # C star efficiency factor (N/A)

Lstar = 100 # Characteristic length (in)

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

g = 9.81 # Gravitational constant (m/s^2)
g0 = 32.174 # Gravitational constant (ft/s^2)
gamma_ox = 1.4 # GOx specific heat ratio

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

lbf_to_N = 4.44822
psi_to_pa = 6894.76
ft_to_m = 0.3048
lb_to_kg = 0.453592
m_to_in = 39.3701

# ──────────────────────────────────────────────────────────────
#  CALCULATIONS
# ──────────────────────────────────────────────────────────────

p_c_pa = p_c * psi_to_pa # Chamber pressure (Pa)
F_N = F * lbf_to_N # Desired thrust (N)

cstar = engine.get_Cstar(Pc = p_c, MR = OF)
cstar_real = cstar * cstar_eff

Cf, IspVac, IspSL = engine.get_PambCf(Pamb=14.7, Pc=p_c, MR=OF, eps=4.43)
Isp = (cstar_real * Cf / g0)

m_dot_lb = F / Isp # Total combined mass flow (lb/s)

m_dot = m_dot_lb * lb_to_kg # Total combined mass flow (kg/s)
m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

p_ox_pa = p_ox * 6894.76 # Oxidizer feed pressure (Pa)
dp_fu_pa = (p_fu - p_c) * 6894.76 # Pressure drop (Pa)

rho_ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_ox_pa), Input.temperature(20)).density 
# Oxygen density at manifold pressure and ambient temperature (kg/m^3)

rho_fu = Fluid(FluidsList.Methane).with_state(Input.pressure(p_fu * 6894.76), Input.temperature(20)).density

# rho_fu = 810 # RP-1 density at ambient temp and pressure (may want to update this to match manifold press), (kg/m^3)

# Compressible flow equation
# A_fu = m_dot_fu / (C_d_fu * np.sqrt(2 * rho_fu * dp_fu_pa)) # Total fuel injection area (m^2)
A_fu = m_dot_fu / (C_d_fu * np.sqrt(rho_fu * p_fu * psi_to_pa * (2 / (gamma_ox + 1)) ** ((gamma_ox + 1)/(gamma_ox - 1))))


d_fu = 2 * np.sqrt((A_fu / fu_orifice_num) / np.pi) # Fuel orifice diameter (m)
d_fu_in = d_fu * 39.3701 # Fuel orifice diameter (in)

# Incompressible flow equation (Choked condition)

A_ox = m_dot_ox / (C_d_ox * np.sqrt(rho_ox * p_ox_pa * (2 / (gamma_ox + 1)) ** ((gamma_ox + 1)/(gamma_ox - 1))))
d_ox = 2 * np.sqrt((A_ox / ox_orifice_num) / np.pi) # Ox orifice diameter (m)
d_ox_min_in = 2 * np.sqrt(A_ox / np.pi) * 39.3701 # Minimum ox inlet diameter (in)

# Nozzle sizing process

A_t = m_dot * (cstar_real * ft_to_m) / p_c_pa

d_t = 2 * np.sqrt(A_t / np.pi)

d_t_in = d_t * m_to_in
d_ox_in = d_ox * m_to_in # Ox orifice diameter (in)

V_chamber = Lstar * A_t * m_to_in ** 2

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

print(f"Specific impulse): {math.floor(Isp)} s")
print(f"Expected characteristic velocity: {math.floor(cstar_real * ft_to_m)} m/s")
print(f"Total mass flow: {m_dot:.3f} kg/s")
print(f"Fuel mass flow: {m_dot_fu:.3f} kg/s")
print(f"Ox mass flow: {m_dot_ox:.3f} kg/s")

print(f"\nThroat diameter: {d_t_in:.2f} in")

print(f"\nTotal fuel injection area: {A_fu:.6f} m^2")
print(f"Fuel orifice diameter: {d_fu_in:.3f} in")

print(f"\nTotal oxidizer injection area: {A_ox:.6f} m^2")
print(f"Oxidizer orifice diameter: {d_ox_in:.3f} in")
print(f"Oxidizer minimum inlet diameter: {d_ox_min_in:.3f} in")

print(f"\nChamber Volume: {V_chamber:.3f} in^3")