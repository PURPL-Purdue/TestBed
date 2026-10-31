# Pintle Fill Time Analysis Code
# Authors: Chanak Gautam
# First Created: 10/29/2026
# Last Updated: 10/29/2026
# Calculations done in SI unit

import math
import pint
from pyfluids import Fluid, FluidsList, Input
from CoolProp.CoolProp import PropsSI

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


# ---------------------------------------------------------------------------
# Inputs (imperial). All pressures are absolute.
# ---------------------------------------------------------------------------
LIQUID = FluidsList.nDodecane               # propellant (n-dodecane is the usual RP-1 surrogate)
GAS = FluidsList.Nitrogen                   # purge gas

V0_imp = Q_(2.0, "in**3")                   # manifold internal volume
P_f_imp = Q_(500.0, "psi")                  # feed pressure upstream of liquid inlet orifice
P_a_imp = Q_(14.7, "psi")                   # pressure downstream of gas outlet
P_amb_imp = Q_(14.7, "psi")                 # initial manifold pressure (purge gas at rest)
T_gas_imp = Q_(70.0, "degF")                # purge gas temperature (held constant, isothermal)
T_liq_imp = Q_(70.0, "degF")                # propellant temperature
D_in_imp = Q_(0.0625, "in")                 # liquid inlet orifice diameter
D_out_imp = Q_(0.0625, "in")                # gas outlet orifice diameter

Cd_liq = 0.61                               # estimate: sharp-edged orifice, liquid
Cd_gas = 0.80                               # estimate: gas through valve/fitting
dt = 1e-3                                   # s


# ---------------------------------------------------------------------------
# Convert inputs to SI
# ---------------------------------------------------------------------------
V0 = V0_imp.m_as("m**3")                    # m^3
P_f = P_f_imp.m_as("Pa")                    # Pa
P_a = P_a_imp.m_as("Pa")                    # Pa
P_amb = P_amb_imp.m_as("Pa")                # Pa
T = T_gas_imp.m_as("K")                     # K
T_liq = T_liq_imp.m_as("K")                 # K
A_in = math.pi / 4 * D_in_imp.m_as("m") ** 2    # m^2
A_out = math.pi / 4 * D_out_imp.m_as("m") ** 2  # m^2

if P_f <= P_amb:
    raise ValueError("Feed pressure must exceed initial manifold pressure for liquid to flow in")


# ---------------------------------------------------------------------------
# Fluid properties (pyfluids takes temperature in degC, pressure in Pa)
# ---------------------------------------------------------------------------
liquid = Fluid(LIQUID).with_state(Input.pressure(P_f), Input.temperature(T_liq - 273.15))
gas = Fluid(GAS).with_state(Input.pressure(P_amb), Input.temperature(T - 273.15))

rho_L = liquid.density                      # kg/m^3, liquid treated as incompressible

R_u = PropsSI("gas_constant", GAS.coolprop_name)    # J/(mol K), universal gas constant
M_g = PropsSI("molar_mass", GAS.coolprop_name)      # kg/mol
R = R_u / M_g                               # J/(kg K), specific gas constant of purge gas
gamma = gas.specific_heat / PropsSI("Cvmass", "T", T, "P", P_amb, GAS.coolprop_name)


# ---------------------------------------------------------------------------
# Initialize (t = 0, manifold full of purge gas)
# rho_g comes from the real-gas EOS rather than P/RT, so Z carries the real-gas correction
# (the ideal-gas rho_g in the flowchart would make Z = 1 by construction).
# ---------------------------------------------------------------------------
V_g = V0                                    # m^3
V_L = 0.0                                   # m^3
P_m = P_amb                                 # Pa
t_tot = 0.0                                 # s

rho_g = gas.density                         # kg/m^3
m_g = rho_g * V0                            # kg
Z = P_m * V_g / (m_g * R * T)               # compressibility factor


# Print out initial state
print(f"Liquid: {LIQUID.name}, rho_L = {Q_(rho_L, 'kg/m**3').m_as('lb/ft**3'):.3f} lb/ft^3")
print(f"Gas:    {GAS.name}, R = {R:.3f} J/(kg K), gamma = {gamma:.4f}, Z = {Z:.5f}")
print(f"        rho_g = {Q_(rho_g, 'kg/m**3').m_as('lb/ft**3'):.5f} lb/ft^3, "
      f"m_g = {Q_(m_g, 'kg').m_as('lb'):.4e} lb")
print(f"A_in = {Q_(A_in, 'm**2').m_as('in**2'):.4e} in^2, A_out = {Q_(A_out, 'm**2').m_as('in**2'):.4e} in^2")
