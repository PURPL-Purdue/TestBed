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

F = 1.6 # Nominal thrust of the torch (lbf)
p_c = 140 # Nominal chamber pressure (psi)
OF = 3 # Nominal OF Ratio

fuel_choice = 'kerosene'
# =fuel_choice = 'methane'
# fuel_choice = 'hydrogen'

p_ox = 600 # Oxidizer feed pressure (psi)
ox_orifice_num = 5 # Number of oxidizer orifices (N/A)
Cd_ox = 0.8 # Oxidizer orifice anticipated C_d (N/A)

p_fu = 740 # Fuel feed pressure (psi)
fu_orifice_num = 1 # Number of fuel orifices (N/A)
Cd_fu = 0.8 # Fuel orifice anticipated C_d (N/A)

cstar_eff = 0.7 # C star efficiency factor (N/A)

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

g0 = 32.174 # Gravitational constant (ft/s^2)
gamma_o2 = 1.4 # Oxygen Specific heat ratio
gamma_h2 = 1.41 # Hydrogen Specific heat ratio
gamma_ch4 = 1.3 # Methane Specific heat ratio
R_o2 = 259.8 # Oxygen Specific gas constant
R_h2 = 4124 # Hydrogen Specific gas constant
R_ch4 = 518.3 # Methane Specific gas constant
T = 293.15 # Ambient temperature (K)

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

lbf_to_N = 4.4482
psi_to_pa = 6894.76
ft_to_m = 0.3048
lb_to_kg = 0.453592
m_to_in = 39.3701

p_ox_pa = p_ox * psi_to_pa # Oxidizer feed pressure (Pa)
p_fu_pa = p_fu * psi_to_pa # Fuel feed pressure (Pa)
p_c_pa = p_c * psi_to_pa # Chamber pressure (Pa)
F_N = F * lbf_to_N # Desired thrust (N)

# ──────────────────────────────────────────────────────────────
#  FUNCTION DEFINITIONS
# ──────────────────────────────────────────────────────────────

def main(fuel_choice):

    # Oxygen density at manifold pressure and ambient temperature (kg/m^3)
    rho_ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_ox_pa), Input.temperature(T-273.15)).density 

    if fuel_choice == 'kerosene':
        engine = CEA_Obj(oxName='GOX', fuelName='JetA')
        rho_fu = 810 # Kerosene density (kg/m^3) @ STP
        m_dot_fu, m_dot_ox, cstar_real = get_mdots(engine,F,p_c,OF,cstar_eff)
        A_fu = get_liquid_inj_area(m_dot_fu, p_fu_pa, p_c_pa, rho_fu, Cd_fu)

    elif fuel_choice == 'methane':
        engine = CEA_Obj(oxName='GOX', fuelName='CH4(g)')
        rho_fu = Fluid(FluidsList.Methane).with_state(Input.pressure(p_fu_pa), Input.temperature(T-273.15)).density
        m_dot_fu, m_dot_ox, cstar_real = get_mdots(engine,F,p_c,OF,cstar_eff)
        A_fu = get_gas_inj_area(m_dot_fu, p_fu_pa, p_c_pa, rho_fu, Cd_fu, gamma_ch4, R_ch4, T)

    elif fuel_choice == 'hydrogen':
        engine = CEA_Obj(oxName='GOX', fuelName='H2(g)')
        rho_fu = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(p_fu_pa), Input.temperature(T-273.15)).density
        m_dot_fu, m_dot_ox, cstar_real = get_mdots(engine,F,p_c,OF,cstar_eff)
        A_fu = get_gas_inj_area(m_dot_fu, p_fu_pa, p_c_pa, rho_fu, Cd_fu, gamma_h2, R_h2, T)

    A_ox = get_gas_inj_area(m_dot_ox, p_ox_pa, p_c_pa, rho_ox, Cd_ox, gamma_o2, R_o2, T) 
    A_t = (m_dot_fu + m_dot_ox) * (cstar_real * ft_to_m) / p_c_pa

    d_fu = 2 * np.sqrt((A_fu / fu_orifice_num) / np.pi) # Fuel orifice diameter (m)
    d_ox = 2 * np.sqrt((A_ox / ox_orifice_num) / np.pi) # Ox orifice diameter (m)
    d_t = 2 * np.sqrt(A_t / np.pi) # Throat diameter (m)

    return cstar_real, m_dot_fu, m_dot_ox, A_fu, d_fu, A_ox, d_ox, d_t

def get_mdots(engine,F,p_c,OF,cstar_eff):

    cstar = engine.get_Cstar(Pc = p_c, MR = OF)
    cstar_real = cstar * cstar_eff
    Cf, IspVac, IspSL = engine.get_PambCf(Pamb=14.7, Pc=p_c, MR=OF, eps=1)
    Isp = (cstar_real * Cf / g0)
    m_dot_lb = F / Isp # Total combined mass flow (lb/s)
    m_dot = m_dot_lb * lb_to_kg # Total combined mass flow (kg/s)
    m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
    m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

    return m_dot_fu, m_dot_ox, cstar_real

def get_gas_inj_area(mdot, p_feed, p_c, rho, Cd, gamma, R, T0):

    p_cr = (2 / (gamma + 1)) ** (gamma / (gamma -1))

    if p_feed / p_c > p_cr: # Choked condition
        area = mdot / (Cd * np.sqrt(rho * p_feed * (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))))
    
    else: # Unchoked condition
        area = mdot / (Cd * p_c * np.sqrt((gamma/(R*T0)) * (2/(gamma-1)) * ((p_feed/p_c)**(2/gamma) - (p_feed/p_c)**((gamma+1)/gamma))))

    return area

def get_liquid_inj_area(mdot, p_feed, p_c, rho, Cd):

    return mdot / (Cd * np.sqrt(2 * rho * (p_feed - p_c))) # Total fuel injection area (m^2)

# ──────────────────────────────────────────────────────────────
#  RESULTS
# ──────────────────────────────────────────────────────────────

cstar_real, m_dot_fu, m_dot_ox, A_fu, d_fu, A_ox, d_ox, d_t = main(fuel_choice)

print(f"Expected characteristic velocity: {math.floor(cstar_real * ft_to_m)} m/s")
print(f"Total mass flow: {m_dot_fu + m_dot_ox:.4f} kg/s")
print(f"Fuel mass flow: {m_dot_fu:.4f} kg/s")
print(f"Ox mass flow: {m_dot_ox:.4f} kg/s")

print(f"\nThroat diameter: {d_t * m_to_in:.3f} in")

print(f"\nTotal fuel injection area: {A_fu:.6f} m^2")
print(f"Fuel orifice diameter: {d_fu * m_to_in:.3f} in")

print(f"\nTotal oxidizer injection area: {A_ox:.6f} m^2")
print(f"Oxidizer orifice diameter: {d_ox * m_to_in:.3f} in")