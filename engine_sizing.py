# Hollistic Engine Sizing Code
# Authors: Dominik Sloup
# First Created: 06/06/2025
# Last Updated: 12/20/2025
# Calculations done in SI units

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from matplotlib.widgets import TextBox
from pyfluids import Fluid, FluidsList, Input
from pathlib import Path
from ruamel.yaml import YAML
from Injector.get_area_from_mdot import get_area_from_mdot
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

# ──────────────────────────────────────────────────────────────
#  INPUTS
# ──────────────────────────────────────────────────────────────

F = data["thrust"] # Nominal thrust of the torch (lbf)
p_c = data["chamber_pressure"] # Nominal chamber pressure (psi)
p_e = data["exhaust_pressure"] # Exhaust pressure (psi)
OF = data["of_ratio"] # Nominal OF Ratio
T_gox = data["gox_temperature"] # Oxidizer temperature (K)
T_fuel = data["fuel_temperature"] # Fuel temperature (K). WARNING CEA DOES NOT USE THIS AS OF RN

# If I size the engine around the maximum OF ratio and then move down I will need to decrease 
# gox pressure and increase fuel pressure to stay at the same chamber pressure, thus gox stiffness
# decreases (gox velocity goes up? because Im also decreasing mass flow so its counteracting that.
# Our gox stiffness is already very high, so this likely isn't an issue. Is there anything else I
# am forgetting? Oh yeah, is an increase in fuel mass flow good for the cavitating venturi? I think 
# so. An increase in mass flow means a smaller pressure ratio.

chamber_dia = data['chamber_diameter'] # Chamber diameter at injector (in)
cA = data['converging_angle']
dA = data['diverging_angle']
cR = chamber_dia / 2
tF = data['throat_fillet']
tR = data['throat_diameter'] / 2

K_ox = data["gox_stiffness"] # Desired Ox injector stiffness (N/A)
ox_orifice_num = data["gox_orifice_number"] # Number of oxidizer orifices (N/A)
Cd_ox = data["gox_discharge_coeff"] # Gox orifice anticipated C_d (N/A)
d_line_ox = data["gox_tube_inner_dia"] # Oxidizer line diameter (in)

K_fu = data["fuel_stiffness"] # Desired Ox injector stiffness (N/A)
fu_orifice_num = data["fuel_orifice_number"] # Number of fuel orifices (N/A)
Cd_fu = data["fuel_discharge_coeff"] # Fuel orifice anticipated C_d (N/A)
d_line_fu = data["fuel_tube_inner_dia"] # Oxidizer line diameter (in)

cstar_eff = data["cstar_eff"] # C star efficiency factor (N/A)
Isp_eff = data["isp_eff"] # Nozzle efficiency factor (N/A)
Lstar = data["Lstar"] # Characteristic length (in)

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

gamma_o2 = 1.4 # GOx specific heat ratio
g = 9.81 # Gravitational constant (m/s^2)
g0 = 32.174 # Gravitational constant (ft/s^2)

# ──────────────────────────────────────────────────────────────
#  UNIT CONVERSIONS
# ──────────────────────────────────────────────────────────────

lbf_to_N = 4.44822
psi_to_pa = 6894.76
ft_to_m = 0.3048
lb_to_kg = 0.453592
m_to_in = 39.3701

# ──────────────────────────────────────────────────────────────
#  FUNCTION DEFINITIONS
# ──────────────────────────────────────────────────────────────

def main():

    engine = CEA_Obj(oxName='GOX', fuelName='C2H5OH')

    p_c_pa = p_c * psi_to_pa # Chamber pressure (Pa)
    p_ox_pa = p_c * (1 + K_ox) * psi_to_pa # Oxidizer injector feed pressure (Pa)
    p_fu_pa = p_c * (1 + K_fu) * psi_to_pa # Fuel injector feed pressure (Pa)

    # Oxygen density at manifold pressure and ambient temperature (kg/m^3)
    rho_ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_ox_pa), Input.temperature(T_gox-273.15)).density 
    rho_fu = Fluid(FluidsList.Ethanol).with_state(Input.pressure(p_fu_pa), Input.temperature(T_fuel-273.15)).density
    m_dot_fu, m_dot_ox, cstar_real, Isp = get_mdots(engine,F,p_c,OF,cstar_eff)
    A_fu = get_area_from_mdot('liquid',m_dot_fu, p_fu_pa, p_c_pa, rho_fu, Cd_fu) # Fuel injection area (m^2)
    A_ox = get_area_from_mdot('gas',m_dot_ox, p_ox_pa, p_c_pa, rho_ox, Cd_ox, gamma_o2) # Oxidizer injection area (m^2)
    A_t = (m_dot_fu + m_dot_ox) * (cstar_real * ft_to_m) / p_c_pa

    molwt, gamma = engine.get_Chamber_MolWt_gamma(Pc=p_c, MR=OF)
    Me = np.sqrt((2/(gamma - 1)) * ((p_c/p_e)**((gamma - 1)/gamma) - 1))
    A_e = A_t * (1/Me) * ((2/(gamma + 1) * (1 + (gamma - 1)/2 * Me**2))**((gamma + 1)/(2*(gamma - 1))))

    d_fu = 2 * np.sqrt((A_fu / fu_orifice_num) / np.pi) # Fuel orifice diameter (m)
    d_ox = 2 * np.sqrt((A_ox / ox_orifice_num) / np.pi) # Ox orifice diameter (m)
    d_t = 2 * np.sqrt(A_t / np.pi) # Throat diameter (m)
    d_e = 2 * np.sqrt(A_e / np.pi) # Exit diameter (m)

    eR = d_e / 2 * m_to_in # Exit radius (in)

    V_chamber = Lstar * A_t * m_to_in ** 2

    total_length = (V_chamber/(np.pi * cR**2) # Calculates the length of the chamber based on the required chamber volume
      + 2*cR/(3*np.tan(cA))
      + eR/np.tan(dA)
      - tF/np.tan(dA)
      - tF/np.tan(cA)
      + tF/np.sin(dA)
      + tF/np.sin(cA)
      - tR/np.tan(dA)
      - tR/np.tan(cA))

    rho_ox_amb = Fluid(FluidsList.Oxygen).with_state(Input.pressure(101325), Input.temperature(20)).density

    V_dot_ox = m_dot_ox / rho_ox_amb # Oxygen volumetric flow rate at STP (m^3/s)

    rho_n2_tank   = Fluid(FluidsList.Nitrogen).with_state(Input.pressure(p_fu_pa), Input.temperature(20)).density
    rho_n2_std = Fluid(FluidsList.Nitrogen).with_state(Input.pressure(101325), Input.temperature(20)).density
    
    V_dot_n2 = m_dot_fu / rho_fu
    m_dot_n2 = rho_n2_tank * V_dot_n2     # kg/s N2 required
    Vdot_n2_std = m_dot_n2 / rho_n2_std

    SCFM_ox =  V_dot_ox / np.pow(ft_to_m, 3) * 60

    SCFM_n2 = Vdot_n2_std / np.pow(ft_to_m, 3) * 60

    CdA_fu = Cd_fu * A_fu * (m_to_in) ** 2
    CdA_ox = Cd_ox * A_ox * (m_to_in) ** 2

    data["cstar"] = int(cstar_real * ft_to_m)
    data["isp"] = int(Isp)

    data["throat_diameter"] = float(np.round(d_t * m_to_in, 3))
    data["throat_area"] = float(np.round(A_t * (m_to_in ** 2), 3))
    data["exit_diameter"] = float(np.round(d_e * m_to_in, 3))
    data["chamber_volume"] = float(np.round(V_chamber, 1))
    data["total_length"] = float(np.round(total_length, 1))
 
    data["fuel_orifice_dia"] = float(np.round(d_fu * m_to_in, 3))
    data["gox_orifice_dia"] = float(np.round(d_ox * m_to_in, 3))
    data["gox_design_mdot"] = float(np.round(m_dot_ox / lb_to_kg, 3))
    data["fuel_design_mdot"] = float(np.round(m_dot_fu / lb_to_kg, 3))
    data["gox_feed_pressure"] = int(p_ox_pa / psi_to_pa)
    data["fuel_CdA"] = float(np.round(CdA_fu, 6))
    data["gox_CdA"] = float(np.round(CdA_ox, 6))
    data["fuel_injector_pressure"] = int(p_fu_pa / psi_to_pa)
    data["gox_SCFM"] = float(np.round(SCFM_ox, 2))
    data["n2_SCFM"] = float(np.round(SCFM_n2, 2))
    
    with open(yaml_path, "w") as f:
        yaml.dump(data, f)

def get_mdots(engine,F,p_c,OF,cstar_eff):

    cstar = engine.get_Cstar(Pc = p_c, MR = OF)
    cstar_real = cstar * cstar_eff
    Cf, IspVac, IspSL = engine.get_PambCf(Pamb=14.7, Pc=p_c, MR=OF, eps=4.43)
    Isp = (cstar_real * Cf / g0) * Isp_eff # Specific impulse (s)
    m_dot = F / Isp * lb_to_kg # Total mass flow (ks/s)
    m_dot_fu = m_dot / (1 + OF) # Fuel mass flow (kg/s)
    m_dot_ox = m_dot - m_dot_fu # Oxidizer mass flow (kg/s)

    return m_dot_fu, m_dot_ox, cstar_real, Isp

main()