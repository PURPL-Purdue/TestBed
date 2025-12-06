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
from matplotlib.widgets import TextBox
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

# ──────────────────────────────────────────────────────────────
#  INPUTS
# ──────────────────────────────────────────────────────────────

F = data["thrust"] # Nominal thrust of the torch (lbf)
p_c = data["chamber_pressure"] # Nominal chamber pressure (psi)
p_e = data["exhaust_pressure"] # Exhaust pressure (psi)
dp_regen = data["regen_pressure_drop"] # Pressure drop across regenerative nozzle 
OF = data["of_ratio"] # Nominal OF Ratio
T = 293.15 # Fuel and oxidizer temperature (K)

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

K_fu = data["rp_stiffness"] # Desired Ox injector stiffness (N/A)
fu_orifice_num = data["rp_orifice_number"] # Number of fuel orifices (N/A)
Cd_fu = data["rp_discharge_coeff"] # Fuel orifice anticipated C_d (N/A)
d_line_fu = data["rp_tube_inner_dia"] # Oxidizer line diameter (in)

cstar_eff = data["cstar_eff"] # C star efficiency factor (N/A)
Isp_eff = data["isp_eff"] # Nozzle efficiency factor (N/A)
Lstar = data["Lstar"] # Characteristic length (in)

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

gamma_o2 = 1.4 # GOx specific heat ratio
g = 9.81 # Gravitational constant (m/s^2)
g0 = 32.174 # Gravitational constant (ft/s^2)
rho_fu = 810 # RP-1 density at ambient temp and pressure (may want to update this to match manifold press), (kg/m^3)

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

    engine = CEA_Obj(oxName='GOX', fuelName='JetA')

    p_c_pa = p_c * psi_to_pa # Chamber pressure (Pa)
    p_ox_pa = p_c * (1 + K_ox) * psi_to_pa # Oxidizer injector feed pressure (Pa)
    p_fu_pa = p_c * (1 + K_fu) * psi_to_pa # Fuel injector feed pressure (Pa)

    # Oxygen density at manifold pressure and ambient temperature (kg/m^3)
    rho_ox = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_ox_pa), Input.temperature(T-273.15)).density 
    m_dot_fu, m_dot_ox, cstar_real, Isp = get_mdots(engine,F,p_c,OF,cstar_eff)
    A_fu = get_liquid_inj_area(m_dot_fu, p_fu_pa, p_c_pa, rho_fu, Cd_fu)
    A_ox = get_gas_inj_area(m_dot_ox, p_ox_pa, p_c_pa, rho_ox, Cd_ox, gamma_o2) 
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

    total_length = (V_chamber/(np.pi * cR**2) # This is some ChatGPT Voodoo magic bro I could not derive this math
      + 2*cR/(3*np.tan(cA))
      + eR/np.tan(dA)
      - tF/np.tan(dA)
      - tF/np.tan(cA)
      + tF/np.sin(dA)
      + tF/np.sin(cA)
      - tR/np.tan(dA)
      - tR/np.tan(cA))

    # Chamber sizing process

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

    V_dot_ox = m_dot_ox / rho_ox_amb # Oxygen volumetric flow rate at STP (m^3/s)

    SCFM_ox =  V_dot_ox / np.pow(ft_to_m, 3) * 60
    SCFM_n2 = (m_dot_fu / rho_fu) * (rho_fu / 1.165) * 35.3147 * 60

    data["cstar"] = int(cstar_real * ft_to_m)
    data["isp"] = int(Isp)

    data["throat_diameter"] = float(np.round(d_t * m_to_in, 3))
    data["exit_diameter"] = float(np.round(d_e * m_to_in, 3))
    data["chamber_volume"] = float(np.round(V_chamber, 1))
    data["total_length"] = float(np.round(total_length, 1))
 
    data["rp_orifice_dia"] = float(np.round(d_fu * m_to_in, 3))
    data["gox_orifice_dia"] = float(np.round(d_ox * m_to_in, 3))
    data["gox_design_mdot"] = float(np.round(m_dot_ox / lb_to_kg, 3))
    data["rp_design_mdot"] = float(np.round(m_dot_fu / lb_to_kg, 3))
    data["gox_feed_pressure"] = int(p_ox_pa / psi_to_pa)
    data["rp_feed_pressure"] = int(p_fu_pa / psi_to_pa + dp_regen)
    data["gox_line_velocity"] = float(np.round(v_ox_ft, 2))
    data["rp_line_velocity"] = float(np.round(v_fu_ft, 2))
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

def get_gas_inj_area(mdot, p_feed, p_c, rho, Cd, gamma):

    p_cr = (2 / (gamma + 1)) ** (gamma / (gamma -1))
    if p_c / p_feed < p_cr: # Choked condition
        area = mdot / (Cd * np.sqrt(gamma * rho * p_feed * (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))))
    else: # Unchoked condition
        area = mdot / (Cd * np.sqrt(2 * rho * p_feed * (gamma / (gamma-1)) * (((p_c/p_feed) ** (2/gamma)) - ((p_c/p_feed) ** ((gamma + 1)/gamma)))))

    return area

def get_liquid_inj_area(mdot, p_feed, p_c, rho, Cd):

    return mdot / (Cd * np.sqrt(2 * rho * (p_feed - p_c))) # Total fuel injection area (m^2)

def clc():
    os.system('cls' if os.name == 'nt' else 'clear')

main()