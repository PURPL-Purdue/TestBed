# This code is meant to size the orifices for both oxidizer and fuel
from pyfluids import Fluid, FluidsList, Input
from CoolProp.CoolProp import PropsSI
import math



# Define basic parameters
mdot_fuel = 0.0057                          # kg/s
mdot_oxid = 0.0228                          # kg/s
Cd = 0.9
chamber_pressure = 932860.66                # Pa
inlet_pressure = 4.137e6                    # Pa
input_temp = 25 + 273.15                    # K
m2_to_mm2 = 1e6                             # Multiply by this value

# Instantiate Fluid objects for oxid and fuel
oxygen = Fluid(FluidsList.Oxygen).with_state(Input.temperature(25), Input.pressure(inlet_pressure))
methane = Fluid(FluidsList.Methane).with_state(Input.temperature(25), Input.pressure(inlet_pressure))

# Determine density of the oxid and fuel
rho_o2 = oxygen.density                 # kg/m^3
rho_ch4 = methane.density               # kg/m^3

# Calculate specific heat ratio values for both at given temp and pressure
gamma_o2 = oxygen.specific_heat / PropsSI("Cvmass", "T", input_temp, "P", inlet_pressure, "Oxygen")
gamma_ch4 = methane.specific_heat / PropsSI("Cvmass", "T", input_temp, "P", inlet_pressure, "Methane")

# Calculate upstream total pressure values for choking
upstream_O2 = chamber_pressure / math.pow(2 / (gamma_o2 + 1), gamma_o2 / (gamma_o2 - 1))
upstream_CH4 = chamber_pressure / math.pow(2 / (gamma_ch4 + 1), gamma_ch4 / (gamma_ch4 - 1))

# Solve for orifice areas and diameters
gamma_exp_o2 = (gamma_o2 + 1) / (gamma_o2 - 1)
gamma_exp_ch4 = (gamma_ch4 + 1) / (gamma_ch4 - 1)
oxid_area = (mdot_oxid / (Cd * math.sqrt(gamma_o2 * rho_o2 * inlet_pressure * math.pow(2 / (gamma_o2 + 1), gamma_exp_o2)))) * m2_to_mm2
fuel_area = (mdot_fuel / (Cd * math.sqrt(gamma_ch4 * rho_ch4 * inlet_pressure * math.pow(2 / (gamma_ch4 + 1), gamma_exp_ch4)))) * m2_to_mm2

oxid_diam = 2 * math.sqrt(oxid_area / math.pi)                      # mm
fuel_diam = 2 * math.sqrt(fuel_area / math.pi)                      # mm


# Print out the results
print(f"O2:  gamma = {gamma_o2:.4f}; area = {oxid_area:.4e} mm^2; diameter = {oxid_diam:.3f} mm")
print(f"CH4: gamma = {gamma_ch4:.4f}; area = {fuel_area:.4e} mm^2; diameter = {fuel_diam:.3f} mm")


