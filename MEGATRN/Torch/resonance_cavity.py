# This code is meant to size the resonance cavity
from pyfluids import Fluid, FluidsList, Input
from CoolProp.CoolProp import PropsSI
from orifice_sizing import oxid_area, oxid_diam, mdot_oxid
import math



# Define basic parameters
Cd = 0.9                                                  
preig_chamber_pressure = 689476                                         # Pa
inlet_pressure = 4.137e6                                                # Pa (600 psia)
input_temp = 25 + 273.15                                                # K
R_universal = 8.314                                                     # J/(mol*K)
m2_to_mm2 = 1e6                                                         # mm^2/m^2

# Instantiate Fluid objects for oxid and fuel
oxygen = Fluid(FluidsList.Oxygen).with_state(Input.temperature(25), Input.pressure(inlet_pressure))  
methane = Fluid(FluidsList.Methane).with_state(Input.temperature(25), Input.pressure(inlet_pressure))

# Determine density of the oxid and fuel
rho_o2 = oxygen.density                                                 # kg/m^3
rho_ch4 = methane.density                                               # kg/m^3

# Determine specific gas constant for oxygen
oxygen_MW = 0.032                                                       # kg/mol
R_ox = R_universal / oxygen_MW                                          # J/(kg*K)

# Calculate specific heat ratio values for both at given temp and pressure
gamma_o2 = oxygen.specific_heat / PropsSI("Cvmass", "T", input_temp, "P", inlet_pressure, "Oxygen")    # Dimensionless (Cp/Cv)
gamma_ch4 = methane.specific_heat / PropsSI("Cvmass", "T", input_temp, "P", inlet_pressure, "Methane")  # Dimensionless (Cp/Cv)


# Calculate oxidizer exit velocity
rho_o2 = oxygen.density                                                 # kg/m^3
v_ox = mdot_oxid / (rho_o2 * (oxid_area / m2_to_mm2))                   # m/s

# Static pressure at the exit of a choked converging nozzle
critical_pressure_ratio = (2 / (gamma_o2 + 1)) ** (gamma_o2 / (gamma_o2 - 1))   # Dimensionless
p_exit = inlet_pressure * critical_pressure_ratio                       # Pa

M_ox = 1.0                                                  
n = p_exit / preig_chamber_pressure                        

# Calculate oxidizer jet axial dimensions (dimensionless quantities)
# Avduevskii et al - https://www.tandfonline.com/doi/epdf/10.1080/00102202.2011.589874?needAccess=true
A_exit = oxid_area / m2_to_mm2                                          # m^2
d_exit = oxid_diam / 1000                                               # m
Xm_d = ((0.8 + 0.085 * (M_ox - 2.1) ** 2) * M_ox * math.sqrt(n - 0.5))   
Xc_Xm = 1.3 + 0.5 * n ** (1 / 3)                       

Xm = Xm_d * d_exit                                                      # m (Mach disk position)
Xc = (Xm_d * Xc_Xm) * d_exit                                            # m (End of first shock cell)

print(f"Cavity Position Range: {Xm * 1000 :.2f} to {Xc * 1000 :.2f}mm")

# Calculate Nozzle to cavity gap
# Using gap at the Hartmann–Trolle estimated cell length - https://iris.polito.it/retrieve/handle/11583/2782256/e384c433-0327-d4b2-e053-9f05fe0a1d67/resonant_igniter.pdf
s_d = 1.8 * ((p_exit / preig_chamber_pressure) - 1.9) ** 0.5            # Dimensionless (s/d)
cavity_dist = s_d * d_exit                                              # m

print(f"Hartmann-Trolle Estimated Cavity Gap: {cavity_dist * 1000:.2f} mm")

# Cavity Entrance Diameter calculation


# Static exit pressure for a choked orifice
exit_pressure_ratio = p_exit / preig_chamber_pressure       

# Estimated diameter at distance x (1) and max diameter (2)
# Try testing both with CFD, but using (2) for DCR
d_jet1 = d_exit * (1 + 1.45 * (1 - math.exp(-0.322 * Xm / d_exit)) * math.log(exit_pressure_ratio)) * 1000  # mm
d_jet2 = d_exit * (1.45 - (0.65 / exit_pressure_ratio ** 0.25)) * math.sqrt(exit_pressure_ratio) * 1000     # mm

D_d1 = d_jet1 / (d_exit * 1000)                            
D_d2 = d_jet2 / (d_exit * 1000)                            

print(f"Jet diameter1: {d_jet1:.3f} mm")
print(f"D/d1 Ratio: {D_d1:.3f}")

print(f"Jet diameter2: {d_jet2:.3f} mm")
print(f"D/d2 Ratio: {D_d2:.3f}")

cavity_entry_diam = D_d2 * oxid_diam                                    # mm
print(f"Cavity Entrance Diameter: {cavity_entry_diam:.3f} mm")


# Picking L/D of 6
L_D = 6                                                    
cavity_len = L_D * cavity_entry_diam                                    # mm
print(f"Cavity Length: {cavity_len:.3f} mm")

# Picking angle of cavity to be 4
# Source - https://www.tandfonline.com/doi/full/10.1080/00102202.2011.589874?scroll=top&needAccess=true
cavity_angle = 4                                                        # degrees
cavity_tip_diam = cavity_entry_diam - 2 * cavity_len * math.tan(math.radians(cavity_angle / 2))  # mm

print(f"Cavity Tip Diameter: {cavity_tip_diam:.3f} mm")