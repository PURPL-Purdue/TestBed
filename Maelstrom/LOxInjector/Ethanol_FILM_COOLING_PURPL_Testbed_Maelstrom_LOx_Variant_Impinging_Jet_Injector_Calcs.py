##### PURPL Testbed Maelstrom LOx Variant Impinging Jet Injector -- Hand Calcs for Unlike Triplet Configuration #####
# use params.yaml from Turbopump for maximum LOx & RP-1 (Ethanol) flow rates (& likely densities of these propellants)
# use params.yaml from Testbed (eventually) for chamber pressure, chamber size, etc.
# main assumptions:
# - incompressible flow for both RP-1 & LOx (both liquids)

### Code to import params.yaml stuff once Testbed makes a Zephyr one
# from distutils.command.install import install
# import pip
# pip.main(['install', 'pyyaml'])

# import yaml

# with open('params.yaml', 'r') as file:
#     params = yaml.safe_load(file)

# # Access values, for example:

# learning_rate = params['train']['learning_rate']
# print(learning_rate)

# Import stuff
import math
import CoolProp.CoolProp as CoolProp
from CoolProp.CoolProp import PropsSI
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

print("##### PURPL Testbed Maelstrom LOx Variant Unlike Triplet Impinging Jet Injector #####\n")
#print("Note: Chamber Pressure, O/F Ratio, and Total Mass Flowrate change from Zephyr Engine Architecture updates. RP1 Pressure Drop may also change based on Regen team.\n")

## Define Parameters
m_to_in = 39.3701; # in/m
m2_to_in2 = 1550; # in^2/m^2
kgs_to_lbm = 2.20462; # lbm/kg
psia_to_Pa = 6894.76; # Pa/psia
in_to_ft = 1/12; # ft/in
ft_to_m = (1/in_to_ft)/m_to_in; # m/ft
g = 9.81; #m/s^2
lbs_to_slug = 32.174; #slug*ft/s^2/lb
N_to_lbf = 0.22481; #lbf/N
impingment_distance_to_orifice_diameter_ratio = 5; #unitless, NASA range suggests 3-7. need a better reason for selecting a value tho
print(f"Impingement Distance to Orifice Diameter Ratio: {impingment_distance_to_orifice_diameter_ratio}")
resultant_angle= 0; #degrees, angle between resultant impinged fluid stream & normal line from injector face
estimated_discharge_coefficient_thick_plate = 0.70; #unitless; assumed based on thick plate, square edge
estimated_discharge_coefficient_square_shoulder = 0.61; #unitless; assumed based on thin plate, square shoulder
#actual discharge coefficient will be determined experimentally - c_d_actual = mdot_measured/mdot_calculated
num_injector_elements = 32; #unitless
print(f"Number of Injector Elements: {num_injector_elements}")

### MAELSTROM ENGINE ARCHITECTURE PARAMETERS -- THESE CAN CHANGE
chamber_pressure_imperial = 250; #psia (TBD based on updates)
chamber_pressure_pa = chamber_pressure_imperial*psia_to_Pa; #Pa
print(f"Chamber Pressure: {chamber_pressure_imperial} psia")
OF_ratio = 1.3; #unitless (TBD based on Maelstrom GOx testing updates)
print(f"O/F Ratio: {OF_ratio}")
#thrust = 500 #lbf (THIS IS MAELSTROM GOx VARIANT THRUST -- WE NEED TO CALCULATE OUR OWN)
##CEA stuff
cea_obj = CEA_Obj(oxName = 'LOX', fuelName = 'Ethanol')
# gets CEA list output
output = cea_obj.get_PambCf(Pamb = 14.7, Pc = chamber_pressure_imperial, MR = OF_ratio, eps = 3.05) #eps = A_exit/A_throat from Maelstrom.yaml
# get force coefficient
cf = output[0]# get c* (need to convert ft/s to m/s)
cstar = (0.846)*cea_obj.get_Cstar(Pc = chamber_pressure_imperial, MR = OF_ratio) * ft_to_m #m/s (84.6% cstar efficiency from GG testing)
# calculate isp
isp = cstar * cf #s
# get exhaust velocity
v_e = isp # m/s (assumed 100% isp efficiency)
print(f"Exhaust velocity: {round(v_e*m_to_in*in_to_ft, 5)} ft/s")
#insert throat area from Maelstrom .yaml
throat_area_imperial = 1.373; #in^2 (from Maelstrom .yaml -- this is a fixed constraint we cannot change (it's Maelstrom's geometry))
throat_area = throat_area_imperial*(1/m_to_in**2); #m^2
print(f"Throat Area: {throat_area_imperial} in^2")
#get m_dot_total based on throat area
m_dot_total = chamber_pressure_pa * throat_area / cstar #kg/s
#get thrust based on that m_dot_total
thrust = m_dot_total * v_e #N
thrust_imperial = thrust*N_to_lbf; #lbf
print(f"Total Mass Flowrate: {round(m_dot_total*kgs_to_lbm, 5)} lbm/s")
print(f"Calculated Thrust: {round(thrust_imperial, 5)} lbf\n")


##Ethanol Parameters
#density_RP1 = 789; #kg/m^3
mdot_Ethanol = m_dot_total/(1 + OF_ratio); #kg/s
print(f"Ethanol Mass Flowrate: {round(mdot_Ethanol*kgs_to_lbm, 5)} lbm/s")
nonFC_mdot_Ethanol = 0.91*mdot_Ethanol; #kg/s (assuming 9% of Ethanol flow is diverted to film cooling)
FC_mdot_Ethanol = 0.09*mdot_Ethanol; #kg/s (9% of Ethanol flow for film cooling)
#mdot_RP1 = 7.4; #lbm/s (from Turbopump pump mdots - CHECK TURBOPUMP PUMPS_PARAMS.YAML FOR CHANGES)
flow_coeff_Ethanol_first_approx = 0.62; #
Ethanol_injector_stiffness = 0.19; #unitless
print(f"Ethanol Injector Stiffness: {Ethanol_injector_stiffness}")
Ethanol_regen_channel_pressure_imperial = chamber_pressure_imperial*(1+Ethanol_injector_stiffness); #psia (from Zephyr regen people)
print(f"Ethanol Feed Pressure: {Ethanol_regen_channel_pressure_imperial} psia")
Ethanol_regen_channel_pressure = Ethanol_regen_channel_pressure_imperial*psia_to_Pa; #Pa
Ethanol_injector_pressure_drop = Ethanol_regen_channel_pressure - chamber_pressure_pa; #units
Ethanol_regen_exit_temp = 400; #K (from James regen code for Maelstrom GOx variant)
density_Ethanol = PropsSI('D', 'T', Ethanol_regen_exit_temp, 'P', Ethanol_regen_channel_pressure, 'Ethanol') #kg/m^3 (using CoolProp Ethanol properties)
print(f"Ethanol Density at {Ethanol_regen_exit_temp} K and {Ethanol_regen_channel_pressure_imperial} psia: {round(density_Ethanol, 2)} kg/m^3")
dynamic_viscosity_Ethanol = PropsSI('V', 'T', Ethanol_regen_exit_temp, 'P', Ethanol_regen_channel_pressure, 'Ethanol') #Pa*s (using CoolProp Ethanol properties at the specified temperature and feed pressure)
print(f"Ethanol Dynamic Viscosity at 400 K and {Ethanol_regen_channel_pressure_imperial} psia: {dynamic_viscosity_Ethanol} Pa*s")
num_Ethanol_orifices = 2*num_injector_elements; #unitless
print(f"Ethanol Orifice Count: {num_Ethanol_orifices}")

##LOx Parameters
#density_LOx = 1237; #kg/m^3 (@ 70 K - around boiling point)
mdot_LOx = mdot_Ethanol * OF_ratio; #kg/s
print(f"LOx Mass Flowrate: {round(mdot_LOx*kgs_to_lbm, 5)} lbm/s")
#mdot_LOx = 15.4; #lbm/s (from Turbopump pump mdots - CHECK TURBOPUMP PUMPS_PARAMS.YAML FOR CHANGES)
flow_coeff_LOx_first_approx = 0.62; #unitless
LOx_feedline_pressure_imperial = 300; #psia (TBD based on updates)
LOX_injector_stiffness = (LOx_feedline_pressure_imperial - chamber_pressure_imperial)/chamber_pressure_imperial; #unitless
print(f"LOx Injector Stiffness: {LOX_injector_stiffness}")
print(f"LOx Feed Pressure: {LOx_feedline_pressure_imperial} psia")
LOx_feedline_pressure = LOx_feedline_pressure_imperial*psia_to_Pa; #Pa
LOx_injector_pressure_drop = LOx_feedline_pressure - chamber_pressure_pa; #units
temp_LOx = 90; #K (around boiling point of 90.19 K) 
density_LOx = PropsSI('D', 'T', temp_LOx, 'P', LOx_feedline_pressure, 'Oxygen') #kg/m^3 (using Oxygen at 90 K and LOx feed pressure)
print(f"LOx Density at 90 K and {LOx_feedline_pressure_imperial} psia: {round(density_LOx, 2)} kg/m^3")
dynamic_viscosity_LOx = PropsSI('V', 'T', temp_LOx, 'P', LOx_feedline_pressure, 'Oxygen') #Pa*s (using Oxygen at 90 K and LOx feed pressure)
print(f"LOx Dynamic Viscosity at 90 K and {LOx_feedline_pressure_imperial} psia: {round(dynamic_viscosity_LOx, 5)} Pa*s")
num_LOx_orifices = num_injector_elements; #unitless
print(f"LOx Orifice Count: {num_LOx_orifices}\n")


'''
ignore this for now - just use for reference if needed later
### First Approximation Calculations -- assuming incompressible flow for liquids (LOx & RP-1)
print("----First Approximation Calculations (Flow Coefficient = 0.62)----")
# RP-1 Orifice Diameter
RP1_total_orifice_area = (mdot_RP1/flow_coeff_RP1_first_approx)*math.sqrt(1/(2*density_RP1*RP1_injector_pressure_drop)); #m^2
RP1_single_orifice_area = RP1_total_orifice_area/num_RP1_orifices; #m^2
RP1_orifice_diameter = math.sqrt((4*RP1_single_orifice_area)/math.pi); #m
print(f"RP1 Orifice Diameter: {round(RP1_orifice_diameter*m_to_in, 5)} in")

# LOx Orifice Diameter
LOx_total_orifice_area = (mdot_LOx/flow_coeff_LOx_first_approx)*math.sqrt(1/(2*density_LOx*LOx_injector_pressure_drop)); #m^2
LOx_single_orifice_area = LOx_total_orifice_area/num_LOx_orifices; #m^2
LOx_orifice_diameter = math.sqrt((4*LOx_single_orifice_area)/math.pi); #m
print(f"LOx Orifice Diameter: {round(LOx_orifice_diameter*m_to_in, 5)} in")

# NASA Recommended Impingement Distance
average_orifice_diameter = (LOx_orifice_diameter + RP1_orifice_diameter)/2; #m
NASA_recommended_impingement_distance = impingment_distance_to_orifice_diameter_ratio*average_orifice_diameter; #m
print(f"NASA Recommended Impingement Distance: {round(NASA_recommended_impingement_distance*m_to_in, 5)} in")

# RP-1 Injector Angle Calcs
theta_RP1 = 45; #degrees
theta_LOx = 0; #degrees (assuming LOx spits straight out, perpendicular to injector face)
orifice_distance = NASA_recommended_impingement_distance*(math.tan(theta_RP1*(math.pi/180)) + math.tan(theta_RP1*(math.pi/180))); #m
print(f"Orifice Distance: {round(orifice_distance*m_to_in, 5)} in\n")
'''


### C_d*A Method Calculations -- assuming incompressible flow for liquids (LOx & RP-1)
# Thick Plate, Square Edge Orifice
print("----C_d*A Method Calculations----")
print(f"Assuming Estimated Discharge Coefficient of {estimated_discharge_coefficient_thick_plate} for a Thick Plate, Square Edge Orifice")
# Ethanol Orifice Diameter
Ethanol_CdA_total = (nonFC_mdot_Ethanol)*math.sqrt(1/(2*density_Ethanol*Ethanol_injector_pressure_drop)); #m^2
print(f"Ethanol Total CdA: {round(Ethanol_CdA_total*m2_to_in2, 5)} in^2")
Ethanol_total_orifice_area = Ethanol_CdA_total/estimated_discharge_coefficient_thick_plate; #m^2
Ethanol_min_inlet_diameter = math.sqrt((4*((mdot_Ethanol)*math.sqrt(1/(2*density_Ethanol*Ethanol_injector_pressure_drop)))/estimated_discharge_coefficient_thick_plate)/(math.pi)); #m
print(f"Ethanol Minimum Inlet Diameter: {round(Ethanol_min_inlet_diameter*m_to_in, 5)} in")
Ethanol_single_orifice_area = Ethanol_total_orifice_area/num_Ethanol_orifices; #m^2
Ethanol_orifice_diameter = math.sqrt((4*Ethanol_single_orifice_area)/math.pi); #m
print(f"Ethanol Orifice Diameter: {round(Ethanol_orifice_diameter*m_to_in, 5)} in")

# LOx Orifice Diameter
LOx_CdA_total = (mdot_LOx)*math.sqrt(1/(2*density_LOx*LOx_injector_pressure_drop)); #m^2
print(f"LOx Total CdA: {round(LOx_CdA_total*m2_to_in2, 5)} in^2")
LOx_total_orifice_area = LOx_CdA_total/estimated_discharge_coefficient_thick_plate; #m^2
LOx_min_inlet_diameter = math.sqrt((4*LOx_total_orifice_area)/(math.pi)); #m
print(f"LOx Minimum Inlet Diameter: {round(LOx_min_inlet_diameter*m_to_in, 5)} in")
LOx_single_orifice_area = LOx_total_orifice_area/num_LOx_orifices; #m^2
LOx_orifice_diameter = math.sqrt((4*LOx_single_orifice_area)/math.pi); #m
print(f"LOx Orifice Diameter: {round(LOx_orifice_diameter*m_to_in, 5)} in")

# NASA Recommended Impingement Distance
average_orifice_diameter = (LOx_orifice_diameter + Ethanol_orifice_diameter)/2; #m
NASA_recommended_impingement_distance = impingment_distance_to_orifice_diameter_ratio*average_orifice_diameter; #m
print(f"NASA Recommended Impingement Distance: {round(NASA_recommended_impingement_distance*m_to_in, 5)} in")

# Ethanol Injector Angle Calcs
theta_Ethanol = 22.5; #degrees
print(f"Ethanol Individual Injector Angle: {theta_Ethanol} degrees")
theta_LOx = 0; #degrees (assuming LOx spits straight out, perpendicular to injector face)
orifice_distance = NASA_recommended_impingement_distance*(math.tan(theta_Ethanol*(math.pi/180))); #m
print(f"Orifice Distance: {round(orifice_distance*m_to_in, 5)} in")



#LOx Lines Pressure Drop Calcs
print("\n----LOx Line Pressure Drop Calculations----")
velocity_single_LOx_line = (mdot_LOx/num_LOx_orifices)/(density_LOx*LOx_single_orifice_area); #m/s
print(f"LOx Velocity in Single Line: {round(velocity_single_LOx_line*m_to_in*in_to_ft, 5)} ft/s")
print(f"LOx Density at 90 K and {LOx_feedline_pressure_imperial} psia: {round(density_LOx*(kgs_to_lbm/(m_to_in*1/12)**3), 3)} lb/ft^3")
print(f"LOx Dynamic Viscosity at 90 K and {LOx_feedline_pressure_imperial} psia: {dynamic_viscosity_LOx} Pa*s")
Re_single_LOx_line = (density_LOx*velocity_single_LOx_line*LOx_orifice_diameter)/dynamic_viscosity_LOx; #unitless
print(f"LOx Reynolds Number in Single Line: {round(Re_single_LOx_line, 2)}")
average_roughness_commerical_steel = 0.0000015; #m, average roughness of commercial steel
print(f"Average Roughness of Commercial Steel: {round(average_roughness_commerical_steel*m_to_in, 5)} in")
relative_roughness_LOx_line = average_roughness_commerical_steel/LOx_orifice_diameter; #unitless, assuming commercial steel pipe with roughness of 0.15 mm
print(f"LOx Relative Roughness in Single Line: {round(relative_roughness_LOx_line, 8)}")
#friction_factor_moody_chart = 0.01875; #unitless, estimated from Moody chart

friction_factor_LOx_Haaland = 1/((-1.8*math.log10((6.9/Re_single_LOx_line)+((relative_roughness_LOx_line/3.7)**1.11)))**2); #unitless, Haaland equation
#Tkachenko, Mileikovskyi calculations
TM_A_0_LOx = -0.79638*(math.log((relative_roughness_LOx_line/8.208)+(7.3357/Re_single_LOx_line))); #unitless
TM_A_1_LOx = (Re_single_LOx_line*relative_roughness_LOx_line) + (9.3120665*TM_A_0_LOx); #unitless
friction_factor_LOx_TM = ((8.128943+TM_A_1_LOx)/((8.128943*TM_A_0_LOx)-(.86859209*TM_A_1_LOx*math.log(TM_A_1_LOx/(3.7099535*Re_single_LOx_line)))))**2; #unitless

#print(f"LOx Friction Factor from Moody Chart: {friction_factor_moody_chart}")
print(f"LOx Friction Factor from Haaland Equation: {round(friction_factor_LOx_Haaland, 5)}")
print(f"LOx Friction Factor from Tkachenko, Mileikovskyi: {round(friction_factor_LOx_TM, 5)}")
length_LOx_line = (0.525)*(1/m_to_in); #0.525in to m (longest line length) (inner ring to outer ring: (.422"+.103"=.525", .381"+.114"=.495", .487", .285"+.204"=.489", .367"+.134"=.501")
print(f"LOx Line Length: {round(length_LOx_line*m_to_in, 5)} in")
#head_loss_LOx_line_Moody = friction_factor_moody_chart*(length_LOx_line/LOx_orifice_diameter)*((velocity_single_LOx_line**2)/(2*g)); #m
head_loss_LOx_line_Haaland = friction_factor_LOx_Haaland*(length_LOx_line/LOx_orifice_diameter)*((velocity_single_LOx_line**2)/(2*g)); #m
head_loss_LOx_line_TM = friction_factor_LOx_TM*(length_LOx_line/LOx_orifice_diameter)*((velocity_single_LOx_line**2)/(2*g)); #m
#print(f"LOx Head Loss in Single Line from Moody Chart: {round(head_loss_LOx_line_Moody*m_to_in, 5)} in")
print(f"LOx Head Loss in Single Line from Haaland Equation: {round(head_loss_LOx_line_Haaland*m_to_in, 5)} in")
print(f"LOx Head Loss in Single Line from Tkachenko, Mileikovskyi: {round(head_loss_LOx_line_TM*m_to_in, 5)} in")
#pressure_loss_LOx_line_Moody = head_loss_LOx_line_Moody * density_LOx * g; #Pa
pressure_loss_LOx_line_Haaland = head_loss_LOx_line_Haaland * density_LOx * g; #Pa
pressure_loss_LOx_line_TM = head_loss_LOx_line_TM * density_LOx * g; #Pa
#print(f"LOx Pressure Loss in Single Line from Moody Chart: {round(pressure_loss_LOx_line_Moody*(1/psia_to_Pa), 5)} psia")
print(f"LOx Pressure Loss in Single Line from Haaland Equation: {round(pressure_loss_LOx_line_Haaland*(1/psia_to_Pa), 5)} psia")
print(f"LOx Pressure Loss in Single Line from Tkachenko, Mileikovskyi: {round(pressure_loss_LOx_line_TM*(1/psia_to_Pa), 5)} psia\n")



#Ethanol Lines Pressure Drop Calcs
print("\n----Ethanol Line Pressure Drop Calculations----")
velocity_single_Ethanol_line = (mdot_Ethanol/num_Ethanol_orifices)/(density_Ethanol*Ethanol_single_orifice_area); #m/s
print(f"Ethanol Velocity in Single Line: {round(velocity_single_Ethanol_line*m_to_in*in_to_ft, 5)} ft/s")
print(f"Ethanol Density at 90 K and {Ethanol_regen_channel_pressure_imperial} psia: {round(density_Ethanol*(kgs_to_lbm/(m_to_in*1/12)**3), 3)} lb/ft^3")
print(f"Ethanol Dynamic Viscosity at 90 K and {Ethanol_regen_channel_pressure_imperial} psia: {dynamic_viscosity_Ethanol} Pa*s")
Re_single_Ethanol_line = (density_Ethanol*velocity_single_Ethanol_line*Ethanol_orifice_diameter)/dynamic_viscosity_Ethanol; #unitless
print(f"Ethanol Reynolds Number in Single Line: {round(Re_single_Ethanol_line, 2)}")
average_roughness_commerical_steel = 0.0000015; #m, average roughness of commercial steel
print(f"Average Roughness of Commercial Steel: {round(average_roughness_commerical_steel*m_to_in, 5)} in")
relative_roughness_Ethanol_line = average_roughness_commerical_steel/Ethanol_orifice_diameter; #unitless, assuming commercial steel pipe with roughness of 0.15 mm
print(f"Ethanol Relative Roughness in Single Line: {round(relative_roughness_Ethanol_line, 8)}")
#friction_factor_moody_chart = 0.01875; #unitless, estimated from Moody chart

friction_factor_Ethanol_Haaland = 1/((-1.8*math.log10((6.9/Re_single_Ethanol_line)+((relative_roughness_Ethanol_line/3.7)**1.11)))**2); #unitless, Haaland equation
#Tkachenko, Mileikovskyi calculations
TM_A_0_ethanol = -0.79638*(math.log((relative_roughness_Ethanol_line/8.208)+(7.3357/Re_single_Ethanol_line))); #unitless
TM_A_1_ethanol = (Re_single_Ethanol_line*relative_roughness_Ethanol_line) + (9.3120665*TM_A_0_ethanol); #unitless
friction_factor_Ethanol_TM = ((8.128943+TM_A_1_ethanol)/((8.128943*TM_A_0_ethanol)-(.86859209*TM_A_1_ethanol*math.log(TM_A_1_ethanol/(3.7099535*Re_single_Ethanol_line)))))**2; #unitless

#print(f"LOx Friction Factor from Moody Chart: {friction_factor_moody_chart}")
print(f"Ethanol Friction Factor from Haaland Equation: {round(friction_factor_Ethanol_Haaland, 5)}")
print(f"Ethanol Friction Factor from Tkachenko, Mileikovskyi: {round(friction_factor_Ethanol_TM, 5)}")
length_Ethanol_line = 0.333*(1/m_to_in); #0.333in to m, longest length of Ethanol feed line from distribution holes to injector face (inner ring to outer ring:.293", .309", .318", .324", .333")
print(f"Ethanol Line Length: {round(length_Ethanol_line*m_to_in, 5)} in")
#head_loss_Ethanol_line_Moody = friction_factor_moody_chart*(length_Ethanol_line/Ethanol_orifice_diameter)*((velocity_single_Ethanol_line**2)/(2*g)); #m
head_loss_Ethanol_line_Haaland = friction_factor_Ethanol_Haaland*(length_Ethanol_line/Ethanol_orifice_diameter)*((velocity_single_Ethanol_line**2)/(2*g)); #m
head_loss_Ethanol_line_TM = friction_factor_Ethanol_TM*(length_Ethanol_line/Ethanol_orifice_diameter)*((velocity_single_Ethanol_line**2)/(2*g)); #m
#print(f"Ethanol Head Loss in Single Line from Moody Chart: {round(head_loss_Ethanol_line_Moody*m_to_in, 5)} in")
print(f"Ethanol Head Loss in Single Line from Haaland Equation: {round(head_loss_Ethanol_line_Haaland*m_to_in, 5)} in")
print(f"Ethanol Head Loss in Single Line from Tkachenko, Mileikovskyi: {round(head_loss_Ethanol_line_TM*m_to_in, 5)} in")
#pressure_loss_Ethanol_line_Moody = head_loss_Ethanol_line_Moody * density_Ethanol * g; #Pa
pressure_loss_Ethanol_line_Haaland = head_loss_Ethanol_line_Haaland * density_Ethanol * g; #Pa
pressure_loss_Ethanol_line_TM = head_loss_Ethanol_line_TM * density_Ethanol * g; #Pa
#print(f"Ethanol Pressure Loss in Single Line from Moody Chart: {round(pressure_loss_Ethanol_line_Moody*(1/psia_to_Pa), 5)} psia")
print(f"Ethanol Pressure Loss in Single Line from Haaland Equation: {round(pressure_loss_Ethanol_line_Haaland*(1/psia_to_Pa), 5)} psia")
print(f"Ethanol Pressure Loss in Single Line from Tkachenko, Mileikovskyi: {round(pressure_loss_Ethanol_line_TM*(1/psia_to_Pa), 5)} psia\n")



#LOx Manifold Inlet Velocity
LOx_manifold_inlet_width = (.474)*(1/m_to_in); #0.474in to m
LOx_manifold_inlet_min_height = (.100)*(1/m_to_in); #0.10in to m
LOx_manifold_inlet_max_height = (.500)*(1/m_to_in); #0.5in to m
LOx_manifold_min_area = LOx_manifold_inlet_width*LOx_manifold_inlet_min_height; #m^2
LOx_manifold_max_area = LOx_manifold_inlet_width*LOx_manifold_inlet_max_height; #m^2
min_velocity_LOx_manifold_inlet = (mdot_LOx)/(density_LOx*LOx_manifold_max_area); #m/s
max_velocity_LOx_manifold_inlet = (mdot_LOx)/(density_LOx*LOx_manifold_min_area); #m/s
print(f"Minimum LOx Manifold Inlet Velocity: {round(min_velocity_LOx_manifold_inlet*m_to_in*in_to_ft, 5)} ft/s")
print(f"Maximum LOx Manifold Inlet Velocity: {round(max_velocity_LOx_manifold_inlet*m_to_in*in_to_ft, 5)} ft/s\n")


## Ethanol Film cooling calcs
print("----Ethanol Film Cooling Calculations----")
num_Ethanol_FC_orifices = 10; #unitless
Ethanol_FC_CdA_total = (FC_mdot_Ethanol)*math.sqrt(1/(2*density_Ethanol*Ethanol_injector_pressure_drop)); #m^2
print(f"Ethanol Total CdA: {round(Ethanol_FC_CdA_total*m2_to_in2, 5)} in^2")
Ethanol_FC_total_orifice_area = Ethanol_FC_CdA_total/estimated_discharge_coefficient_thick_plate; #m^2
Ethanol_FC_single_orifice_area = Ethanol_FC_total_orifice_area/num_Ethanol_FC_orifices; #m^2
Ethanol_FC_orifice_diameter = math.sqrt((4*Ethanol_FC_single_orifice_area)/math.pi); #m
print(f"Ethanol Film Cooling Orifice Diameter: {round(Ethanol_FC_orifice_diameter*m_to_in, 5)} in")




'''
# Square Shoulder, Thin Plate Orifice
print("----C_d*A Method Calculations----")
print(f"Assuming Estimated Discharge Coefficient of {estimated_discharge_coefficient_square_shoulder} for a Square Shoulder, Thin Plate Orifice")
# RP-1 Orifice Diameter
RP1_CdA_total = (mdot_RP1)*math.sqrt(1/(2*density_RP1*RP1_injector_pressure_drop)); #m^2
print(f"RP1 Total CdA: {round(RP1_CdA_total*m2_to_in2, 5)} in^2")
RP1_total_orifice_area = RP1_CdA_total/estimated_discharge_coefficient_square_shoulder; #m^2
RP1_single_orifice_area = RP1_total_orifice_area/num_RP1_orifices; #m^2
RP1_orifice_diameter = math.sqrt((4*RP1_single_orifice_area)/math.pi); #m
print(f"RP1 Orifice Diameter: {round(RP1_orifice_diameter*m_to_in, 5)} in")

# LOx Orifice Diameter
LOx_CdA_total = (mdot_LOx)*math.sqrt(1/(2*density_LOx*LOx_injector_pressure_drop)); #m^2
print(f"LOx Total CdA: {round(LOx_CdA_total*m2_to_in2, 5)} in^2")
LOx_total_orifice_area = LOx_CdA_total/estimated_discharge_coefficient_square_shoulder; #m^2
LOx_single_orifice_area = LOx_total_orifice_area/num_LOx_orifices; #m^2
LOx_orifice_diameter = math.sqrt((4*LOx_single_orifice_area)/math.pi); #m
print(f"LOx Orifice Diameter: {round(LOx_orifice_diameter*m_to_in, 5)} in")

# NASA Recommended Impingement Distance
average_orifice_diameter = (LOx_orifice_diameter + RP1_orifice_diameter)/2; #m
NASA_recommended_impingement_distance = impingment_distance_to_orifice_diameter_ratio*average_orifice_diameter; #m
print(f"NASA Recommended Impingement Distance: {round(NASA_recommended_impingement_distance*m_to_in, 5)} in")

# RP-1 Injector Angle Calcs
theta_RP1 = 30; #degrees
print(f"RP1 Individual Injector Angle: {theta_RP1} degrees")
theta_LOx = 0; #degrees (assuming LOx spits straight out, perpendicular to injector face)
orifice_distance = NASA_recommended_impingement_distance*(math.tan(theta_RP1*(math.pi/180))); #m
print(f"Orifice Distance: {round(orifice_distance*m_to_in, 5)} in")
'''

# Number of Monte Carlo trials
n = 100000

# Assumption:
sigma = 3
Ethanol_orifice_tolerance = 0.0005   # in
LOx_orifice_tolerance = 0.0005        # in
flow_tolerance_percent = 2 # %
orifice_angle_tolerance = 0.1 # degrees
Weber_mismatch_limit = 3 # %
Ethanol_LOx_Weber_mismatch_limit = 2 # %

# Convert diameters from meters to inches
Ethanol_orifice_diameter_in = (Ethanol_orifice_diameter * m_to_in)
LOx_orifice_diameter_in = (LOx_orifice_diameter * m_to_in)

# Generate random manufactured orifice diameters
Ethanol_orifice_samples_in = np.random.normal(Ethanol_orifice_diameter_in, Ethanol_orifice_tolerance / sigma, (n, num_Ethanol_orifices))
LOx_orifice_samples_in = np.random.normal(LOx_orifice_diameter_in, LOx_orifice_tolerance / sigma, (n, num_LOx_orifices))

# Convert sample diameters from inches to meters
Ethanol_orifice_samples_m = Ethanol_orifice_samples_in / m_to_in
LOx_orifice_samples_m = LOx_orifice_samples_in / m_to_in

# Calculate actual area of each orifice
Ethanol_orifice_areas = (np.pi * Ethanol_orifice_samples_m**2 / 4)
LOx_orifice_areas = (np.pi * LOx_orifice_samples_m**2 / 4)

# Calculate CdA for each manufactured injector
Ethanol_orifice_CdA = (estimated_discharge_coefficient_thick_plate * Ethanol_orifice_areas)
LOx_orifice_CdA = (estimated_discharge_coefficient_thick_plate * LOx_orifice_areas)

# Calculate resulting mass flow
Ethanol_orifice_mdot = (Ethanol_orifice_CdA * np.sqrt(2 * density_Ethanol * Ethanol_injector_pressure_drop))
LOx_orifice_mdot = (LOx_orifice_CdA * np.sqrt(2 * density_LOx * LOx_injector_pressure_drop))

# Sum injector orifices
Ethanol_total_mdot = np.sum(Ethanol_orifice_mdot, axis=1)
LOx_total_mdot = np.sum(LOx_orifice_mdot, axis=1)

# Calculate main-injector O/F
main_injector_OF = (LOx_total_mdot / Ethanol_total_mdot)

# Print results
print("\nEthanol Orifice Diameter:")
print(f"Nominal: {Ethanol_orifice_diameter*m_to_in:.5f} in")
print(f"Mean:    {np.mean(Ethanol_orifice_samples_in):.5f} in")
print(f"Std Dev: {np.std(Ethanol_orifice_samples_in):.6f} in")
print(f"Min:     {np.min(Ethanol_orifice_samples_in):.5f} in")
print(f"Max:     {np.max(Ethanol_orifice_samples_in):.5f} in")

print("\nLOx Orifice Diameter:")
print(f"Nominal: {LOx_orifice_diameter*m_to_in:.5f} in")
print(f"Mean:    {np.mean(LOx_orifice_samples_in):.5f} in")
print(f"Std Dev: {np.std(LOx_orifice_samples_in):.6f} in")
print(f"Min:     {np.min(LOx_orifice_samples_in):.5f} in")
print(f"Max:     {np.max(LOx_orifice_samples_in):.5f} in")

print("\nMixture Ratio:")
print(f"Nominal: {OF_ratio:.5f}")
print(f"Mean:    {np.mean(main_injector_OF):.5f}")
print(f"Std Dev: {np.std(main_injector_OF):.5f}")
print(f"Min:     {np.min(main_injector_OF):.5f}")
print(f"Max:     {np.max(main_injector_OF):.5f}")

# Plot mixture ratio distribution
plt.figure(figsize=(10,5))
plt.hist(main_injector_OF,bins=100,color="#281E78")
plt.axvline(OF_ratio/0.91,color="red",linestyle="--",label="Nominal O/F Ratio")
plt.xlabel("O/F Ratio")
plt.ylabel("Number of Trials")
plt.title("Monte Carlo O/F Ratio Distribution")
plt.legend()
plt.show()

# Individual orifice flow statistics

# Ethanol
ethanol_mean_mdot = np.mean(Ethanol_orifice_mdot, axis=0)
ethanol_std_mdot = np.std(Ethanol_orifice_mdot, axis=0)
ethanol_min_mdot = np.min(Ethanol_orifice_mdot, axis=0)
ethanol_max_mdot = np.max(Ethanol_orifice_mdot, axis=0)
ethanol_cv = ethanol_std_mdot / ethanol_mean_mdot

# LOx
lox_mean_mdot = np.mean(LOx_orifice_mdot, axis=0)
lox_std_mdot = np.std(LOx_orifice_mdot, axis=0)
lox_min_mdot = np.min(LOx_orifice_mdot, axis=0)
lox_max_mdot = np.max(LOx_orifice_mdot, axis=0)
lox_cv = lox_std_mdot / lox_mean_mdot

# Group the two Ethanol orifices into each injector group
Ethanol_group_mdot = Ethanol_orifice_mdot.reshape(n, num_injector_elements, 2)

# Sum the two Ethanol orifices in each group
Ethanol_group_mdot = np.sum(Ethanol_group_mdot, axis=2)

# Each LOx orifice corresponds to one injector group
LOx_group_mdot = LOx_orifice_mdot

# Calculate combined mass flow for each injector group
group_mdot = Ethanol_group_mdot + LOx_group_mdot

# Calculate expected Ethanol and LOx flow for one group based on the target O/F
expected_Ethanol_group_mdot = nonFC_mdot_Ethanol / num_injector_elements
expected_LOx_group_mdot = mdot_LOx / num_injector_elements

# Calculate expected combined flow for one injector group
expected_group_mdot = expected_Ethanol_group_mdot + expected_LOx_group_mdot

# Set fixed engineering flow tolerance
#flow_tolerance_percent = 2 #%

# Calculate fixed flow limits
group_upper_flow_limit = expected_group_mdot * (1 + flow_tolerance_percent / 100)
group_lower_flow_limit = expected_group_mdot * (1 - flow_tolerance_percent / 100)

# Find groups above and below the fixed flow limits
group_above_flow_limit = group_mdot > group_upper_flow_limit
group_below_flow_limit = group_mdot < group_lower_flow_limit
group_outside_flow_limit = group_above_flow_limit | group_below_flow_limit

# Calculate number of groups above and below the limits for every Monte Carlo trial
groups_above_flow_limit = np.sum(group_above_flow_limit, axis=1)
groups_below_flow_limit = np.sum(group_below_flow_limit, axis=1)

# Calculate total number of group observations
total_group_observations = n * num_injector_elements

# Calculate percentage of group observations above and below the flow limits
percent_groups_above_flow_limit = np.sum(group_above_flow_limit) / total_group_observations * 100
percent_groups_below_flow_limit = np.sum(group_below_flow_limit) / total_group_observations * 100
percent_groups_outside_flow_limit = np.sum(group_outside_flow_limit) / total_group_observations * 100

# Print group flow results
print("\nInjector Group Flow Analysis:")
print(f"Number of Injector Groups: {num_injector_elements}")
print(f"Expected Ethanol Flow per Group: {expected_Ethanol_group_mdot:.6e} kg/s")
print(f"Expected LOx Flow per Group: {expected_LOx_group_mdot:.6e} kg/s")
print(f"Expected Combined Group Flow: {expected_group_mdot:.6e} kg/s")
print(f"Flow Tolerance: +/-{flow_tolerance_percent:.2f}%")
print(f"Flow Lower Limit: {group_lower_flow_limit:.6e} kg/s")
print(f"Flow Upper Limit: {group_upper_flow_limit:.6e} kg/s")
print(f"Groups Above Flow Limit: {np.sum(group_above_flow_limit)}")
print(f"Groups Below Flow Limit: {np.sum(group_below_flow_limit)}")
print(f"Percent Above Flow Limit: {percent_groups_above_flow_limit:.4f}%")
print(f"Percent Below Flow Limit: {percent_groups_below_flow_limit:.4f}%")
print(f"Percent Outside Flow Limits: {percent_groups_outside_flow_limit:.4f}%")

# Plot combined group flow distribution
plt.figure(figsize=(10,5))
plt.hist(group_mdot.flatten(), bins=100, color="#281E78")
plt.axvline(expected_group_mdot, color="black", linestyle="--", label="Expected Group Flow")
plt.axvline(group_upper_flow_limit, color="red", linestyle="--", label="+5% Flow Limit")
plt.axvline(group_lower_flow_limit, color="red", linestyle="--", label="-5% Flow Limit")
plt.xlabel("Combined Group Mass Flow Rate (kg/s)")
plt.ylabel("Number of Groups")
plt.title("Monte Carlo Injector Group Flow Distribution")
plt.legend()
plt.show()

# Investigate different tolerances and find the minimum acceptable tolerance for the orifices.
# Check impingement distance to see if it varies from other groups
# Check impingement angle
# Focus on a single group and investigate asymmetric collision

# Define the nominal orifice angle and manufacturing tolerance
orifice_angle = 60
#orifice_angle_tolerance = 0.1 # degrees look at start of code for assumptions

# Generate normally distributed angle variations for the two opposing orifices in each group
orifice_angle_left_samples = np.random.normal(orifice_angle, orifice_angle_tolerance / sigma, (n, num_injector_elements))
orifice_angle_right_samples = np.random.normal(orifice_angle, orifice_angle_tolerance / sigma, (n, num_injector_elements))

# Calculate the angular deviation of each orifice from the nominal 60 degree angle
orifice_angle_left_deviation = orifice_angle_left_samples - orifice_angle
orifice_angle_right_deviation = orifice_angle_right_samples - orifice_angle

# Calculate the resulting included impingement angle for each group
group_impingement_angle = orifice_angle + orifice_angle_left_deviation - orifice_angle_right_deviation

# Calculate the expected impingement angle
expected_impingement_angle = orifice_angle

# Calculate the standard deviation of the resulting impingement angle
group_impingement_angle_std = np.std(group_impingement_angle)

# Set fixed engineering impingement angle tolerance
angle_tolerance = 0.5 #degrees

# Calculate fixed impingement angle limits
group_impingement_angle_upper_limit = expected_impingement_angle + angle_tolerance
group_impingement_angle_lower_limit = expected_impingement_angle - angle_tolerance

# Find groups above and below the fixed angle limits
group_impingement_angle_above_limit = group_impingement_angle > group_impingement_angle_upper_limit
group_impingement_angle_below_limit = group_impingement_angle < group_impingement_angle_lower_limit
group_impingement_angle_outside_limit = group_impingement_angle_above_limit | group_impingement_angle_below_limit

# Calculate total number of group observations
total_angle_group_observations = n * num_injector_elements

# Calculate percentage of groups above and below the angle limits
percent_angle_groups_above_limit = np.sum(group_impingement_angle_above_limit) / total_angle_group_observations * 100
percent_angle_groups_below_limit = np.sum(group_impingement_angle_below_limit) / total_angle_group_observations * 100
percent_angle_groups_outside_limit = np.sum(group_impingement_angle_outside_limit) / total_angle_group_observations * 100

# Print angle results
print("\nInjector Group Impingement Angle Analysis:")
print(f"Expected Impingement Angle: {expected_impingement_angle:.4f} degrees")
print(f"Impingement Angle Std Dev: {group_impingement_angle_std:.4f} degrees")
print(f"Angle Tolerance: +/-{angle_tolerance:.4f} degrees")
print(f"Lower Angle Limit: {group_impingement_angle_lower_limit:.4f} degrees")
print(f"Upper Angle Limit: {group_impingement_angle_upper_limit:.4f} degrees")
print(f"Groups Above Angle Limit: {np.sum(group_impingement_angle_above_limit)}")
print(f"Groups Below Angle Limit: {np.sum(group_impingement_angle_below_limit)}")
print(f"Percent Above Angle Limit: {percent_angle_groups_above_limit:.4f}%")
print(f"Percent Below Angle Limit: {percent_angle_groups_below_limit:.4f}%")
print(f"Percent Outside Angle Limits: {percent_angle_groups_outside_limit:.4f}%")

# Plot impingement angle distribution
plt.figure(figsize=(10,5))
plt.hist(group_impingement_angle.flatten(), bins=100, color="#281E78")
plt.axvline(expected_impingement_angle, color="black", linestyle="--", label="Expected Impingement Angle")
plt.axvline(group_impingement_angle_upper_limit, color="red", linestyle="--", label="+0.5° Angle Limit")
plt.axvline(group_impingement_angle_lower_limit, color="red", linestyle="--", label="-0.5° Angle Limit")
plt.xlabel("Impingement Angle (degrees)")
plt.ylabel("Number of Groups")
plt.title("Monte Carlo Injector Group Impingement Angle Distribution")
plt.legend()
plt.show()

# Identify groups with unacceptable mass flow or impingement angle
group_outside_flow_or_angle = group_outside_flow_limit | group_impingement_angle_outside_limit

# Calculate total number of group observations outside either flow or angle limit
groups_outside_flow_or_angle = np.sum(group_outside_flow_or_angle)

# Calculate percentage of groups outside either flow or angle limit
percent_groups_outside_flow_or_angle = groups_outside_flow_or_angle / total_group_observations * 100

# Calculate groups outside both flow and angle limits
group_outside_both_flow_and_angle = group_outside_flow_limit & group_impingement_angle_outside_limit

# Calculate total number of group observations outside both limits
groups_outside_both_flow_and_angle = np.sum(group_outside_both_flow_and_angle)

# Calculate percentage of group observations outside both limits
percent_groups_outside_both_flow_and_angle = groups_outside_both_flow_and_angle / total_group_observations * 100

# Print combined flow and angle results
print("\nCombined Injector Group Flow and Angle Analysis:")
print(f"Groups Outside Flow or Angle Limits: {groups_outside_flow_or_angle}")
print(f"Percent Outside Flow or Angle Limits: {percent_groups_outside_flow_or_angle:.4f}%")
print(f"Groups Outside Both Flow and Angle Limits: {groups_outside_both_flow_and_angle}")
print(f"Percent Outside Both Flow and Angle Limits: {percent_groups_outside_both_flow_and_angle:.4f}%")

# Ethanol surface tension
Ethanol_surface_tension = PropsSI('I', 'T', Ethanol_regen_exit_temp, 'Q', 0, 'Ethanol')  # N/m - CoolProp saturated-liquid surface tension at the fuel temperature

# LOx surface tension
LOx_surface_tension = 0.0132  # N/m 

# Calculate nominal Ethanol velocity using the same orifice flow model
velocity_single_Ethanol_orifice = estimated_discharge_coefficient_thick_plate * np.sqrt(2 * Ethanol_injector_pressure_drop / density_Ethanol)

# Calculate nominal Ethanol Weber number
Nominal_weber = density_Ethanol * velocity_single_Ethanol_orifice**2 * Ethanol_orifice_diameter / Ethanol_surface_tension

print(f"Nominal Ethanol Weber Number: {Nominal_weber:.2f}")

# Calculate velocity through each Ethanol orifice
Ethanol_orifice_velocity = Ethanol_orifice_mdot / (density_Ethanol * Ethanol_orifice_areas)

#Calculate velocity through LOx orifice
LOx_orifice_velocity = LOx_orifice_mdot / (density_LOx * LOx_orifice_areas)

# Calculate Weber number for each Ethanol orifice
Ethanol_orifice_Weber = (density_Ethanol * Ethanol_orifice_velocity**2 * Ethanol_orifice_samples_m) / Ethanol_surface_tension

#Calculate Weber number for each LOx orifice
LOx_orifice_Weber = (density_LOx * LOx_orifice_velocity**2 * LOx_orifice_samples_m) / LOx_surface_tension

# Group the two Ethanol orifices into 32 injector groups
Ethanol_group_Weber = Ethanol_orifice_Weber.reshape(n, num_injector_elements, 2)

# Group LOx orifices
LOx_group_Weber = LOx_orifice_Weber

# Separate the two Ethanol Weber numbers in each injector group
Ethanol_Weber_1 = Ethanol_group_Weber[:, :, 0]
Ethanol_Weber_2 = Ethanol_group_Weber[:, :, 1]

# Calculate Weber number mismatch between the two Ethanol jets
Weber_mismatch = np.abs(Ethanol_Weber_1 - Ethanol_Weber_2) / ((Ethanol_Weber_1 + Ethanol_Weber_2) / 2) * 100

# Set the maximum acceptable Weber mismatch
#Weber_mismatch_limit = 5

# Identify injector groups with excessive Weber number mismatch
Weber_mismatch_outside_limit = Weber_mismatch > Weber_mismatch_limit

# Calculate percentage of injector groups with excessive Weber number mismatch
percent_Weber_mismatch_outside_limit = np.sum(Weber_mismatch_outside_limit) / Weber_mismatch_outside_limit.size * 100

# Ethanol Weber average
Ethanol_average_Weber = (Ethanol_Weber_1 + Ethanol_Weber_2) / 2

#Ethanol LOx Mismatch
Ethanol_LOx_Weber_mismatch = np.abs(Ethanol_average_Weber - LOx_group_Weber) / ((Ethanol_average_Weber + LOx_group_Weber) / 2) * 100

#Nominal Ethanol LOx Weber Mismatch
Nominal_LOx_velocity = estimated_discharge_coefficient_thick_plate * np.sqrt(2 * LOx_injector_pressure_drop / density_LOx)
Nominal_LOx_weber = density_LOx * Nominal_LOx_velocity**2 * LOx_orifice_diameter / LOx_surface_tension
Nominal_Ethanol_LOx_Weber_mismatch = np.abs(Nominal_weber - Nominal_LOx_weber) / ((Nominal_weber + Nominal_LOx_weber) / 2) * 100

#Calculate percent of injectors with ethanol LOx Weber mismatch outside of 2% of nominal
#Ethanol_LOx_Weber_mismatch_limit = 2 # %
Ethanol_LOx_Weber_mismatch_outside_upper_limit = Ethanol_LOx_Weber_mismatch > (Nominal_Ethanol_LOx_Weber_mismatch + (Ethanol_LOx_Weber_mismatch_limit))
Ethanol_LOx_Weber_mismatch_outside_lower_limit = Ethanol_LOx_Weber_mismatch < (Nominal_Ethanol_LOx_Weber_mismatch - (Ethanol_LOx_Weber_mismatch_limit))
Ethanol_LOx_Weber_mismatch_percent = (np.sum(Ethanol_LOx_Weber_mismatch_outside_upper_limit) + np.sum(Ethanol_LOx_Weber_mismatch_outside_lower_limit)) / Ethanol_LOx_Weber_mismatch.size * 100

# Print Weber number results
print("\nWeber Number Analysis:")
print(f"Mean Ethanol Weber Number: {np.mean(Ethanol_orifice_Weber):.2f}")
print(f"Std Dev Ethanol Weber Number: {np.std(Ethanol_orifice_Weber):.2f}")
print(f"Minimum Ethanol Weber Number: {np.min(Ethanol_orifice_Weber):.2f}")
print(f"Maximum Ethanol Weber Number: {np.max(Ethanol_orifice_Weber):.2f}")
print(f"Mean Weber Number Mismatch: {np.mean(Weber_mismatch):.2f}%")
print(f"Maximum Weber Number Mismatch: {np.max(Weber_mismatch):.2f}%")
print(f"Weber Mismatch Limit: {Weber_mismatch_limit:.2f}%")
print(f"Percent of Groups Above Weber Mismatch Limit: {percent_Weber_mismatch_outside_limit:.4f}%")
print()
print(f"Mean LOx Weber Number: {np.mean(LOx_group_Weber):.2f}")
print(f"Std Dev LOx Weber Number: {np.std(LOx_group_Weber):.2f}")
print(f"Minimum LOx Weber Number: {np.min(LOx_group_Weber):.2f}")
print(f"Maximum LOx Weber Number: {np.max(LOx_group_Weber):.2f}")
print(f"Percent of Groups with Ethanol-LOx Weber Mismatch Outside Limits: {Ethanol_LOx_Weber_mismatch_percent:.2f}%")
# Plot Weber number distribution

plt.figure(figsize=(10,5))
plt.hist(Ethanol_orifice_Weber.flatten(), bins=100, color="#281E78")
plt.axvline(np.mean(Ethanol_orifice_Weber), color="black", linestyle="--", label="Mean Weber Number")
plt.xlabel("Ethanol Weber Number")
plt.ylabel("Number of Orifices")
plt.title("Monte Carlo Ethanol Weber Number Distribution")
plt.legend()
plt.show()

plt.figure(figsize=(10,5))
plt.hist(LOx_group_Weber.flatten(), bins=100, color="#281E78")
plt.axvline(np.mean(LOx_group_Weber), color="black", linestyle="--", label="Mean LOx Weber Number")
plt.xlabel("LOx Weber Number")
plt.ylabel("Number of Injector Groups")
plt.title("Monte Carlo LOx Weber Number Distribution")
plt.legend()
plt.show()

# Plot Weber number mismatch distribution
plt.figure(figsize=(10,5))
plt.hist(Weber_mismatch.flatten(), bins=100, color="#281E78")
plt.axvline(Weber_mismatch_limit, color="red", linestyle="--", label="Maximum Acceptable Mismatch")
plt.axvline(np.mean(Weber_mismatch), color="black", linestyle="--", label="Mean Weber Mismatch")
plt.xlabel("Weber Number Mismatch (%)")
plt.ylabel("Number of Injector Groups")
plt.title("Monte Carlo Ethanol Weber Number Mismatch")
plt.legend()
plt.show()
#Ethanol LOx Weber Mismatch Distribution
plt.figure()
plt.hist(Ethanol_LOx_Weber_mismatch.flatten(), bins=50)
plt.axvline(Nominal_Ethanol_LOx_Weber_mismatch, color = "black", linestyle='--', linewidth=2, label='Nominal Ethanol-LOx mismatch')
plt.axvline(Nominal_Ethanol_LOx_Weber_mismatch + 2, color = "red", linestyle='--', linewidth=2, label='+2% Ethanol-LOx mismatch Upper Limit')
plt.axvline(Nominal_Ethanol_LOx_Weber_mismatch - 2, color = "red", linestyle='--', linewidth=2, label='-2% Ethanol-LOx mismatch Lower Limit')
plt.xlabel("Ethanol-LOx Weber Mismatch (%)")
plt.ylabel("Count")
plt.title("Ethanol-LOx Weber Mismatch")
plt.legend()
plt.show()

# Determine whether each complete injector violates any flow, angle, or Weber mismatch limit
injector_flow_outside_limit = np.any(group_outside_flow_limit, axis=1)
injector_angle_outside_limit = np.any(group_impingement_angle_outside_limit, axis=1)
injector_Weber_outside_limit = np.any(Weber_mismatch_outside_limit + Ethanol_LOx_Weber_mismatch_outside_upper_limit + Ethanol_LOx_Weber_mismatch_outside_lower_limit, axis=1)

# Determine whether each complete injector violates at least one limit
injector_outside_any_limit = injector_flow_outside_limit | injector_angle_outside_limit | injector_Weber_outside_limit

# Calculate overall injector pass/fail percentages
percent_injectors_outside_any_limit = np.sum(injector_outside_any_limit) / n * 100
percent_injectors_within_all_limits = 100 - percent_injectors_outside_any_limit
percent_injectors_outside_flow = np.sum(injector_flow_outside_limit) / n * 100
percent_injectors_outside_angle = np.sum(injector_angle_outside_limit) / n * 100
percent_injectors_outside_Weber = np.sum(injector_Weber_outside_limit) / n * 100

# Print overall injector performance results
print("\n-----Overall Injector Performance Analysis-----")
print(f"Percent of Injectors Outside Any Limit: {percent_injectors_outside_any_limit:.2f}%")
print(f"Percent of Injectors Within All Limits: {percent_injectors_within_all_limits:.2f}%")
print(f"Injectors with Flow Limit Violation: {percent_injectors_outside_flow:.2f}%")
print(f"Injectors with Angle Limit Violation: {percent_injectors_outside_angle:.2f}%")
print(f"Injectors with Weber Mismatch Violation: {percent_injectors_outside_Weber:.2f}%")

# Calculate number of limit categories violated by each complete injector
injector_number_of_violations = injector_flow_outside_limit.astype(int) + injector_angle_outside_limit.astype(int) + injector_Weber_outside_limit.astype(int)

# Plot number of limit categories violated by each complete injector
plt.figure(figsize=(10,5))
plt.hist(injector_number_of_violations,bins=np.arange(-0.5,4.5,1),rwidth=0.8)
plt.xlabel("Number of Limit Categories Violated")
plt.ylabel("Number of Injectors")
plt.title("Complete Injector Limit Violations")
plt.xticks([0,1,2,3])
plt.show()

# Find out what the weber number means
# Find how it changes with Ethanol

