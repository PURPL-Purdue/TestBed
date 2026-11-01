##### PURPL Testbed Maelstrom LOx Variant Impinging Jet Injector -- Hand Calcs for Unlike Triplet Configuration #####
# use params.yaml from Turbopump for maximum LOx & RP-1 (Kerosene) flow rates (& likely densities of these propellants)
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
cea_obj = CEA_Obj(oxName = 'LOX', fuelName = 'Kerosene')
# gets CEA list output
output = cea_obj.get_PambCf(Pamb = 14.7, Pc = chamber_pressure_imperial, MR = OF_ratio, eps = 3.05) #eps = A_exit/A_throat from Maelstrom.yaml
# get force coefficient
cf = output[0]
# get c* (need to convert ft/s to m/s)
cstar = (0.846)*cea_obj.get_Cstar(Pc = chamber_pressure_imperial, MR = OF_ratio) * ft_to_m #m/s (84.6% cstar efficiency from GG testing)
print(f"Characteristic velocity (C star): {round(cstar*m_to_in*in_to_ft, 5)} ft/s")
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


##Kerosene Parameters
#density_RP1 = 789; #kg/m^3
mdot_Kerosene = m_dot_total/(1 + OF_ratio); #kg/s
print(f"Kerosene Mass Flowrate: {round(mdot_Kerosene*kgs_to_lbm, 5)} lbm/s")
#mdot_RP1 = 7.4; #lbm/s (from Turbopump pump mdots - CHECK TURBOPUMP PUMPS_PARAMS.YAML FOR CHANGES)
flow_coeff_Kerosene_first_approx = 0.62; #
Kerosene_injector_stiffness = 0.19; #unitless
print(f"Kerosene Injector Stiffness: {Kerosene_injector_stiffness}")
Kerosene_regen_channel_pressure_imperial = chamber_pressure_imperial*(1+Kerosene_injector_stiffness); #psia (from Zephyr regen people)
print(f"Kerosene Feed Pressure: {Kerosene_regen_channel_pressure_imperial} psia")
Kerosene_regen_channel_pressure = Kerosene_regen_channel_pressure_imperial*psia_to_Pa; #Pa
Kerosene_injector_pressure_drop = Kerosene_regen_channel_pressure - chamber_pressure_pa; #units
Kerosene_regen_exit_temp = 400; #K (from James regen code for Maelstrom GOx variant)
density_Kerosene = PropsSI('D', 'T', Kerosene_regen_exit_temp, 'P', Kerosene_regen_channel_pressure, 'Dodecane') #kg/m^3 (using Dodecane to estimate Kerosene)
print(f"Kerosene Density at {Kerosene_regen_exit_temp} K and {Kerosene_regen_channel_pressure_imperial} psia: {round(density_Kerosene, 2)} kg/m^3")
dynamic_viscosity_Kerosene = PropsSI('V', 'T', Kerosene_regen_exit_temp, 'P', Kerosene_regen_channel_pressure, 'Dodecane') #Pa*s (using Dodecane at 400 K and Kerosene feed pressure)
print(f"Kerosene Dynamic Viscosity at 400 K and {Kerosene_regen_channel_pressure_imperial} psia: {dynamic_viscosity_Kerosene} Pa*s")
num_Kerosene_orifices = 2*num_injector_elements; #unitless
print(f"Kerosene Orifice Count: {num_Kerosene_orifices}")
prandtl_Kerosene = PropsSI('Prandtl', 'T', Kerosene_regen_exit_temp, 'P', Kerosene_regen_channel_pressure, 'Dodecane')
print(f"Kerosene Prandtl Number at 90 K and {Kerosene_regen_channel_pressure_imperial} psia: {round(prandtl_Kerosene, 5)}")


##LOx Parameters
#density_LOx = 1237; #kg/m^3 (@ 70 K - around boiling point)
mdot_LOx = mdot_Kerosene * OF_ratio; #kg/s
print(f"LOx Mass Flowrate: {round(mdot_LOx*kgs_to_lbm, 5)} lbm/s")
#mdot_LOx = 15.4; #lbm/s (from Turbopump pump mdots - CHECK TURBOPUMP PUMPS_PARAMS.YAML FOR CHANGES)
flow_coeff_LOx_first_approx = 0.62; #unitless
LOx_feedline_pressure_imperial = 300 ; #psia (TBD based on updates)
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
print(f"LOx Orifice Count: {num_LOx_orifices}")
average_roughness_commerical_steel = 0.0000015; #m, average roughness of commercial steel
print(f"Average Roughness of Commercial Steel: {round(average_roughness_commerical_steel*m_to_in, 5)} in")
prandtl_LOx = PropsSI('Prandtl', 'T', temp_LOx, 'P', LOx_feedline_pressure, 'Oxygen')
print(f"LOx Prandtl Number at 90 K and {LOx_feedline_pressure_imperial} psia: {round(prandtl_LOx, 5)}")

stoich_ratio = 37/2 # The stoiciometric ratio of O2 to Dodecane in combustion, with Kerosene approximated as Dodecane (2 C12H26 + 37 O2 --> 24 CO2 + 26 H2O)
kero_molar_mass =  170.34
LOx_molar_mass = 31.999
equivalence_ratio = stoich_ratio * LOx_molar_mass/kero_molar_mass / OF_ratio
print(f"Equivalence ratio: {round(equivalence_ratio, 5)}\n")


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
# Kerosene Orifice Diameter
Kerosene_CdA_total = (mdot_Kerosene)*math.sqrt(1/(2*density_Kerosene*Kerosene_injector_pressure_drop)); #m^2
print(f"Kerosene Total CdA: {round(Kerosene_CdA_total*m2_to_in2, 5)} in^2")
Kerosene_total_orifice_area = Kerosene_CdA_total/estimated_discharge_coefficient_thick_plate; #m^2
Kerosene_min_inlet_diameter = math.sqrt((4*Kerosene_total_orifice_area)/(math.pi)); #m
print(f"Kerosene Minimum Inlet Diameter: {round(Kerosene_min_inlet_diameter*m_to_in, 5)} in")
Kerosene_single_orifice_area = Kerosene_total_orifice_area/num_Kerosene_orifices; #m^2
Kerosene_orifice_diameter = math.sqrt((4*Kerosene_single_orifice_area)/math.pi); #m
print(f"Kerosene Orifice Diameter: {round(Kerosene_orifice_diameter*m_to_in, 5)} in")

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
average_orifice_diameter = (LOx_orifice_diameter + Kerosene_orifice_diameter)/2; #m
NASA_recommended_impingement_distance = impingment_distance_to_orifice_diameter_ratio*average_orifice_diameter; #m
print(f"NASA Recommended Impingement Distance: {round(NASA_recommended_impingement_distance*m_to_in, 5)} in")

# Kerosene Injector Angle Calcs
theta_Kerosene = 22.5; #degrees
print(f"Kerosene Individual Injector Angle: {theta_Kerosene} degrees")
theta_LOx = 0; #degrees (assuming LOx spits straight out, perpendicular to injector face)
orifice_distance = NASA_recommended_impingement_distance*(math.tan(theta_Kerosene*(math.pi/180))); #m
print(f"Orifice Distance: {round(orifice_distance*m_to_in, 5)} in")



#LOx Lines Pressure Drop Calcs

# Created as a function for future calculations that involve finding velocity and frictional losses in other plate distribution regions
# diameter (m): Diameter of single orifice in region of interest
# num_orifices (unitless): Number of orifices in region of interest
# region [string]: Name of region for printing purposes; not used in calculations
def LOx_analysis(diameter, num_orifices, region):
    print(f"\n----LOx {region} Calculations----")
    single_orifice_area = math.pi*(diameter/2)**2

    single_line_velocity = (mdot_LOx/num_orifices)/(density_LOx*single_orifice_area); #m/s
    print(f"LOx Velocity in {region}: {round(single_line_velocity*m_to_in*in_to_ft, 5)} ft/s")
    Re = (density_LOx*single_line_velocity*diameter)/dynamic_viscosity_LOx; #unitless
    print(f"LOx Reynolds Number in {region}: {round(Re, 2)}")

    relative_roughness_LOx_line = average_roughness_commerical_steel/diameter; #unitless, assuming commercial steel pipe with roughness of 0.15 mm
    print(f"LOx Relative Roughness in {region}: {round(relative_roughness_LOx_line, 8)}")
    #friction_factor_moody_chart = 0.01875; #unitless, estimated from Moody chart

    # Haaland calculations
    friction_factor_LOx_Haaland = 1/((-1.8*math.log10((6.9/Re)+((relative_roughness_LOx_line/3.7)**1.11)))**2); #unitless, Haaland equation
    #print(f"LOx Friction Factor from Moody Chart: {friction_factor_moody_chart}")
    print(f"LOx Friction Factor from Haaland Equation for {region}: {round(friction_factor_LOx_Haaland, 5)}")

    # Tkachenko, Mileikovskyi calculations
    TM_A_0_LOx = -0.79638*(math.log((relative_roughness_LOx_line/8.208)+(7.3357/Re))); #unitless
    TM_A_1_LOx = (Re*relative_roughness_LOx_line) + (9.3120665*TM_A_0_LOx); #unitless
    friction_factor_LOx_TM = ((8.128943+TM_A_1_LOx)/((8.128943*TM_A_0_LOx)-(.86859209*TM_A_1_LOx*math.log(TM_A_1_LOx/(3.7099535*Re)))))**2; #unitless
    print(f"LOx Friction Factor from Tkachenko, Mileikovskyi for {region}: {round(friction_factor_LOx_TM, 5)}")

    return single_line_velocity, Re, friction_factor_LOx_Haaland, friction_factor_LOx_TM

output = LOx_analysis(LOx_orifice_diameter, num_LOx_orifices, "Single Line")
velocity_single_LOx_line = output[0]
Re_single_LOx_line = output[1]
friction_factor_LOx_Haaland = output[2]
friction_factor_LOx_TM = output[3]

print(f"LOx Density at 90 K and {LOx_feedline_pressure_imperial} psia: {round(density_LOx*(kgs_to_lbm/(m_to_in*1/12)**3), 3)} lb/ft^3")
print(f"LOx Dynamic Viscosity at 90 K and {LOx_feedline_pressure_imperial} psia: {dynamic_viscosity_LOx} Pa*s")
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



#Kerosene Lines Pressure Drop Calcs

# Created as a function for future calculations that involve finding velocity and frictional losses in other plate distribution regions
# diameter (m): Diameter of single orifice in region of interest
# num_orifices (unitless): Number of orifices in region of interest
# region [string]: Name of region for printing purposes; not used in calculations
def Kerosene_analysis(diameter, num_orifices, region):
    print(f"\n----Kerosene {region} Calculations----")
    single_orifice_area = math.pi*(diameter/2)**2

    single_line_velocity = (mdot_Kerosene/num_orifices)/(density_Kerosene*single_orifice_area); #m/s
    print(f"Kerosene Velocity in {region}: {round(single_line_velocity*m_to_in*in_to_ft, 5)} ft/s")
    Re = (density_Kerosene*single_line_velocity*diameter)/dynamic_viscosity_Kerosene; #unitless
    print(f"Kerosene Reynolds Number in {region}: {round(Re, 2)}")

    relative_roughness_Kerosene_line = average_roughness_commerical_steel/diameter; #unitless, assuming commercial steel pipe with roughness of 0.15 mm
    print(f"Kerosene Relative Roughness in {region}: {round(relative_roughness_Kerosene_line, 8)}")
    #friction_factor_moody_chart = 0.01875; #unitless, estimated from Moody chart

    # Haaland calculations
    friction_factor_Kerosene_Haaland = 1/((-1.8*math.log10((6.9/Re)+((relative_roughness_Kerosene_line/3.7)**1.11)))**2); #unitless, Haaland equation
    #print(f"Kerosene Friction Factor from Moody Chart: {friction_factor_moody_chart}")
    print(f"Kerosene Friction Factor from Haaland Equation for {region}: {round(friction_factor_Kerosene_Haaland, 5)}")

    # Tkachenko, Mileikovskyi calculations
    TM_A_0_Kerosene = -0.79638*(math.log((relative_roughness_Kerosene_line/8.208)+(7.3357/Re))); #unitless
    TM_A_1_Kerosene = (Re*relative_roughness_Kerosene_line) + (9.3120665*TM_A_0_Kerosene); #unitless
    friction_factor_Kerosene_TM = ((8.128943+TM_A_1_Kerosene)/((8.128943*TM_A_0_Kerosene)-(.86859209*TM_A_1_Kerosene*math.log(TM_A_1_Kerosene/(3.7099535*Re)))))**2; #unitless
    print(f"Kerosene Friction Factor from Tkachenko, Mileikovskyi for {region}: {round(friction_factor_Kerosene_TM, 5)}")

    return single_line_velocity, Re, friction_factor_Kerosene_Haaland, friction_factor_Kerosene_TM

output = Kerosene_analysis(Kerosene_orifice_diameter, num_Kerosene_orifices, "Single Line")
velocity_single_Kerosene_line = output[0]
Re_single_Kerosene_line = output[1]
friction_factor_Kerosene_Haaland = output[2]
friction_factor_Kerosene_TM = output[3]

length_Kerosene_line = 0.333*(1/m_to_in); #0.333in to m, longest length of Kerosene feed line from distribution holes to injector face (inner ring to outer ring:.293", .309", .318", .324", .333")
print(f"Kerosene Line Length: {round(length_Kerosene_line*m_to_in, 5)} in")
#head_loss_Kerosene_line_Moody = friction_factor_moody_chart*(length_Kerosene_line/Kerosene_orifice_diameter)*((velocity_single_Kerosene_line**2)/(2*g)); #m
head_loss_Kerosene_line_Haaland = friction_factor_Kerosene_Haaland*(length_Kerosene_line/Kerosene_orifice_diameter)*((velocity_single_Kerosene_line**2)/(2*g)); #m
head_loss_Kerosene_line_TM = friction_factor_Kerosene_TM*(length_Kerosene_line/Kerosene_orifice_diameter)*((velocity_single_Kerosene_line**2)/(2*g)); #m
#print(f"Kerosene Head Loss in Single Line from Moody Chart: {round(head_loss_Kerosene_line_Moody*m_to_in, 5)} in")
print(f"Kerosene Head Loss in Single Line from Haaland Equation: {round(head_loss_Kerosene_line_Haaland*m_to_in, 5)} in")
print(f"Kerosene Head Loss in Single Line from Tkachenko, Mileikovskyi: {round(head_loss_Kerosene_line_TM*m_to_in, 5)} in")
#pressure_loss_Kerosene_line_Moody = head_loss_Kerosene_line_Moody * density_Kerosene * g; #Pa
pressure_loss_Kerosene_line_Haaland = head_loss_Kerosene_line_Haaland * density_Kerosene * g; #Pa
pressure_loss_Kerosene_line_TM = head_loss_Kerosene_line_TM * density_Kerosene * g; #Pa
#print(f"Kerosene Pressure Loss in Single Line from Moody Chart: {round(pressure_loss_Kerosene_line_Moody*(1/psia_to_Pa), 5)} psia")
print(f"Kerosene Pressure Loss in Single Line from Haaland Equation: {round(pressure_loss_Kerosene_line_Haaland*(1/psia_to_Pa), 5)} psia")
print(f"Kerosene Pressure Loss in Single Line from Tkachenko, Mileikovskyi: {round(pressure_loss_Kerosene_line_TM*(1/psia_to_Pa), 5)} psia\n")


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


### Calculate heat transfer coefficient

# Created as a function to retrieve the heat transfer at all different regions of the plate
# propellant [string]: Choice of propellant (LOx or Kerosene)
# diameter (m): Diameter of single orifice in region of interest
# num_orifices (unitless): Number of orifices in region of interest
# region [string]: Name of region for printing purposes; not used in calculations

def heat_transfer(propellant, diameter, num_orifices, region):
    
    if propellant == "LOx":
        # LOx paramters
        prandtl = prandtl_LOx
        thermal_conductivity = PropsSI('L', 'T', temp_LOx, 'P', LOx_feedline_pressure, 'Oxygen')
        # Reynold's number and frictional loss calculations
        output = LOx_analysis(diameter, num_orifices, region)
        Re = output[1]
        friction_factor = output[3] # Using TM approximation for Darcy friction factor
    elif propellant == "Kerosene":
        # Kerosene paramters
        prandtl = prandtl_Kerosene
        thermal_conductivity = PropsSI('L', 'T', Kerosene_regen_exit_temp, 'P', Kerosene_regen_channel_pressure, 'Dodecane')
        # Reynold's number and frictional loss calculations
        output = Kerosene_analysis(diameter, num_orifices, region)
        Re = output[1]
        friction_factor = output[3]

    # Gnielinski's Correlation using Darcy friction factor
    nusselt = (friction_factor / 8)*(Re - 1000)*prandtl / (1 + 12.7*(friction_factor/8)**(1/2)*(prandtl**(2/3) - 1))  
    heat_transfer_coeff = nusselt * thermal_conductivity / diameter
    return heat_transfer_coeff, prandtl


# Calculations for LOx

# Impinging lines
heat_transfer_coeff_LOx = heat_transfer("LOx", LOx_orifice_diameter, num_LOx_orifices, "Single Line")[0]
print(f"Heat transfer coefficient of LOx @ Impinging Lines: {round(heat_transfer_coeff_LOx*(1/m2_to_in2), 5)} W/(in^2*K)")
# Distribution Lines (angled sections); Redoing friction and Reynold's Number calculations with differing area
LOx_angled_dia_imperical = 0.092 # in
LOx_angled_dia = LOx_angled_dia_imperical / m_to_in # m
heat_transfer_coeff_LOx_angled_line = heat_transfer("LOx", LOx_angled_dia, num_LOx_orifices, "Angled Line")[0]
print(f"Heat transfer coefficient of LOx @ Angled Distribution Lines: {round(heat_transfer_coeff_LOx_angled_line*(1/m2_to_in2), 5)} W/(in^2*K)")


# Calculations for Kerosene

# Impinging lines
heat_transfer_coeff_Kerosene = heat_transfer("Kerosene", Kerosene_orifice_diameter, num_Kerosene_orifices, "Single Line")[0]
print(f"Heat transfer coefficient of Kerosene @ Impinging Lines: {round(heat_transfer_coeff_Kerosene*(1/m2_to_in2), 5)} W/(in^2*K)")
# Distribution Lines (larger drill section); Redoing friction and Reynold's Number calculations with differing area
num_distribution_orifices = 16 # unitless
Kerosene_large_dist_dia_imperical = .15748 # in
Kerosene_large_dist_dia = Kerosene_large_dist_dia_imperical / m_to_in # m
heat_transfer_coeff_Kerosene_large_dist = heat_transfer("Kerosene", Kerosene_large_dist_dia, num_distribution_orifices, "Large Distribution Line")[0]
print(f"Heat transfer coefficient of Kerosene @ Larger Distribution Lines: {round(heat_transfer_coeff_Kerosene_large_dist*(1/m2_to_in2), 5)} W/(in^2*K)")
# Distribution Lines (smaller drill section); Redoing friction and Reynold's Number calculations with differing area
Kerosene_small_dist_dia_imperical = .1378 # in
Kerosene_small_dist_dia = Kerosene_small_dist_dia_imperical / m_to_in # m
heat_transfer_coeff_Kerosene_small_dist = heat_transfer("Kerosene", Kerosene_small_dist_dia, num_distribution_orifices, "Small Distribution Line")[0]
print(f"Heat transfer coefficient of Kerosene @ Smaller Distribution Lines: {round(heat_transfer_coeff_Kerosene_small_dist*(1/m2_to_in2), 5)} W/(in^2*K)\n")

# Metric units
print(f"Heat transfer coefficient of LOx @ Impinging Lines: {round(heat_transfer_coeff_LOx, 5)} W/(m^2*K)")
print(f"Heat transfer coefficient of LOx @ Angled Distribution Lines: {round(heat_transfer_coeff_LOx_angled_line, 5)} W/(m^2*K)")
print(f"Heat transfer coefficient of Kerosene @ Impinging Lines: {round(heat_transfer_coeff_Kerosene, 5)} W/(m^2*K)")
print(f"Heat transfer coefficient of Kerosene @ Larger Distribution Lines: {round(heat_transfer_coeff_Kerosene_large_dist, 5)} W/(m^2*K)") 
print(f"Heat transfer coefficient of Kerosene @ Smaller Distribution Lines: {round(heat_transfer_coeff_Kerosene_small_dist, 5)} W/(m^2*K)\n")



'''
cea_obj.get_Chamber_Transport(Pc = chamber_pressure_imperial, MR = OF_ratio)
print(f"Viscosity of Gas Mixture: {round(output[0]*mP_to_Pa_s, 10)} Pa*s")
print(f"Prandtl number: {round(output[3], 5)}\n")
'''


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