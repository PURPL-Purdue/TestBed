from rocketcea.cea_obj import CEA_Obj
import numpy as np

# Constants
g = 9.81
#mix_eff = 0.9
while True:
    mix_eff = input("Enter mix efficiency as a decimal: ")
    try:
        mix_eff = float(mix_eff)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

p_amb = 14.7
#Lstar = 1 #You can change this to anything between 0.8 and 3.0, I just picked 1.1 so that it wouldnt be too skinny
#1.02 - 1.27 m is a common range for L*
while True:
    Lstar = input("Enter L* value: ")
    try:
        Lstar = float(Lstar)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

# Conversion Ratios
feetToMeters = 0.3048
psiToPa = 6894.76
lbftoN = 4.448
inchestometers = 0.0254


# ________Setpoints_________

# Running Stoich RP1/LOx:
#o_f = 2.0 #Run the OF/Temp graph to approximate a good OF value using your chosen pressure
#p_c = 1000 # psia, You can change this as well, I chose 1000 since I thought a chamber pressure between 500 and 1500 was good
#thrust = 5000 * lbftoN # Our thrust converted to Newtons
while True:
    o_f = input("Enter O/F value: ")
    try:
        o_f = float(o_f)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    p_c = input("Enter chamber pressure in psi: ")
    try:
        p_c = float(p_c)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    thrust = input("Enter thrust in lbf: ")
    try:
        thrust = float(thrust)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

thrust = thrust * lbftoN # Our thrust converted to Newtons

# Function to translate of, pc, thrust to a throat diameter
def RP1LOx_throat_sizing_function(of, pc, F):
    # Create a CEA Rocket Problem
    cea = CEA_Obj(oxName = "LOX", fuelName = "RP1")

    # Get Cstar
    cstar = cea.get_Cstar(Pc=pc, MR=of) * feetToMeters

    # Use CEA to get exhaust velocity
    # Can't use get_Throat_Isp() because it returns vacuum isp, not ambient
    # Instead, we will use Coefficient of Thrust (CF)
    # v_e (m/s) = c* x CF

    CF = cea.get_PambCf(Pamb=p_amb, Pc=pc, MR=of, eps=8)[0]

    v_e = cstar * CF * mix_eff  # Apply mixture efficiency
    isp = v_e / g
    print(f"ISP(s): {isp}")
    print(f"Exhaust Velocity (m/s): {v_e}")
    #Empty space for formatting
    print()

    # Calculate total mass flow rate
    mdot = F/v_e
    print(f"Total Mass Flow Rate (kg/s):{mdot}")
    print()

    # Size the throat
    # mdot = pc * A_t/cstar => mdot * cstar/pc = A_t
    # Need to convert pc to pascals for this equation
    pc_pa = pc * psiToPa
    A_t_meters = mdot * cstar/pc_pa
    print(f"Throat Area (m^2): {A_t_meters}")
    # Convert to diameter in inches
    D_t_meters = np.sqrt( 4* A_t_meters/np.pi)


    D_t_inches = D_t_meters / inchestometers 
    print(f"Throat Diameter (m): {D_t_meters}")

    V_c = Lstar * A_t_meters #Formula for Lstar rearranged for Volume of the chamber
    D_c_meters = (2*D_t_meters) #Throat diameter * 2
    A_c = (np.pi * D_c_meters**2 )/4 #Cross-sectional Area
    L_c = V_c / A_c #Chamber Length in meters

    print(f"Chamber Length (m) : {L_c}")
    print(f"Chamber Diameter (m): {D_c_meters}")



def main():
    RP1LOx_throat_sizing_function(o_f, p_c, thrust)

if __name__ == "__main__":
    main()