from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

# Conversion Ratio
rankineToKelvin = 0.555556

# Chamber Pressure (psi)
#p_c = 1000 # Test 300 to 3000 with 300 increments
while True:
    p_c = input("Enter chamber pressure in psi: ")
    try:
        p_c = float(p_c)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

# OF Inputs
# OF_min = 0.1                       # minimum O/F ratio to test
# OF_max = 20.1                      # maximum O/F ratio to test
# OF_inc = 0.1                         # increment for test
while True:
    OF_min = input("Enter minimum O/F ratio to test: ")
    try:
        OF_min = float(OF_min)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    OF_max = input("Enter maximum O/F ratio to test. (Note: value used in code will be max + increment): ")
    try:
        OF_max = float(OF_max)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    OF_inc = input("Enter increment for test: ")
    try:
        OF_inc = float(OF_inc)  # convert to int
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

OF_max = OF_max + OF_inc  # Add increment to max to include max in test

# O/F ratios to test
OFs = np.arange(OF_min, OF_max, OF_inc)  # Test 0.1 to 20.0 with 0.1 increment

# In order to make graphs, we will pass in the OFs as the x axis, and 
# use a function to convert the OF input to a temperature output
# CEA Rocket Problem functions taking OF input
def RP1LOxOFToTemp(of):
    cea = CEA_Obj(oxName = "LOX", fuelName = "RP-1")
    return cea.get_Tcomb(p_c, of) * rankineToKelvin

# Chamber Temperatures
# We run the function for value in OFs and put it into a new array of
# temperatures
RP1t_c = [RP1LOxOFToTemp(of) for of in OFs]

#---- INTERPOLATION TO FIND O/F FOR TARGET TEMPERATURE ----#
# Suppose you want the O/F value for a target temperature
while True:
    T_target = input("Enter target chamber temperature to find O/F at Temp in K: ")
    try:
        T_target = float(T_target)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

