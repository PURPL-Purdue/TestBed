#

#



import numpy as np

# --- Input Gas Properties (Throat) ---
GasTemp = 2350.0                # K
GasDensity = 0.808              # kg/m^3
GasVelocity = 3555.9            # m/s
GasDynamicViscosity = 7.7529e-5 # Pa-s
GasSpecificHeat = 2526.1        # J/kg-K
GasThermalConductivity = 0.3635 # W/m-K
GasPr = 0.5384                  # From CEA
Gamma = 1.2243                  # From CEA

# --- Input Coolant Properties (Water) ---
WaterTemp = 295                     # K
WaterDensity = 997.8                # kg/m^3
WaterDynamicViscosity = 0.0009566   # Pa-s
WaterSpecificHeat = 4180.0          # J/kg-K
WaterThermalConductivity = 0.44955  # W/m-K

# --- Engine & Channel Geometry ---
RadOfCurve = 0.0395             # m (Nozzle throat radius of curvature)
ChannelCount = 40
PumpRate = 4500                 # gph
ChannelDia = 0.002              # m
ThroatDia = 0.0527              # m
WallThickness = 0.003          # m
WallConductivity = 177          # W/m-K (7175)
Pc = 1.379e+6                   # Pa (200 psi)
Cstar = 1577.6                  # m/s (from CEA)
WaterFeed = 2.413e+6            # Pa (350 psi)
AmbPressure = 101352.9          # Pa (14.7 psi)


WaterMassflow = PumpRate * WaterDensity * 0.0000010515 
#WaterMassflow = 
#WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (ChannelDia/2)**2) * WaterDensity) #Use for pump-fed system
WaterVelocity = (2 * (WaterFeed - AmbPressure) / WaterDensity)**(0.5)                        #Use for pressure-fed system

ReWater = WaterDensity * WaterVelocity * ChannelDia / WaterDynamicViscosity
WaterPr = WaterSpecificHeat * WaterDynamicViscosity / WaterThermalConductivity
ReGas = GasDensity * GasVelocity * ThroatDia / GasDynamicViscosity

#Dittus-Boelter
# n=0.4 for heating the fluid
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / ChannelDia

#BARTZ ITERATION LOOP
# Initial guess for wall temperature (average of gas and water)
HotWallTemp = (GasTemp + WaterTemp) / 2.0
PrevHotWallTemp = 0.0
Tolerance = 0.1 # Kelvin
IterationCount = 0

# Bartz specific constant
omega = 0.6 

while abs(HotWallTemp - PrevHotWallTemp) > Tolerance:
    PrevHotWallTemp = HotWallTemp
    
    sigma = ( (0.5 * HotWallTemp / GasTemp * (1 + (Gamma - 1) / 2 * 1**2) + 0.5)**(0.8 - omega / 5.0) *  (1 + (Gamma - 1) / 2 * 1**2)**(omega/5.0) )**-1
    

    GasTransferCoeff = ( (0.026 / ThroatDia**0.2) * (GasDynamicViscosity**0.2 * GasSpecificHeat / (GasPr**0.6)) * (Pc / Cstar)**0.8 * (ThroatDia / RadOfCurve)**0.1 * sigma )
    
    flux = (GasTemp - WaterTemp) / (1.0 / GasTransferCoeff + WallThickness / WallConductivity + 1.0 / WaterTransferCoeff)
    HotWallTemp = GasTemp - (flux / GasTransferCoeff)
    
    IterationCount += 1
    if IterationCount > 1000: 
        break

ColdWallTemp = WaterTemp + (flux / WaterTransferCoeff)






# =====================================================
# === OUTPUTS =========================================
# =====================================================

print(f"--- Iteration Converged in {IterationCount} steps ---")
print(f"Hot Wall Temp: {HotWallTemp:.2f} K")
print(f"Cold Wall Temp: {ColdWallTemp:.2f} K")
print(f"Heat Flux: {flux:.2f} W/m^2")
print(f"Water Velocity: {WaterVelocity:.2f} m/s")