##



##
import numpy as np
import math

GasTemp = 2350.0        #K
GasDensity = 
GasVelocity = 
GasDynamicViscosity = 
GasSpecificHeat = 
GasThermalConductivity = 

WaterTemp = 295         #K
WaterDensity = 997.8       #kg/m^3
#WaterVelocity = 
WaterDynamicViscosity = 0.0009566  #Pa-s
WaterSpecificHeat = 4.18 * 1000   #J/kg-K
WaterThermalConductivity = 449.55 / 1000   #W/m-K


ChannelCount = 28
PumpRate = 3000 #gph
ChannelDia = 0.01 #m
ThroatDia = 0.0527 #m
WaterMassflow = PumpRate * WaterDensity * 0.0000010515 #kg/s
WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (ChannelDia/2)**2) * WaterDensity) #m/s
WallThickness = 0.005 #m
WallConductivity = 177 #W/m-K

#Calculations
################################################################################

#Prandtl Numbers
WaterPr = WaterSpecificHeat * WaterDynamicViscosity / WaterThermalConductivity
GasPr = GasSpecificHeat * GasDynamicViscosity / GasThermalConductivity

#Reynolds Numbers
ReWater = WaterDensity * WaterVelocity * ChannelDia / WaterDynamicViscosity
ReGas = GasDensity * GasVelocity * ThroatDia / GasDynamicViscosity


#Bartz
GasNusselts = 0.026 * ReGas**0.8 * GasPr**0.4

#Dittus Boelter
WaterNusselts = .023 * ReWater**0.8 * WaterPr**0.4 


GasTransferCoeff = GasNusselts * GasThermalConductivity / ThroatDia
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / ChannelDia



flux = (GasTemp - WaterTemp) / (1 / GasTransferCoeff + WallThickness / WallConductivity + 1 / WaterTransferCoeff)