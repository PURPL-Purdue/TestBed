import numpy as np
import matplotlib.pyplot as plt


########################################
#geometric properties
########################################

throatID = 
throatThick = 
throatChanDia = 
chamberID = 
chamberThick = 
chamberChanDia = 
exitID = 

eps = exitID / throatID
ctr = chamberID / throatID
roc = 
engLen = 
channelCount = 54


########################################
#System parameters
########################################
resolution = 300
step = engLen / resolution
PumpRate = 4500 #GPH
ox = "LOX"
fuel = "ETHANOL"
thrust = 900 #lbf
Pc = 200 #psi
of = 1.0
ambient = 14.7 #psi


########################################
#Conversions
########################################

gphtoMdotWater = 



########################################
#Fluid and Material Property Equations
########################################

#7175 yield as function of time

#water dynamic visco as a funtion of T

#7175 conductivity as funtion of T



########################################
#Local Fluid Elements
########################################

#local gas Element
gasElement = []

#local fluid element
waterElement = []




########################################
#CEA for gas properties
########################################
def fac_CEA(ox, fuel, of, ctr, eps, pc, location):
    




    return {
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 
        'mdot': mdot, 'localtemp': localtemp, 'localdensity': localrho, 
        'localvel': localvel, 'localcp': localcp, 'localmu': localmu, 
        'local_k': localK, 'localPR': localPR, 'ch_gamma': ch_gamma,
    }


########################################
#1Dlocal Thermal
########################################
def localThermal():
    
    
    
    return{
    }



########################################
#Displays
########################################
def graphsNdata():
    #graph of engine contour with pressure display

    #graph of engine contour with hotwall and coldwall T

    #graph of engine contour with fluid temperature

    #graph of engine contour with heatflux through channels
    
    return{
    }


########################################
#Main
########################################
print("Starting analysis...")
location = 0
x = []
GasTemp = []
ColdwallT = []
HotwallT = []
heatflux = []



while location < engLen:
    fluidprops = localwater()
    gasprops = fac_CEA()
    localthermal = localThermal(fluidprops, gasprops)
    x.append(location)
    GasTemp.append(T)
    ColdwallT.append(Tcw)
    HotwallT.append(Thw)
    heatflux.append(q)
    location += step
print("")