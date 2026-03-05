import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from rocketcea.cea_obj import CEA_Obj
import CoolProp.CoolProp as CP
import math


import math


def calculate_convection_coeff(ox, fuel, of_ratio, p_chamber_pa, D_inner, v_gas, eps):
    """
    Uses RocketCEA to get gas properties and calculates convection coefficient.
    p_chamber_pa: Chamber pressure in Pascals
    """
    #Velocity of gas at chamber
    #v_gas = cea.get_SonicVelocities(Pc=p_chamber_pa / 6894.76, MR=of_ratio, eps=eps)[1] * 0.3048 # Approximate velocity at chamber (half of sonic velocity)
    
    # Convert Pa to PSI for RocketCEA
    pc_psi = p_chamber_pa / 6894.76
    
    # Create CEA object
    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    
    # 1. Get Combustion Properties
    # Tc: Flame Temp, M: Mol Weight, gamma: Ratio of specific heats
    # We get these at the chamber (cp)
    I_sp = cea.get_Isp(Pc=pc_psi, MR=of_ratio, eps=eps, frozen=1, frozenAtThroat=1) # Placeholder for structure
    # Alternatively, get more detailed properties:
    T_k = 5/9 *cea.get_Temperatures(Pc=pc_psi, MR=of_ratio, eps=eps, frozen=1)[1]
    
    # 2. Get Transport Properties
    # (Cp in kJ/kg-K, Viscosity in millipoise, Pr)
    cp_kj, viscosity_mp, pr, k_g = cea.get_Chamber_Transport(Pc=pc_psi, MR=of_ratio, eps=eps)
    
    # Convert units to SI
    cp_g = cp_kj * 1000.0          # J/kg-K
    mu_g = viscosity_mp * 0.0001   # Pa-s (millipoise to kg/m-s)
    #k_g = (pr != 0) and (cp_g * mu_g / pr) or 0.06 # Thermal conductivity
    
    # 3. Density (Ideal Gas Law)
    # R_specific = R_universal / M
    R_univ = 8314.46 # J/(kmol-K)
    M_val = cea.get_IvacCstrTc_ChmMwGam(Pc=pc_psi, MR=of_ratio, eps=1.0)[3] # kg/kmol
    rho_g = p_chamber_pa / ((R_univ / M_val) * T_k)
    
    # 4. Reynolds Number & Gnielinski Correlation
    re = (rho_g * v_gas * D_inner) / mu_g
    
    if re < 2300:
        nu = 3.66
    else:
        f = (0.79 * np.log(re) - 1.64)**(-2)
        nu = ((f/8) * (re - 1000) * pr) / (1 + 12.7 * (f/8)**0.5 * (pr**(2/3) - 1))
    
    h = (nu * k_g) / D_inner
    return h, T_k, re

# # =================== PARAMETERS ====================
# ri = 0.5 * 0.0254            # m
# ro =  0.6* 0.0254            # m
# length = 0.02 * 0.0254       # m
# thickness = ro - ri

# # ---------- Combustion Setup ----------
# OX = "LOX"
# FUEL = "RP1"
# OF_RATIO = 2.25
# P_CHAMBER = 1.724e+6        # 250 psi in Pa
# V_GAS = 50.0                 # m/s

# # ---------- Material properties ----------
# k_metal = 167      # W/(m·K)  
# c_metal = 896       # J/(kg·K)   
# rho_metal = 2700
# alpha = k_metal / (rho_metal * c_metal)
# T_melt = 582.0 + 273.15  # Celsius to Kelvin

# # ---------- Gas Properties (2500C Air) ----------
# #T_gas = 2500.0 + 273.15  # Celsius to Kelvin
# h_conv = calculate_convection_coeff(ox=OX, fuel=FUEL, of_ratio=OF_RATIO, p_chamber_pa=P_CHAMBER, D_inner=2*ri, v_gas=V_GAS)[0]    # Convection coefficient (W/m^2*K) - calculated from Re/Nu
# T_gas = (5/9) *calculate_convection_coeff(ox=OX, fuel=FUEL, of_ratio=OF_RATIO, p_chamber_pa=P_CHAMBER, D_inner=2*ri, v_gas=V_GAS)[1]    # Gas temperature (K) - from CEA
# T_initial = 25.0 + 273.15  # Initial metal temp (Celsius to Kelvin)

# ==================== TRANSIENT SOLVER ====================
def heat_equation(t,T, r, dr, n_nodes, alpha, h_conv, T_gas, k_metal):
    dTdt = np.zeros_like(T)
    
    for i in range(n_nodes):
        if i == 0: # INNER WALL (Convection Boundary)
            # Energy balance: h(Tgas - T[0]) + k(dT/dr) = 0
            # Using second-order discretization for the boundary:
            dT_dr = (T[1] - T[0]) / dr
            heat_from_gas = h_conv * (T_gas - T[i])
            # This is a simplified surface node balance:
            dTdt[i] = (alpha * (2/dr**2) * (T[i+1] - T[i] + (dr*h_conv/k_metal)*(T_gas - T[i])))
            
        elif i == n_nodes - 1: # OUTER WALL (Insulated for this example)
            dTdt[i] = alpha * (2/dr**2) * (T[i-1] - T[i])
            
        else: # INTERIOR NODES
            # d2T/dr2 + (1/r)dT/dr
            term1 = (T[i+1] - 2*T[i] + T[i-1]) / dr**2
            term2 = (1/r[i]) * (T[i+1] - T[i-1]) / (2*dr)
            dTdt[i] = alpha * (term1 + term2)
            
    return dTdt

# # ==================== SOLVE ====================
# n_nodes = 50
# r_nodes = np.linspace(ri, ro, n_nodes)
# dr = r_nodes[1] - r_nodes[0]
# T0 = np.full(n_nodes, T_initial)
# t_final = 10 # 5 minutes

# sol = solve_ivp(heat_equation, [0, t_final], T0, args=(r_nodes, dr, n_nodes, alpha, h_conv, T_gas, k_metal), 
#                 method='BDF', t_eval=np.linspace(0, t_final, 500))

# # ==================== RESULTS & PLOTTING ====================
# T_res = sol.y
# times = sol.t


# inner_wall_temp = T_res[5, :]

# # Find melting time
# melt_indices = np.where(inner_wall_temp >= T_melt)[0]
# if len(melt_indices) > 0:
#     t_melt = times[melt_indices[0]]
#     print(f"CRITICAL: Inner wall melts at {t_melt:.3f} seconds.")
# else:
#     t_melt = -1
#     print("Wall did not melt within the time frame.")


# plt.figure(figsize=(10, 6))
# plt.plot(times, inner_wall_temp, label='Inner Wall (at Gas Interface)')
# plt.plot(times, T_res[-1, :], label='Outer Wall')
# plt.axhline(T_melt, color='red', linestyle='--', label='Melting Point')
# plt.xlabel('Time (s)')
# plt.ylabel('Temperature (K)')
# plt.title('Cylinder Temperature Response')
# plt.legend()
# plt.grid(True)
# plt.show()
