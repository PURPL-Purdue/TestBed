import numpy as np
import math
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

# =====================================================
# === USER CONFIGURATION & CONSTANTS ==================
# =====================================================
OXIDIZER = "LOX"
BASE_FUEL = "ETHANOL"
CHAMBER_PRESSURE_PSI = 200.0   
OF_TARGET = 1.0                
AMBIENT_P_PSI = 14.7
LSTAR = 1.0922

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822

MILLIPOISE_TO_PAS = 0.0001        
MCAL_CM_K_S_TO_W_MK = 0.4184      
CAL_G_K_TO_J_KG_K = 4184.0        

# =====================================================
# === FIXED GEOMETRY INPUTS ===========================
# =====================================================
FIXED_THROAT_DIA_MM = 25.0       
FIXED_EXPANSION_RATIO = 3.5      # Ae / At
FIXED_CONTRACTION_RATIO = 4.0    # Ac / At

# =====================================================
# === ENTHALPY-BASED INJECTOR EFFICIENCY ==============
# =====================================================
injector_efficiency = 0.95
run_mode = input("Enter run mode ('ideal' or 'optimistic'): ").strip().lower()

# Standard liquid ethanol heat of formation is approx -66370 cal/mol
IDEAL_HEX = -66370.0 

FUEL_NAME = BASE_FUEL
if run_mode == "optimistic":
    FUEL_NAME = "EFF_ENTHALPY_ETHANOL"
    # Scale down the chemical enthalpy release to reflect combustion efficiency 
    scaled_enthalpy = IDEAL_HEX * injector_efficiency
    
    custom_fuel_card = f"""
    fuel {FUEL_NAME} C 2 H 6 O 1 wt%=100.00 h,cal={scaled_enthalpy:.1f} t(k)=298.15
    """
    add_new_fuel(FUEL_NAME, custom_fuel_card)

# =====================================================
# === PROPERTY EXTRACTION ENGINE (FIXED GEOMETRY) =====
# =====================================================
def compute_cea_properties(pc_psi, of, eps, contraction_ratio, dt_mm):
    Dt = dt_mm / 1000.0  
    At = (math.pi / 4.0) * (Dt ** 2)
    A_chamber = At * contraction_ratio
    D_chamber = math.sqrt(4.0 * A_chamber / math.pi)

    cea_fac = CEA_Obj(oxName=OXIDIZER, fuelName=FUEL_NAME, fac_CR=contraction_ratio)

    ch_temp, th_temp, _ = cea_fac.get_Temperatures(Pc=pc_psi, MR=of, eps=eps)
    ch_rho_raw, th_rho_raw, _ = cea_fac.get_Densities(Pc=pc_psi, MR=of, eps=eps)
    ch_sonic_raw, th_sonic_raw, _ = cea_fac.get_SonicVelocities(Pc=pc_psi, MR=of, eps=eps)

    ch_mach = cea_fac.get_Chamber_MachNumber(Pc=pc_psi, MR=of, fac_CR=contraction_ratio)
    th_mach = 1.0

    ch_vel_m_s = (ch_sonic_raw * ch_mach) * 0.3048
    ch_rho_m_s = ch_rho_raw * 16.0185
    th_vel_m_s = th_sonic_raw * th_mach * 0.3048 
    th_rho_m_s = th_rho_raw * 16.0185

    ch_trans = cea_fac.get_Chamber_Transport(Pc=pc_psi, MR=of)
    ch_cp = ch_trans[0] * CAL_G_K_TO_J_KG_K
    ch_mu = ch_trans[1] * MILLIPOISE_TO_PAS
    ch_k  = ch_trans[2] * MCAL_CM_K_S_TO_W_MK
    ch_pr = ch_trans[3]

    th_trans = cea_fac.get_Throat_Transport(Pc=pc_psi, MR=of)
    th_cp = th_trans[0] * CAL_G_K_TO_J_KG_K
    th_mu = th_trans[1] * MILLIPOISE_TO_PAS
    th_k  = th_trans[2] * MCAL_CM_K_S_TO_W_MK
    th_pr = th_trans[3]

    _, ch_gamma = cea_fac.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=of, eps=eps)
    _, th_gamma = cea_fac.get_Throat_MolWt_gamma(Pc=pc_psi, MR=of, eps=eps)

    cstar = cea_fac.get_Cstar(Pc=pc_psi, MR=of) * 0.3048

    return {
        'cstar': cstar, 'th_Dia': Dt, 'ch_Dia': D_chamber,
        'ch_temp': ch_temp, 'ch_density': ch_rho_m_s, 'ch_vel': ch_vel_m_s,
        'ch_cp': ch_cp, 'ch_mu': ch_mu, 'ch_k': ch_k, 'ch_pr': ch_pr, 'ch_gamma': ch_gamma,
        'th_temp': th_temp, 'th_density': th_rho_m_s, 'th_vel': th_vel_m_s,
        'th_cp': th_cp, 'th_mu': th_mu, 'th_k': th_k, 'th_pr': th_pr, 'th_gamma': th_gamma,
    }

# =====================================================
# === DIRECT THERMAL INPUT ASSIGNMENT =================
# =====================================================
cea_data = compute_cea_properties(
    CHAMBER_PRESSURE_PSI, 
    OF_TARGET, 
    FIXED_EXPANSION_RATIO, 
    FIXED_CONTRACTION_RATIO, 
    FIXED_THROAT_DIA_MM
)

th_GasTemp = cea_data['th_temp']
th_GasDensity = cea_data['th_density']
th_GasVelocity = cea_data['th_vel']
th_GasDynamicViscosity = cea_data['th_mu']
th_GasSpecificHeat = cea_data['th_cp']
th_GasThermalConductivity = cea_data['th_k']
th_GasPr = cea_data['th_pr']
th_Gamma = cea_data['th_gamma']
th_Dia = cea_data['th_Dia']

ch_GasTemp = cea_data['ch_temp']
ch_GasDensity = cea_data['ch_density']
ch_GasVelocity = cea_data['ch_vel']
ch_GasDynamicViscosity = cea_data['ch_mu']
ch_GasSpecificHeat = cea_data['ch_cp']
ch_GasThermalConductivity = cea_data['ch_k']
ch_GasPr = cea_data['ch_pr']
ch_Gamma = cea_data['ch_gamma']
ch_Dia = cea_data['ch_Dia']

# --- Input Coolant Properties (Water) ---
WaterTemp = 295                     
WaterDensity = 997.8                
WaterDynamicViscosity = 0.0009566   
WaterSpecificHeat = 4180.0          
WaterThermalConductivity = 0.44955  

# --- Engine & Channel Geometry ---
RadOfCurve = 0.0395             
ChannelCount = 54
PumpRate = 4500                 
THChannelDia = 0.00125           
CHChannelDia = 0.003            

wall_thickness_chamber_mm = 2.5 
wall_thickness_throat_mm = 1.25  

WallConductivity = 142          
Pc = CHAMBER_PRESSURE_PSI * PSI_TO_PA 
Cstar = cea_data['cstar']
YieldTemp = 650                 
watersattemp = 462.235          

th_WallThickness = wall_thickness_throat_mm / 1000.0   
ch_WallThickness = wall_thickness_chamber_mm / 1000.0  

WaterMassflow = PumpRate * WaterDensity * 0.0000010515 
WaterPr = WaterSpecificHeat * WaterDynamicViscosity / WaterThermalConductivity
Tolerance = 0.1 
omega = 0.6 

# =====================================================
# === NOMINAL RUN: THROAT STATION (BARTZ) =============
# =====================================================
GasTemp = th_GasTemp
GasDensity = th_GasDensity
GasVelocity = th_GasVelocity
GasDynamicViscosity = th_GasDynamicViscosity
GasSpecificHeat = th_GasSpecificHeat
GasThermalConductivity = th_GasThermalConductivity
GasPr = th_GasPr
Gamma = th_Gamma
LocalDia = th_Dia
Mach_local = 1.0
LocalChannelDia = THChannelDia

WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (LocalChannelDia/2)**2) * WaterDensity)
ReWater = WaterDensity * WaterVelocity * LocalChannelDia / WaterDynamicViscosity
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / LocalChannelDia

HotWallTemp = (GasTemp + WaterTemp) / 2.0
PrevHotWallTemp = 0.0
IterationCount = 0

while abs(HotWallTemp - PrevHotWallTemp) > Tolerance:
    PrevHotWallTemp = HotWallTemp
    sigma = ( (0.5 * HotWallTemp / GasTemp * (1 + (Gamma - 1) / 2 * Mach_local**2) + 0.5)**(0.8 - omega / 5.0) * (1 + (Gamma - 1) / 2 * Mach_local**2)**(omega/5.0) )**-1
    GasTransferCoeff = ( (0.026 / LocalDia**0.2) * (GasDynamicViscosity**0.2 * GasSpecificHeat / (GasPr**0.6)) * (Pc / Cstar)**0.8 * (LocalDia / RadOfCurve)**0.1 * sigma )
    th_flux = (GasTemp - WaterTemp) / (1.0 / GasTransferCoeff + th_WallThickness / WallConductivity + 1.0 / WaterTransferCoeff)
    HotWallTemp = GasTemp - (th_flux / GasTransferCoeff)
    IterationCount += 1
    if IterationCount > 1000: break

th_HotWallTemp = HotWallTemp
th_ColdWallTemp = WaterTemp + (th_flux / WaterTransferCoeff)
th_WaterVelocity = WaterVelocity

# =====================================================
# === NOMINAL RUN: CHAMBER STATION (DITTUS-BOELTER) ===
# =====================================================
GasTemp = ch_GasTemp
GasDensity = ch_GasDensity
GasVelocity = ch_GasVelocity
GasDynamicViscosity = ch_GasDynamicViscosity
GasSpecificHeat = ch_GasSpecificHeat
GasThermalConductivity = ch_GasThermalConductivity
GasPr = ch_GasPr
LocalDia = ch_Dia
LocalChannelDia = CHChannelDia

WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (LocalChannelDia/2)**2) * WaterDensity)
ReWater = WaterDensity * WaterVelocity * LocalChannelDia / WaterDynamicViscosity
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / LocalChannelDia

ReGas = GasDensity * GasVelocity * LocalDia / GasDynamicViscosity
GasNusselts = 0.023 * ReGas**0.8 * GasPr**0.3
GasTransferCoeff = GasNusselts * GasThermalConductivity / LocalDia

ch_flux = (GasTemp - WaterTemp) / (1.0 / GasTransferCoeff + ch_WallThickness / WallConductivity + 1.0 / WaterTransferCoeff)
ch_HotWallTemp = GasTemp - (ch_flux / GasTransferCoeff)
ch_ColdWallTemp = WaterTemp + (ch_flux / WaterTransferCoeff)
ch_WaterVelocity = WaterVelocity

print(f"\n--- Nominal Converged States ({run_mode.upper()} mode) ---")
print(f"Throat  [Bartz] (ChDia={THChannelDia*1000:.1f}mm | Wall={wall_thickness_throat_mm:.2f}mm) -> Coolant Vel: {th_WaterVelocity:.2f} m/s | Hot Wall: {th_HotWallTemp:.2f} K | Cold Wall: {th_ColdWallTemp:.2f} K | Flux: {th_flux/1e6:.3f} MW/m^2")
print(f"Chamber [D-B]   (ChDia={CHChannelDia*1000:.1f}mm | Wall={wall_thickness_chamber_mm:.2f}mm) -> Coolant Vel: {ch_WaterVelocity:.2f} m/s | Hot Wall: {ch_HotWallTemp:.2f} K | Cold Wall: {ch_ColdWallTemp:.2f} K | Flux: {ch_flux/1e6:.3f} MW/m^2\n")

# =====================================================
# === PLOTS ===========================================
# =====================================================
ThicknessSweep = np.linspace(0.0005, 0.003, 100)

# --- Loop 1A: Process Throat Curve (Bartz Loop) ---
GasTemp, GasDynamicViscosity, GasSpecificHeat, GasPr, Gamma, LocalDia, Mach_local, LocalChannelDia = th_GasTemp, th_GasDynamicViscosity, th_GasSpecificHeat, th_GasPr, th_Gamma, th_Dia, 1.0, THChannelDia
WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (LocalChannelDia/2)**2) * WaterDensity)
ReWater = WaterDensity * WaterVelocity * LocalChannelDia / WaterDynamicViscosity
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / LocalChannelDia

th_SweepHW, th_SweepCW = [], []
for CurrentThickness in ThicknessSweep:
    TempHW, PrevTempHW, LoopCount = (GasTemp + WaterTemp) / 2.0, 0.0, 0
    while abs(TempHW - PrevTempHW) > Tolerance:
        PrevTempHW = TempHW
        sigma_sweep = ( (0.5 * TempHW / GasTemp * (1 + (Gamma - 1) / 2 * Mach_local**2) + 0.5)**(0.8 - omega / 5.0) * (1 + (Gamma - 1) / 2 * Mach_local**2)**(omega/5.0) )**-1
        GasCoeff_sweep = ( (0.026 / LocalDia**0.2) * (GasDynamicViscosity**0.2 * GasSpecificHeat / (GasPr**0.6)) * (Pc / Cstar)**0.8 * (LocalDia / RadOfCurve)**0.1 * sigma_sweep )
        flux_sweep = (GasTemp - WaterTemp) / (1.0 / GasCoeff_sweep + CurrentThickness / WallConductivity + 1.0 / WaterTransferCoeff)
        TempHW = GasTemp - (flux_sweep / GasCoeff_sweep)
        if LoopCount > 1000: break
    th_SweepHW.append(TempHW)
    th_SweepCW.append(WaterTemp + (flux_sweep / WaterTransferCoeff))

# --- Loop 1B: Process Chamber Curve (Dittus-Boelter Direct) ---
GasTemp, GasDensity, GasVelocity, GasDynamicViscosity, GasSpecificHeat, GasThermalConductivity, GasPr, LocalDia, LocalChannelDia = ch_GasTemp, ch_GasDensity, ch_GasVelocity, ch_GasDynamicViscosity, ch_GasSpecificHeat, ch_GasThermalConductivity, ch_GasPr, ch_Dia, CHChannelDia
WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (LocalChannelDia/2)**2) * WaterDensity)
ReWater = WaterDensity * WaterVelocity * LocalChannelDia / WaterDynamicViscosity
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / LocalChannelDia

ReGas = GasDensity * GasVelocity * LocalDia / GasDynamicViscosity
GasNusselts = 0.023 * ReGas**0.8 * GasPr**0.3
GasCoeff_sweep = GasNusselts * GasThermalConductivity / LocalDia

ch_SweepHW, ch_SweepCW = [], []
for CurrentThickness in ThicknessSweep:
    flux_sweep = (GasTemp - WaterTemp) / (1.0 / GasCoeff_sweep + CurrentThickness / WallConductivity + 1.0 / WaterTransferCoeff)
    ch_SweepHW.append(GasTemp - (flux_sweep / GasCoeff_sweep))
    ch_SweepCW.append(WaterTemp + (flux_sweep / WaterTransferCoeff))

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ax1.plot(ThicknessSweep * 1000, th_SweepHW, label='Throat (Bartz) - Hot Wall', color='crimson', linewidth=3)
ax1.plot(ThicknessSweep * 1000, ch_SweepHW, label='Chamber (D-B) - Hot Wall', color='darkred', linewidth=3)
ax1.axhline(y=YieldTemp, color='purple', linestyle=':', linewidth=2.5, label=f'Material Yield ({YieldTemp} K)')
ax1.plot(wall_thickness_throat_mm, th_HotWallTemp, marker='o', markersize=10, color='black')
ax1.plot(wall_thickness_chamber_mm, ch_HotWallTemp, marker='s', markersize=10, color='blue')
ax1.set_title('Thickness vs. Hot Wall Temp', fontsize=14, fontweight='bold')
ax1.grid(True, linestyle='--')
ax1.legend()

ax2.plot(ThicknessSweep * 1000, th_SweepCW, label='Throat (Bartz) - Cold Wall', color='royalblue', linewidth=3)
ax2.plot(ThicknessSweep * 1000, ch_SweepCW, label='Chamber (D-B) - Cold Wall', color='navy', linewidth=3)
ax2.axhline(y=watersattemp, color='darkorange', linestyle='--', linewidth=2.5, label='Water Saturation Temp')
ax2.plot(wall_thickness_throat_mm, th_ColdWallTemp, marker='o', markersize=10, color='black')
ax2.plot(wall_thickness_chamber_mm, ch_ColdWallTemp, marker='s', markersize=10, color='blue')
ax2.set_title('Thickness vs. Cold Wall Temp', fontsize=14, fontweight='bold')
ax2.grid(True, linestyle='--')
ax2.legend()
plt.tight_layout()
plt.show()