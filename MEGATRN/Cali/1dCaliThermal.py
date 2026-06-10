import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# === RENAMED STATION PROPERTIES =====================
# =====================================================

# --- Input Gas Properties (Throat) ---
th_GasTemp = 3798.41                # K
th_GasDensity = 0.7970              # kg/m^3
th_GasVelocity = 1083.82            # m/s
th_GasDynamicViscosity = 7.1822e-05 # Pa-s
th_GasSpecificHeat = 2407.9         # J/kg-K
th_GasThermalConductivity = 0.3046  # W/m-K
th_GasPr = 0.5678                   # From CEA
th_Gamma = 1.2332                   # From CEA
th_Dia = 0.0527                     # m

# --- Input Gas Properties (Chamber) ---
ch_GasTemp = 4221.04                # K
ch_GasDensity = 1.2671              # kg/m^3
ch_GasVelocity = 170.43             # m/s
ch_GasDynamicViscosity = 7.7391e-05 # Pa-s
ch_GasSpecificHeat = 2521.8         # J/kg-K
ch_GasThermalConductivity = 0.3621  # W/m-K
ch_GasPr = 0.5389                   # From CEA
ch_Gamma = 1.2244                   # From CEA
ch_Dia = 0.107                      # m 

# --- Input Coolant Properties (Water) ---
WaterTemp = 295                     # K
WaterDensity = 997.8                # kg/m^3
WaterDynamicViscosity = 0.0009566   # Pa-s
WaterSpecificHeat = 4180.0          # J/kg-K
WaterThermalConductivity = 0.44955  # W/m-K

# --- Engine & Channel Geometry ---
RadOfCurve = 0.0395             # m (Nozzle throat radius of curvature)
ChannelCount = 54
PumpRate = 4500                 # gph
THChannelDia = 0.00125           # m (Throat coolant channel diameter)
CHChannelDia = 0.003            # m (Chamber coolant channel diameter)
WallThickness = 0.00125          # m
WallConductivity = 177          # W/m-K (7175 Aluminum)
Pc = 1.379e+6                   # Pa (200 psi)
Cstar = 2113.8                  # m/s (from CEA)
YieldTemp = 561               # K (Material structural yielding point parameter)
watersattemp = 462.235          # K (Water saturation temperature at 200 psi)


#pr/t = 7830.18 yield required with 1.3x SF 73 KSI is tensile yield for 7075, so minimum yield strength percent has to be 10.7 % now check graph for temp, T = 620F/600K at 30min of exposure (1.75mm wall)
#pr/t = 11010.2 yield required with 1.3x SF 73 KSI is tensile yield for 7075, so minimum yield strength percent has to be 15.1 % now check graph for temp, T = 550F/561K at 30min of exposure (1.25mm wall)

# Global Constants & Iteration Setup
WaterMassflow = PumpRate * WaterDensity * 0.0000010515 
WaterPr = WaterSpecificHeat * WaterDynamicViscosity / WaterThermalConductivity
Tolerance = 0.1 # Kelvin
omega = 0.6 

# =====================================================
# === NOMINAL RUN: THROAT STATION (BARTZ) =============
# =====================================================
# Assign temporary variables for Throat
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

# Station Specific Coolant Math
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
    th_flux = (GasTemp - WaterTemp) / (1.0 / GasTransferCoeff + WallThickness / WallConductivity + 1.0 / WaterTransferCoeff)
    HotWallTemp = GasTemp - (th_flux / GasTransferCoeff)
    IterationCount += 1
    if IterationCount > 1000: break

th_HotWallTemp = HotWallTemp
th_ColdWallTemp = WaterTemp + (th_flux / WaterTransferCoeff)
th_WaterVelocity = WaterVelocity

# =====================================================
# === NOMINAL RUN: CHAMBER STATION (DITTUS-BOELTER) ===
# =====================================================
# Assign temporary variables for Chamber
GasTemp = ch_GasTemp
GasDensity = ch_GasDensity
GasVelocity = ch_GasVelocity
GasDynamicViscosity = ch_GasDynamicViscosity
GasSpecificHeat = ch_GasSpecificHeat
GasThermalConductivity = ch_GasThermalConductivity
GasPr = ch_GasPr
LocalDia = ch_Dia
LocalChannelDia = CHChannelDia

# Station Specific Coolant Math
WaterVelocity = WaterMassflow / (ChannelCount * (np.pi * (LocalChannelDia/2)**2) * WaterDensity)
ReWater = WaterDensity * WaterVelocity * LocalChannelDia / WaterDynamicViscosity
WaterNusselts = 0.023 * ReWater**0.8 * WaterPr**0.4 
WaterTransferCoeff = WaterNusselts * WaterThermalConductivity / LocalChannelDia

# Core Gas-Side Convective Math via Dittus-Boelter (n=0.3 for gas cooling)
ReGas = GasDensity * GasVelocity * LocalDia / GasDynamicViscosity
GasNusselts = 0.023 * ReGas**0.8 * GasPr**0.3
GasTransferCoeff = GasNusselts * GasThermalConductivity / LocalDia

# Direct formulation (No iteration loop required for static Dittus-Boelter properties)
ch_flux = (GasTemp - WaterTemp) / (1.0 / GasTransferCoeff + WallThickness / WallConductivity + 1.0 / WaterTransferCoeff)
ch_HotWallTemp = GasTemp - (ch_flux / GasTransferCoeff)
ch_ColdWallTemp = WaterTemp + (ch_flux / WaterTransferCoeff)
ch_WaterVelocity = WaterVelocity

print(f"--- Nominal Converged States ---")
print(f"Throat  [Bartz] (ChDia={THChannelDia*1000:.1f}mm) -> Coolant Vel: {th_WaterVelocity:.2f} m/s | Hot Wall: {th_HotWallTemp:.2f} K | Cold Wall: {th_ColdWallTemp:.2f} K | Flux: {th_flux/1e6:.3f} MW/m^2")
print(f"Chamber [D-B]   (ChDia={CHChannelDia*1000:.1f}mm) -> Coolant Vel: {ch_WaterVelocity:.2f} m/s | Hot Wall: {ch_HotWallTemp:.2f} K | Cold Wall: {ch_ColdWallTemp:.2f} K | Flux: {ch_flux/1e6:.3f} MW/m^2\n")


# =====================================================
# === PLOTS ===========================================
# =====================================================

# --- Data Generation Sweeps ---
ThicknessSweep = np.linspace(0.0005, 0.005, 100)

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


# --- Figure 1: Wall Thickness vs. Hot Wall Temperatures ---
plt.figure(1, figsize=(10, 5.5))
plt.plot(ThicknessSweep * 1000, th_SweepHW, label='Throat (Bartz) - Hot Wall', color='crimson', linewidth=3)
plt.plot(ThicknessSweep * 1000, ch_SweepHW, label='Chamber (D-B) - Hot Wall', color='darkred', linewidth=3)
plt.axhline(y=YieldTemp, color='purple', linestyle=':', linewidth=2.2, label=f'Material Yield Temperature ({YieldTemp:.1f} K)')
#point mapping
plt.plot(WallThickness * 1000, th_HotWallTemp, marker='o', markersize=8, color='black', label=f'Current Wall Thickness ({WallThickness * 1000:.2f} mm)')
plt.title('Wall Thickness vs. Structural Hot Wall Temperatures', fontsize=12, fontweight='bold')
plt.xlabel('Wall Thickness (mm)', fontsize=10)
plt.ylabel('Hot Wall Temperature (K)', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left')
plt.tight_layout()


# --- Figure 2: Wall Thickness vs. Cold Wall Temperatures ---
plt.figure(2, figsize=(10, 5.5))
plt.plot(ThicknessSweep * 1000, th_SweepCW, label='Throat (Bartz) - Cold Wall', color='royalblue', linewidth=2)
plt.plot(ThicknessSweep * 1000, ch_SweepCW, label='Chamber (D-B) - Cold Wall', color='navy', linewidth=2)
plt.axhline(y=watersattemp, color='darkorange', linestyle='--', linewidth=2.2, label='Water Saturation Temp (200 psi)')
plt.title('Wall Thickness vs. Fluid Interface Cold Wall Temperatures', fontsize=12, fontweight='bold')
plt.xlabel('Wall Thickness (mm)', fontsize=10)
plt.ylabel('Cold Wall Temperature (K)', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper right')  # Moved to top right to clear the lines
plt.ylim(375, 490)
plt.tight_layout()


# --- Figure 3: Channel Diameter Sweep vs. Water Velocity ---
DiaSweep = np.linspace(0.001, 0.006, 100) 
VelocitySweep = WaterMassflow / (ChannelCount * (np.pi * (DiaSweep/2)**2) * WaterDensity)

plt.figure(3, figsize=(9, 5))
plt.plot(DiaSweep * 1000, VelocitySweep, color='teal', linewidth=3, label='Velocity Curve')
plt.axvline(x=THChannelDia*1000, color='crimson', linestyle=':', linewidth=2.2, label=f'Throat Design Channel ({THChannelDia*1000:.1f} mm)')
plt.axvline(x=CHChannelDia*1000, color='darkred', linestyle='--', linewidth=2.2, label=f'Chamber Design Channel ({CHChannelDia*1000:.1f} mm)')
plt.title('Channel Diameter vs. Water Velocity', fontsize=12, fontweight='bold')
plt.xlabel('Channel Diameter (mm)', fontsize=10)
plt.ylabel('Water Velocity (m/s)', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()


# --- Figure 4: Hot Wall Temp vs. Required Water Velocity ---
HotWallTargetSweep = np.linspace(450, 850, 200) 

# --- Loop 3A: Throat Velocity Requirements (Bartz Loop Back-calc) ---
GasTemp, GasDynamicViscosity, GasSpecificHeat, GasPr, Gamma, LocalDia, Mach_local, LocalChannelDia = th_GasTemp, th_GasDynamicViscosity, th_GasSpecificHeat, th_GasPr, th_Gamma, th_Dia, 1.0, THChannelDia
th_RequiredVels, th_ValidHW = [], []
for TargetHW in HotWallTargetSweep:
    sigma_req = ( (0.5 * TargetHW / GasTemp * (1 + (Gamma - 1) / 2 * Mach_local**2) + 0.5)**(0.8 - omega / 5.0) * (1 + (Gamma - 1) / 2 * Mach_local**2)**(omega/5.0) )**-1
    GasCoeff_req = ( (0.026 / LocalDia**0.2) * (GasDynamicViscosity**0.2 * GasSpecificHeat / (GasPr**0.6)) * (Pc / Cstar)**0.8 * (LocalDia / RadOfCurve)**0.1 * sigma_req )
    flux_req = GasCoeff_req * (GasTemp - TargetHW)
    ColdWall_req = TargetHW - flux_req * (WallThickness / WallConductivity)
    if ColdWall_req > WaterTemp:
        WaterCoeff_req = flux_req / (ColdWall_req - WaterTemp)
        Nusselt_req = WaterCoeff_req * LocalChannelDia / WaterThermalConductivity
        if Nusselt_req > 0:
            Re_req = (Nusselt_req / (0.023 * WaterPr**0.4))**(1.0 / 0.8)
            th_RequiredVels.append(Re_req * WaterDynamicViscosity / (WaterDensity * LocalChannelDia))
            th_ValidHW.append(TargetHW)

# --- Loop 3B: Chamber Velocity Requirements (Dittus-Boelter Direct Back-calc) ---
GasTemp, GasDensity, GasVelocity, GasDynamicViscosity, GasSpecificHeat, GasThermalConductivity, GasPr, LocalDia, LocalChannelDia = ch_GasTemp, ch_GasDensity, ch_GasVelocity, ch_GasDynamicViscosity, ch_GasSpecificHeat, ch_GasThermalConductivity, ch_GasPr, ch_Dia, CHChannelDia
ch_RequiredVels, ch_ValidHW = [], []

ReGas = GasDensity * GasVelocity * LocalDia / GasDynamicViscosity
GasCoeff_req = (0.023 * ReGas**0.8 * GasPr**0.3) * GasThermalConductivity / LocalDia

for TargetHW in HotWallTargetSweep:
    flux_req = GasCoeff_req * (GasTemp - TargetHW)
    ColdWall_req = TargetHW - flux_req * (WallThickness / WallConductivity)
    if ColdWall_req > WaterTemp:
        WaterCoeff_req = flux_req / (ColdWall_req - WaterTemp)
        Nusselt_req = WaterCoeff_req * LocalChannelDia / WaterThermalConductivity
        if Nusselt_req > 0:
            Re_req = (Nusselt_req / (0.023 * WaterPr**0.4))**(1.0 / 0.8)
            ch_RequiredVels.append(Re_req * WaterDynamicViscosity / (WaterDensity * LocalChannelDia))
            ch_ValidHW.append(TargetHW)

plt.figure(4, figsize=(11, 6))
plt.plot(th_ValidHW, th_RequiredVels, color='forestgreen', linewidth=3, label='Required Velocity - Throat (Bartz Model)')
plt.plot(ch_ValidHW, ch_RequiredVels, color='darkolivegreen', linewidth=3, label='Required Velocity - Chamber (D-B Model)')
plt.axvline(x=YieldTemp, color='purple', linestyle='--', linewidth=2.2, label=f'Material Yield Temperature ({YieldTemp:.1f} K)')
plt.yscale('log') 
plt.title('Target Hot Wall Temperature vs. Required Water Velocity', fontsize=12, fontweight='bold')
plt.xlabel('Target Hot Wall Temperature (K)', fontsize=10)
plt.ylabel('Required Water Velocity (m/s) [Log Scale]', fontsize=10)
plt.grid(True, which="both", linestyle='--', alpha=0.6)
plt.legend(loc='best')
plt.tight_layout()

plt.show()