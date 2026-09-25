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

#Sub-component Wall Thicknesses
wall_thickness_chamber_mm = 2.5 # mm
wall_thickness_throat_mm = 1.25  # mm


#throat pr/t = 5502.44 yield required with 1.3x SF 75400 pSI is tensile yield for 7175, so minimum yield strength percent has to be 7.29 % now check graph for temp, T = 650F/616.5k at 30min of exposure (1.25mm wall)
#Chamber pr/t = 4638 yield required with 1.3x SF 75400 pSI is tensile yield for 7175, so minimum yield strength percent has to be 6.15 % now check graph for temp, T = 650F/616.5k at 30min of exposure (2.5mm wall)

WallConductivity = 142          # W/m-K (7175 Aluminum)
Pc = 1.379e+6                   # Pa (200 psi)
Cstar = 2113.8                  # m/s (from CEA)
YieldTemp = 650                 # K (Material structural yielding point parameter)
watersattemp = 462.235          # K (Water saturation temperature at 200 psi)

# Global Constants & Iteration Setup
th_WallThickness = wall_thickness_throat_mm / 1000.0   # Convert to meters
ch_WallThickness = wall_thickness_chamber_mm / 1000.0  # Convert to meters

WaterMassflow = PumpRate * WaterDensity * 0.0000010515 
WaterPr = WaterSpecificHeat * WaterDynamicViscosity / WaterThermalConductivity
Tolerance = 0.1 # Kelvin
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

print(f"--- Nominal Converged States ---")
print(f"Throat  [Bartz] (ChDia={THChannelDia*1000:.1f}mm | Wall={wall_thickness_throat_mm:.2f}mm) -> Coolant Vel: {th_WaterVelocity:.2f} m/s | Hot Wall: {th_HotWallTemp:.2f} K | Cold Wall: {th_ColdWallTemp:.2f} K | Flux: {th_flux/1e6:.3f} MW/m^2")
print(f"Chamber [D-B]   (ChDia={CHChannelDia*1000:.1f}mm | Wall={wall_thickness_chamber_mm:.2f}mm) -> Coolant Vel: {ch_WaterVelocity:.2f} m/s | Hot Wall: {ch_HotWallTemp:.2f} K | Cold Wall: {ch_ColdWallTemp:.2f} K | Flux: {ch_flux/1e6:.3f} MW/m^2\n")
print(f"Bartz Heat Transfer Coeff: {GasTransferCoeff:.2f} W/m^2-K")

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


# --- Figure 3: Channel Diameter Sweep vs. Water Velocity ---
DiaSweep = np.linspace(0.001, 0.006, 100) 
VelocitySweep = WaterMassflow / (ChannelCount * (np.pi * (DiaSweep/2)**2) * WaterDensity)


# --- Figure 4: Hot Wall Temp vs. Required Water Velocity ---
HotWallTargetSweep = np.linspace(450, 850, 200) 

# --- Loop 3A: Throat Velocity Requirements (Bartz Loop Back-calc) ---
GasTemp, GasDynamicViscosity, GasSpecificHeat, GasPr, Gamma, LocalDia, Mach_local, LocalChannelDia = th_GasTemp, th_GasDynamicViscosity, th_GasSpecificHeat, th_GasPr, th_Gamma, th_Dia, 1.0, THChannelDia
th_RequiredVels, th_ValidHW = [], []
for TargetHW in HotWallTargetSweep:
    sigma_req = ( (0.5 * TargetHW / GasTemp * (1 + (Gamma - 1) / 2 * Mach_local**2) + 0.5)**(0.8 - omega / 5.0) * (1 + (Gamma - 1) / 2 * Mach_local**2)**(omega/5.0) )**-1
    GasCoeff_req = ( (0.026 / LocalDia**0.2) * (GasDynamicViscosity**0.2 * GasSpecificHeat / (GasPr**0.6)) * (Pc / Cstar)**0.8 * (LocalDia / RadOfCurve)**0.1 * sigma_req )
    flux_req = GasCoeff_req * (GasTemp - TargetHW)
    ColdWall_req = TargetHW - flux_req * (th_WallThickness / WallConductivity)
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
    ColdWall_req = TargetHW - flux_req * (ch_WallThickness / WallConductivity)
    if ColdWall_req > WaterTemp:
        WaterCoeff_req = flux_req / (ColdWall_req - WaterTemp)
        Nusselt_req = WaterCoeff_req * LocalChannelDia / WaterThermalConductivity
        if Nusselt_req > 0:
            Re_req = (Nusselt_req / (0.023 * WaterPr**0.4))**(1.0 / 0.8)
            ch_RequiredVels.append(Re_req * WaterDynamicViscosity / (WaterDensity * LocalChannelDia))
            ch_ValidHW.append(TargetHW)



#PLOTS (METRIC)


TITLE_FONT = 20
LABEL_FONT = 16
TICK_FONT = 14

# --- Combined Metric Window 1: Structural & Fluid Wall Temperatures ---
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Subplot 1 (Old Fig 1)
ax1.plot(ThicknessSweep * 1000, th_SweepHW, label='Throat (Bartz) - Hot Wall', color='crimson', linewidth=3)
ax1.plot(ThicknessSweep * 1000, ch_SweepHW, label='Chamber (D-B) - Hot Wall', color='darkred', linewidth=3)
ax1.axhline(y=YieldTemp, color='purple', linestyle=':', linewidth=2.5, label=f'Material Yield Temperature ({YieldTemp:.1f} K)')
ax1.plot(wall_thickness_throat_mm, th_HotWallTemp, marker='o', markersize=10, color='black', label=f'Current Throat Wall ({wall_thickness_throat_mm:.2f} mm)')
ax1.plot(wall_thickness_chamber_mm, ch_HotWallTemp, marker='s', markersize=10, color='blue', label=f'Current Chamber Wall ({wall_thickness_chamber_mm:.2f} mm)')
ax1.set_title('Wall Thickness vs. Structural Hot Wall Temperatures (metric)', fontsize=TITLE_FONT, fontweight='bold')
ax1.set_xlabel('Wall Thickness (mm)', fontsize=LABEL_FONT)
ax1.set_ylabel('Hot Wall Temperature (K)', fontsize=LABEL_FONT)
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend(loc='upper left', fontsize=TICK_FONT)
ax1.tick_params(axis='both', labelsize=TICK_FONT)

# Subplot 2 (Old Fig 2)
ax2.plot(ThicknessSweep * 1000, th_SweepCW, label='Throat (Bartz) - Cold Wall', color='royalblue', linewidth=3)
ax2.plot(ThicknessSweep * 1000, ch_SweepCW, label='Chamber (D-B) - Cold Wall', color='navy', linewidth=3)
ax2.axhline(y=watersattemp, color='darkorange', linestyle='--', linewidth=2.5, label='Water Saturation Temp (200 psi)')
ax2.plot(wall_thickness_throat_mm, th_ColdWallTemp, marker='o', markersize=10, color='black', label=f'Current Throat Wall ({wall_thickness_throat_mm:.2f} mm)')
ax2.plot(wall_thickness_chamber_mm, ch_ColdWallTemp, marker='s', markersize=10, color='blue', label=f'Current Chamber Wall ({wall_thickness_chamber_mm:.2f} mm)')
ax2.set_title('Wall Thickness vs. Fluid Interface Cold Wall Temperatures (metric)', fontsize=TITLE_FONT, fontweight='bold')
ax2.set_xlabel('Wall Thickness (mm)', fontsize=LABEL_FONT)
ax2.set_ylabel('Cold Wall Temperature (K)', fontsize=LABEL_FONT)
ax2.set_ylim(375, 490)
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.legend(loc='upper right', fontsize=TICK_FONT)
ax2.tick_params(axis='both', labelsize=TICK_FONT)

plt.tight_layout()

# --- Combined Metric Window 2: Fluid Dynamics & Velocity Curves ---
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(14, 6))

# Subplot 3 (Old Fig 3)
ax3.plot(DiaSweep * 1000, VelocitySweep, color='teal', linewidth=3, label='Velocity Curve')
ax3.axvline(x=THChannelDia*1000, color='crimson', linestyle=':', linewidth=2.5, label=f'Throat Design Channel ({THChannelDia*1000:.1f} mm)')
ax3.axvline(x=CHChannelDia*1000, color='darkred', linestyle='--', linewidth=2.5, label=f'Chamber Design Channel ({CHChannelDia*1000:.1f} mm)')
ax3.set_title('Channel Diameter vs. Water Velocity (metric)', fontsize=TITLE_FONT, fontweight='bold')
ax3.set_xlabel('Channel Diameter (mm)', fontsize=LABEL_FONT)
ax3.set_ylabel('Water Velocity (m/s)', fontsize=LABEL_FONT)
ax3.grid(True, linestyle='--', alpha=0.6)
ax3.legend(fontsize=TICK_FONT)
ax3.tick_params(axis='both', labelsize=TICK_FONT)

# Subplot 4 (Old Fig 4)
ax4.plot(th_ValidHW, th_RequiredVels, color='forestgreen', linewidth=3, label='Required Velocity - Throat (Bartz)')
ax4.plot(ch_ValidHW, ch_RequiredVels, color='darkolivegreen', linewidth=3, label='Required Velocity - Chamber (D-B)')
ax4.axvline(x=YieldTemp, color='purple', linestyle='--', linewidth=2.5, label=f'Material Yield Temperature ({YieldTemp:.1f} K)')
ax4.set_yscale('log') 
ax4.set_title('Target Hot Wall Temperature vs. Required Water Velocity (metric)', fontsize=TITLE_FONT, fontweight='bold')
ax4.set_xlabel('Target Hot Wall Temperature (K)', fontsize=LABEL_FONT)
ax4.set_ylabel('Required Water Velocity (m/s) [Log Scale]', fontsize=LABEL_FONT)
ax4.grid(True, which="both", linestyle='--', alpha=0.6)
ax4.legend(loc='best', fontsize=TICK_FONT)
ax4.tick_params(axis='both', labelsize=TICK_FONT)

plt.tight_layout()



#PLOTS (IMPERIAL)

# Dynamic Unit Conversion Handlers
def k_to_f(k): return (k * 1.8) - 459.67
def mm_to_in(mm): return mm / 25.4
def m_to_ft(m): return m * 3.2808399

# Vectorized conversions of upstream computed arrays
ThicknessSweep_in = mm_to_in(ThicknessSweep * 1000)
th_SweepHW_F = k_to_f(np.array(th_SweepHW))
ch_SweepHW_F = k_to_f(np.array(ch_SweepHW))
th_SweepCW_F = k_to_f(np.array(th_SweepCW))
ch_SweepCW_F = k_to_f(np.array(ch_SweepCW))
DiaSweep_in = mm_to_in(DiaSweep * 1000)
VelocitySweep_fts = m_to_ft(VelocitySweep)
th_ValidHW_F = k_to_f(np.array(th_ValidHW))
ch_ValidHW_F = k_to_f(np.array(ch_ValidHW))
th_RequiredVels_fts = m_to_ft(np.array(th_RequiredVels))
ch_RequiredVels_fts = m_to_ft(np.array(ch_RequiredVels))

# Conversion of tracking markers and key static points
YieldTemp_F = k_to_f(YieldTemp)
watersattemp_F = k_to_f(watersattemp)
wt_throat_in = mm_to_in(wall_thickness_throat_mm)
wt_chamber_in = mm_to_in(wall_thickness_chamber_mm)
th_HW_F = k_to_f(th_HotWallTemp)
ch_HW_F = k_to_f(ch_HotWallTemp)
th_CW_F = k_to_f(th_ColdWallTemp)
ch_CW_F = k_to_f(ch_ColdWallTemp)
th_ch_dia_in = mm_to_in(THChannelDia * 1000)
ch_ch_dia_in = mm_to_in(CHChannelDia * 1000)

# --- Combined Imperial Window 1: Structural & Fluid Wall Temperatures ---
fig3, (ax1_imp, ax2_imp) = plt.subplots(1, 2, figsize=(14, 6))

# Subplot 1 (Imperial Fig 1)
ax1_imp.plot(ThicknessSweep_in, th_SweepHW_F, label='Throat (Bartz) - Hot Wall', color='crimson', linewidth=3)
ax1_imp.plot(ThicknessSweep_in, ch_SweepHW_F, label='Chamber (D-B) - Hot Wall', color='darkred', linewidth=3)
ax1_imp.axhline(y=YieldTemp_F, color='purple', linestyle=':', linewidth=2.5, label=f'Material Yield Temperature ({YieldTemp_F:.1f} °F)')
ax1_imp.plot(wt_throat_in, th_HW_F, marker='o', markersize=10, color='black', label=f'Current Throat Wall ({wt_throat_in:.3f} in)')
ax1_imp.plot(wt_chamber_in, ch_HW_F, marker='s', markersize=10, color='blue', label=f'Current Chamber Wall ({wt_chamber_in:.3f} in)')
ax1_imp.set_title('Wall Thickness vs. Structural Hot Wall Temperatures', fontsize=TITLE_FONT, fontweight='bold')
ax1_imp.set_xlabel('Wall Thickness (in)', fontsize=LABEL_FONT)
ax1_imp.set_ylabel('Hot Wall Temperature (°F)', fontsize=LABEL_FONT)
ax1_imp.grid(True, linestyle='--', alpha=0.6)
ax1_imp.legend(loc='upper left', fontsize=TICK_FONT)
ax1_imp.tick_params(axis='both', labelsize=TICK_FONT)

# Subplot 2 (Imperial Fig 2)
ax2_imp.plot(ThicknessSweep_in, th_SweepCW_F, label='Throat (Bartz) - Cold Wall', color='royalblue', linewidth=3)
ax2_imp.plot(ThicknessSweep_in, ch_SweepCW_F, label='Chamber (D-B) - Cold Wall', color='navy', linewidth=3)
ax2_imp.axhline(y=watersattemp_F, color='darkorange', linestyle='--', linewidth=2.5, label='Water Saturation Temp (200 psi)')
ax2_imp.plot(wt_throat_in, th_CW_F, marker='o', markersize=10, color='black', label=f'Current Throat Wall ({wt_throat_in:.3f} in)')
ax2_imp.plot(wt_chamber_in, ch_CW_F, marker='s', markersize=10, color='blue', label=f'Current Chamber Wall ({wt_chamber_in:.3f} in)')
ax2_imp.set_title('Wall Thickness vs. Fluid Interface Cold Wall Temperatures', fontsize=TITLE_FONT, fontweight='bold')
ax2_imp.set_xlabel('Wall Thickness (in)', fontsize=LABEL_FONT)
ax2_imp.set_ylabel('Cold Wall Temperature (°F)', fontsize=LABEL_FONT)
ax2_imp.set_ylim(k_to_f(375), k_to_f(490))
ax2_imp.grid(True, linestyle='--', alpha=0.6)
ax2_imp.legend(loc='upper right', fontsize=TICK_FONT)
ax2_imp.tick_params(axis='both', labelsize=TICK_FONT)

plt.tight_layout()

# --- Combined Imperial Window 2: Fluid Dynamics & Velocity Curves ---
fig4, (ax3_imp, ax4_imp) = plt.subplots(1, 2, figsize=(14, 6))

# Subplot 3 (Imperial Fig 3)
ax3_imp.plot(DiaSweep_in, VelocitySweep_fts, color='teal', linewidth=3, label='Velocity Curve')
ax3_imp.axvline(x=th_ch_dia_in, color='crimson', linestyle=':', linewidth=2.5, label=f'Throat Design Channel ({th_ch_dia_in:.3f} in)')
ax3_imp.axvline(x=ch_ch_dia_in, color='darkred', linestyle='--', linewidth=2.5, label=f'Chamber Design Channel ({ch_ch_dia_in:.3f} in)')
ax3_imp.set_title('Channel Diameter vs. Water Velocity', fontsize=TITLE_FONT, fontweight='bold')
ax3_imp.set_xlabel('Channel Diameter (in)', fontsize=LABEL_FONT)
ax3_imp.set_ylabel('Water Velocity (ft/s)', fontsize=LABEL_FONT)
ax3_imp.grid(True, linestyle='--', alpha=0.6)
ax3_imp.legend(fontsize=TICK_FONT)
ax3_imp.tick_params(axis='both', labelsize=TICK_FONT)

# Subplot 4 (Imperial Fig 4)
ax4_imp.plot(th_ValidHW_F, th_RequiredVels_fts, color='forestgreen', linewidth=3, label='Required Velocity - Throat (Bartz)')
ax4_imp.plot(ch_ValidHW_F, ch_RequiredVels_fts, color='darkolivegreen', linewidth=3, label='Required Velocity - Chamber (D-B)')
ax4_imp.axvline(x=YieldTemp_F, color='purple', linestyle='--', linewidth=2.5, label=f'Material Yield Temperature ({YieldTemp_F:.1f} °F)')
ax4_imp.set_yscale('log') 
ax4_imp.set_title('Target Hot Wall Temperature vs. Required Water Velocity', fontsize=TITLE_FONT, fontweight='bold')
ax4_imp.set_xlabel('Target Hot Wall Temperature (°F)', fontsize=LABEL_FONT)
ax4_imp.set_ylabel('Required Water Velocity (ft/s) [Log Scale]', fontsize=LABEL_FONT)
ax4_imp.grid(True, which="both", linestyle='--', alpha=0.6)
ax4_imp.legend(loc='best', fontsize=TICK_FONT)
ax4_imp.tick_params(axis='both', labelsize=TICK_FONT)

plt.tight_layout()
plt.show()