from math import pi, sqrt, isnan, atan2, cos, sin, acos, asin, ceil
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

import get_oring
import cv_manifold
import orb

in_to_m = 0.0254
lb_to_kg = 0.4536
N_to_lbf = 0.2248
lbf_to_lbins = 386.09

def get_throttle_target(Pc, OF, ep):
    rhoe = cea.get_Densities(Pc, OF, ep)[2] / (12 ** 3)
    mach = cea.get_MachNumber(Pc, OF, ep)
    sonic = cea.get_SonicVelocities(Pc, OF, ep)[2] * 12
    Pe = Pc / cea.get_PcOvPe(Pc, OF, ep)

    actual = Pe + (rhoe * (mach * sonic) ** 2) / lbf_to_lbins

    return actual


stiffness = 0.15
sleeve_thickness = 0.100
pintle_wall = 0.100
blockage_factor = 0.85

Cdr = 0.61
Cda = 0.98

Pa = 14.695

max_thrust = 1060
Pc = 225
of_ratio = 1.5
Mdot = 4.36
Mf = Mdot / (1 + of_ratio)
Mo = Mf * of_ratio
rhof = 785 * (in_to_m ** 3) / lb_to_kg
rhoo = 1205 / lb_to_kg * (in_to_m ** 3)

expansion_ratio = 3.274
chamber_diameter = (0.107) / in_to_m
exit_area = expansion_ratio * ((0.05309 / 2) / in_to_m)**2 * pi

Dp = 0.2 * chamber_diameter
Dbore = Dp - 2 * pintle_wall

friction = 0.8
servo_encoder_count = 6400

#======== END OF PARAMETERS ========

Aslot = Mo / sqrt(2 * rhoo * stiffness * Pc * lbf_to_lbins) / Cdr
print(f"Slot Area: {Aslot:.6f}")
slot_height = Aslot / (pi * Dbore * blockage_factor)
slot_circum = Aslot / slot_height

Rmomentum = Mo ** 2 / rhoo / Aslot

Aann = Mf ** 2 / Rmomentum / rhof
Dann = sqrt((4 * Aann) / pi + (Dp + 2 * sleeve_thickness)** 2)

Pdf = (Mf ** 2 / (2 * Cda ** 2 * Aann ** 2 * rhof)) / lbf_to_lbins
ox_feed = (1 + stiffness) * Pc
fuel_feed = Pdf + Pc

min_rp1_diam = 2 * sqrt(Aann / pi)
min_lox_diam = 2 * sqrt(Aslot / pi)

print(f"Pintle Tip Outer Diameter (in): {Dp:.3f}")
print(f"Pintle Tip Inner Diameter (in): {Dbore:.3f}")
print(f"Pintle Slot Height (in): {slot_height:.3f}")
print(f"Pintle Slot Used Circumference (in): {slot_circum:.3f}")
print(f"Annulus Outer Diameter (in): {Dann:.3f}")
print(f"Annulur Gap (in): {(Dann - Dp - 2 * sleeve_thickness)/2:.3f}")
print(f"Fuel Stiffness (%): {100 * Pdf/Pc:.1f}")
print(f"LOx Feed Pressure (psia): {ox_feed: .0f}")
print(f"RP1 Feed Pressure (psia): {fuel_feed: .0f}")
print(f"Mininum LOx Inlet Diameter (in): {min_lox_diam:.3f}")
print(f"Mininum RP1 Inlet Diameter (in): {min_rp1_diam:.3f}")

cea = CEA_Obj(oxName='LOX', fuelName='Ethanol')

throttles = np.linspace(1.0, 0.01, 100)
chamber_pressures = np.empty_like(throttles)
for i in range(len(throttles)):
    throttle = throttles[i]

    desired = throttle * max_thrust / exit_area + Pa
    
    converged = False
    chamber_pressure = Pc * throttle
    while not converged:
        actual = get_throttle_target(chamber_pressure, Mo / Mf, expansion_ratio)
        actual -= desired
        if abs(actual / desired) < 0.0000001:
            converged = True
        else:
            epsilon = get_throttle_target(chamber_pressure + 0.01, Mo / Mf, expansion_ratio)
            epsilon -= desired

            slope = (epsilon - actual) / 0.01
            chamber_pressure -= actual / slope

    chamber_pressures[i] = chamber_pressure

mass_flows = np.empty_like(throttles)
for i in range(len(throttles)):
    rhoe = cea.get_Densities(chamber_pressures[i], Mo / Mf, expansion_ratio)[2] / (12 ** 3)
    mach = cea.get_MachNumber(chamber_pressures[i], Mo / Mf, expansion_ratio)
    sonic = cea.get_SonicVelocities(chamber_pressures[i], Mo / Mf, expansion_ratio)[2] * 12

    mass_flows[i] = rhoe * mach * sonic * exit_area

ox_fraction = Mo / (Mo + Mf)
fuel_fraction = Mf / (Mo + Mf)
lox_areas = np.empty_like(throttles)
fuel_areas = np.empty_like(throttles)
lox_stiffnesses = np.empty_like(throttles)
fuel_stiffnesses = np.empty_like(throttles)
for i in range(len(throttles)):
    lox_areas[i] = (ox_fraction * mass_flows[i]) / sqrt(2 * rhoo * (ox_feed - chamber_pressures[i]) * lbf_to_lbins) / Cdr
    fuel_areas[i] = (fuel_fraction * mass_flows[i]) / sqrt(2 * rhof * (fuel_feed - chamber_pressures[i]) * lbf_to_lbins) / Cda
    lox_stiffnesses[i] = (ox_feed - chamber_pressures[i]) / chamber_pressures[i]
    fuel_stiffnesses[i] = (fuel_feed - chamber_pressures[i]) / chamber_pressures[i]

tmrs = np.empty_like(throttles)
lmrs = np.empty_like(throttles)
spray_angles = np.empty_like(throttles)
for i in range(len(throttles)):
    lox_mom = (ox_fraction * mass_flows[i]) ** 2 / rhoo / lox_areas[i]
    fuel_mom = (fuel_fraction * mass_flows[i]) ** 2 / rhof / fuel_areas[i]
    tmrs[i] = lox_mom / fuel_mom
    lmrs[i] = tmrs[i] / blockage_factor
    spray_angles[i] = (105.5 - 24.5 * tmrs[i]) * sqrt(lmrs[i])


ispvacs = np.empty_like(throttles)
ispambs = np.empty_like(throttles)
separation_point = 0.0
for i in range(len(throttles)):
    ispvacs[i] = cea.get_Isp(chamber_pressures[i], Mo / Mf, expansion_ratio)
    ispambs[i], mode = cea.estimate_Ambient_Isp(chamber_pressures[i], Mo / Mf, expansion_ratio, Pa)
    if mode == "Separated" and throttles[i] > separation_point:
        separation_point = throttles[i]
        print(f"Separation at: {separation_point}")

slotheights = np.empty_like(throttles)
annulargaps = np.empty_like(throttles)
alignment_limit = 0.004
min_slotheight = None
min_throttle = None
for i in range(len(throttles)):
    slotheights[i] = lox_areas[i] / (pi * Dbore * blockage_factor)
    annulargaps[i] = (Dann - sqrt(Dann**2 - 4 * fuel_areas[i] / pi)) / 2
    if annulargaps[i] <= alignment_limit and min_slotheight is None:
        min_slotheight = slotheights[i]
        min_throttle = throttles[i]

sleeve_max = (Dann - Dp) / 2 - np.min(annulargaps)
print(f"Sleeve Thickness at 0% Throttle (in): {sleeve_max:.3f}")
print(f"Reliable Minimum Throttle (%): {min_throttle * 100:.0f}")

get_oring.load_orings()

radial_rings = get_oring.radial_min_rod_od(sleeve_max*2 + Dp)
radial_gland = get_oring.get_radial_gland_size(radial_rings['num'])
print(f"Middle Plate Radial Rings: 2x #{radial_rings['num']:03d}")

fig1, ax1 = plt.subplots(2, 3)
fig2, ax2 = plt.subplots(2, 2)
ax1[0,0].plot(throttles, chamber_pressures)
ax1[0,0].set_xlim(0.0, 1.0)
ax1[0,0].set_ylim(0.0, Pc)
ax1[0,0].set_xlabel("Throttle Level (actual/full)")
ax1[0,0].set_ylabel("Chamber Pressure (psia)")
ax1[0,0].set_title("Chamber Pressure vs. Throttle")
ax1[0,0].grid()

ax1[1,0].plot(throttles, mass_flows)
ax1[1,0].set_xlim(0.0, 1.0)
ax1[1,0].set_ylim(0.0, Mo + Mf)
ax1[1,0].set_xlabel("Throttle Level (actual/full)")
ax1[1,0].set_ylabel("Total Mass Flow (lbm / s)")
ax1[1,0].set_title("Total Mass Flow vs. Throttle")
ax1[1,0].grid()

ax1[0,1].plot(chamber_pressures, mass_flows)
ax1[0,1].set_xlim(np.min(chamber_pressures), Pc)
ax1[0,1].set_ylim(0.0, Mo + Mf)
ax1[0,1].set_xlabel("Chamber Pressure (psia)")
ax1[0,1].set_ylabel("Total Mass Flow (lbm / s)")
ax1[0,1].set_title("Total Mass Flow vs. Chamber Pressure")
ax1[0,1].grid()

ax2[0,0].plot(throttles, lox_areas)
ax2[0,0].plot(throttles, fuel_areas)
ax2[0,0].set_xlim(0.0, 1.0)
ax2[0,0].set_ylim(0.0, max(Aann, Aslot))
ax2[0,0].set_xlabel("Throttle Level (actual/full)")
ax2[0,0].set_ylabel("Orifice Areas (sq.in)")
ax2[0,0].set_title("Orifice Areas vs. Throttle Level")
ax2[0,0].legend(["LOx", "Ethanol"])
ax2[0,0].grid()

ax1[1,1].plot(throttles, tmrs)
ax1[1,1].plot(throttles, lmrs)
ax1[1,1].set_xlim(0.0, 1.0)
ax1[1,1].set_ylim(0.0, 1.05 * max(np.max(tmrs), np.max(lmrs)))
ax1[1,1].set_xlabel("Throttle Level (actual/full)")
ax1[1,1].set_ylabel("TMR and LMR (radial/annular momentum)")
ax1[1,1].set_title("TMR and LMR vs. Throttle Level")
ax1[1,1].legend(["TMR", "LMR"])
ax1[1,1].grid()

ax2[0,1].plot(throttles, 100.0 * lox_stiffnesses)
ax2[0,1].plot(throttles, 100.0 * fuel_stiffnesses)
ax2[0,1].set_xlim(0.0, 1.0)
ax2[0,1].set_ylim(0.0, 100.0 * max(np.max(lox_stiffnesses), np.max(fuel_stiffnesses)))
ax2[0,1].set_xlabel("Throttle Level (actual/full)")
ax2[0,1].set_ylabel("Stiffness (%)")
ax2[0,1].set_title("Stiffness vs. Throttle Level")
ax2[0,1].legend(["LOx", "Ethanol"])
ax2[0,1].grid()


ax1[0,2].plot(throttles, spray_angles)
ax1[0,2].set_xlim(0.0, 1.0)
ax1[0,2].set_ylim(np.min(spray_angles), np.max(spray_angles))
ax1[0,2].set_xlabel("Throttle Level (actual/full)")
ax1[0,2].set_ylabel("Spray Angle (degrees)")
ax1[0,2].set_title("Spray Angle vs. Throttle Level")
ax1[0,2].grid()

ax1[1,2].plot(throttles, ispvacs)
ax1[1,2].plot(throttles, ispambs)
ax1[1,2].vlines(separation_point, 0.0, 1.05 * np.max(ispvacs), linestyles='dashed', label='Flow Separation')
ax1[1,2].set_xlim(0.0, 1.0)
ax1[1,2].set_ylim(0.0, 1.05 * np.max(ispvacs))
ax1[1,2].set_xlabel("Thrust Level (actual/full)")
ax1[1,2].set_ylabel("ISP (s)")
ax1[1,2].set_title("Vacuum and Sea Level ISP vs. Throttle Level")
ax1[1,2].legend(["Vacuum", "Sea Level"])
ax1[1,2].grid()

ax2[1,0].plot(throttles, slotheights)
ax2[1,0].plot(throttles, annulargaps)
ax2[1,0].set_xlim(0.0, 1.0)
ax2[1,0].set_ylim(0.0, max(np.max(slotheights), np.max(annulargaps)))
ax2[1,0].vlines(min_throttle, 0.0, max(np.max(slotheights), np.max(annulargaps)), linestyles='dashed', label='Alignment Limit')
ax2[1,0].set_xlabel("Thrust Level (actual/full)")
ax2[1,0].set_ylabel("Size (in)")
ax2[1,0].set_title("Slot Height and Annular Gap vs. Throttle Level")
ax2[1,0].legend(["Slot Height", "Annular Gap", "Reliable Minimum Throttle"])
ax2[1,0].grid()

ax2[1,1].plot(slotheights, annulargaps)
ax2[1,1].set_xlim(np.min(slotheights), np.max(slotheights))
ax2[1,1].set_ylim(0.0, np.max(annulargaps))
ax2[1,1].vlines(min_slotheight, 0.0, np.max(annulargaps), linestyles='dashed', label='Alignment Limit')
ax2[1,1].set_xlabel("Slot Height (in)")
ax2[1,1].set_ylabel("Annular Gap (in)")
ax2[1,1].set_title("Annular Gap vs. Slot Height")
ax2[1,1].legend(["Annular Gap", "Reliable Minimum Throttle"])
ax2[1,1].grid()

plt.show()
