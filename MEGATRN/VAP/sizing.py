from math import pi, sqrt, isnan, atan2, cos, sin, acos, asin, ceil, floor
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

stiffness = 0.20
sleeve_thickness = 0.100
pintle_wall = 0.100
blockage_factor = 0.85

Cdr = 0.61
Cda = 0.98

Pa = 14.695

max_thrust = 900 # lbf
Pc = 200 # psi
of_ratio = 1
Mdot = 1.90997 / lb_to_kg
Mf = Mdot / (1 + of_ratio)
Mo = Mf * of_ratio
rhof = 785 * (in_to_m ** 3) / lb_to_kg # Ethanol
rhoo = 1205 / lb_to_kg * (in_to_m ** 3) # LOx

expansion_ratio = 2.704
throat_diam = 0.05271 / in_to_m
chamber_diameter = 4.2 # in
chamber_od = 5.0 # in
exit_area = expansion_ratio * throat_diam**2 * pi * 0.25

Dp = 0.2 * chamber_diameter

mill_diam = 0.040

friction = 0.8
servo_encoder_count = 6400

#======== END OF PARAMETERS ========
Dbore = floor(16.0 * (Dp - 2.0 * pintle_wall)) / 16.0

slot_count = round(blockage_factor * pi * Dbore / mill_diam / 2.0) * 2.0

Aslot = Mo / sqrt(2 * rhoo * stiffness * Pc * lbf_to_lbins) / Cdr
print(f"Slot Area: {Aslot:.6f}")
slot_height = Aslot / (slot_count * mill_diam) + (1 - pi / 4) * mill_diam

Rmomentum = Mo ** 2 / rhoo / Aslot

Aann = Mf ** 2 / (1.5 * Rmomentum) / rhof
Dann = sqrt((4 * Aann) / pi + (Dp + 2 * sleeve_thickness)** 2)

Pdf = (Mf ** 2 / (2 * Cda ** 2 * Aann ** 2 * rhof)) / lbf_to_lbins
ox_feed = (1 + stiffness) * Pc
fuel_feed = Pdf + Pc

min_rp1_diam = 2 * sqrt(Aann / pi)
min_lox_diam = 2 * sqrt(Aslot / pi)

delta_ann = (Dann - Dp - 2 * sleeve_thickness)/2

print(f"Pintle Tip Outer Diameter (in): {Dp:.3f}")
print(f"Pintle Tip Inner Diameter (in): {Dbore:.3f}")
print(f"Pintle Slot Count: {slot_count}")
print(f"Pintle Slot Height (in): {slot_height:.3f}")
print(f"Actual Blockage Factor: {slot_count * mill_diam / pi / Dbore:.3f}")
print(f"Annulus Outer Diameter (in): {Dann:.3f}")
print(f"Annulur Gap (in): {delta_ann:.3f}")
print(f"Fuel Stiffness (%): {100 * Pdf/Pc:.1f}")
print(f"LOx Feed Pressure (psia): {ox_feed: .0f}")
print(f"Fuel Feed Pressure (psia): {fuel_feed: .0f}")
print(f"Mininum LOx Inlet Diameter (in): {min_lox_diam:.3f}")
print(f"Mininum Fuel Inlet Diameter (in): {min_rp1_diam:.3f}")

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
lox_vels = np.empty_like(throttles)
fuel_vels = np.empty_like(throttles)
for i in range(len(throttles)):
    lox_mom = (ox_fraction * mass_flows[i]) ** 2 / rhoo / lox_areas[i]
    fuel_mom = (fuel_fraction * mass_flows[i]) ** 2 / rhof / fuel_areas[i]
    tmrs[i] = lox_mom / fuel_mom
    lmrs[i] = tmrs[i] / blockage_factor
    spray_angles[i] = (105.5 - 24.5 * tmrs[i]) * sqrt(lmrs[i])
    lox_vels[i] = (ox_fraction * mass_flows[i]) / rhoo / lox_areas[i]
    fuel_vels[i] = (fuel_fraction * mass_flows[i]) / rhof / fuel_areas[i]



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
alignment_limit = 0.002
min_slotheight = None
min_throttle = None
sleeve_50 = None
min_pc = None


cap1_limit = 1 - 0.125 * slot_count * pi * mill_diam**2 / Aslot
cap2_limit = 0.125 * slot_count * pi * mill_diam**2 / Aslot
cap1_height = slot_height - mill_diam * 0.5
cap2_height = mill_diam * 0.5
cap1_gap = None
cap2_gap = None
for i in range(len(throttles)):
    annulargaps[i] = (Dann - sqrt(Dann**2 - 4 * fuel_areas[i] / pi)) / 2
    if throttles[i] >= cap1_limit:
        theta = 1.0
        err = 1.0
        while abs(err/(throttles[i]*Aslot/slot_count)) > 0.00000001:
            err = Aslot / slot_count * (1.0 - throttles[i]) - mill_diam**2 * 0.125 * (theta - sin(theta))
            theta -= err / (-1.0 * mill_diam**2 * 0.125 * (1.0 - cos(theta)))
        slotheights[i] = slot_height + mill_diam * 0.5 * (cos(theta * 0.5) - 1.0)
    elif throttles[i] > cap2_limit:
        slotheights[i] = throttles[i] * Aslot / slot_count / mill_diam + (0.5 - pi * 0.125) * mill_diam
        if cap1_height is None:
            cap1_height = slotheights[i]
            cap1_gap = annulargaps[i]
    else:
        theta = 1.0
        err = 1.0
        while abs(err/(throttles[i]*Aslot/slot_count)) > 0.00000001:
            err = 0.125 * mill_diam**2 * (theta - sin(theta)) - throttles[i] * Aslot / slot_count
            theta -= err / (0.125 * mill_diam**2 * (1.0 - cos(theta)))
        slotheights[i] = mill_diam * 0.5 * (1.0 - cos(theta * 0.5))
    if annulargaps[i] <= alignment_limit and min_slotheight is None:
        min_slotheight = slotheights[i]
        min_throttle = throttles[i]
        min_pc = chamber_pressures[i]
    if throttles[i] == 0.50:
        sleeve_50 = (Dann - Dp) / 2 - annulargaps[50]

with open("sleeve.txt", "w") as file:
    for i in range(len(throttles)):
        file.write(f"{annulargaps[i]-delta_ann:f},0.0,{slot_height-slotheights[i]:f}\n")

sleeve_max = (Dann - Dp) / 2 - np.min(annulargaps)

optimal_arm = chamber_od / 2.0


print(f"Sleeve Thickness at 50% Throttle (in): {sleeve_50:.3f}")
print(f"Sleeve Thickness at 0% Throttle (in): {sleeve_max:.3f}")
print(f"Reliable Minimum Throttle (%): {min_throttle * 100:.0f}")
print(f"Minimum Expected Chamber Pressure (psia): {min_pc:.0f}")
print(f"Approx. Throttle Accuracy (%): {100.0 * (chamber_od * 0.5 * 2 * pi / servo_encoder_count) / slot_height:.1f}")
print(f"Downwards Force (lbf): {fuel_feed * pi * 0.25 * (Dann**2 - (Dp+sleeve_thickness*2)**2):.1f}")

get_oring.load_orings()

foring = get_oring.inner_diam(chamber_diameter - 0.5)
foring_gland = get_oring.get_face_gland_size(foring['num'])

pintle_ring = get_oring.radial_min_rod_od(Dp)
pintle_gland = get_oring.get_radial_gland_size(pintle_ring['num'])

sleeve_ring = get_oring.radial_min_rod_od(
    pintle_gland['id'] + pintle_gland['depth'] + pintle_gland['spacing']
)
sleeve_gland = get_oring.get_radial_gland_size(sleeve_ring['num'])

print(f"Face O-Ring: #{foring['num']:03d}")
print(f"Sleeve O-Ring: #{sleeve_ring['num']:03d}")
print(f"Pintle O-Ring: #{pintle_ring['num']:03d}")

cv_min_diam = min_rp1_diam / 2 # Vibes
cv_areas = cv_manifold.calc_cv(
    Mf,
    in_to_m * min_rp1_diam / 2,
    rhof,
    4,
    10000 * pi * cv_min_diam**2 / 4 * in_to_m**2
)
cv_max_diam = sqrt(4 * cv_areas[0] / pi)

pintle_orb = orb.orb_min_thread(Dbore)
fuel_orb = orb.orb_min_diam(min_rp1_diam)
pt_orb = orb.dash_num(2)

print(f"Pintle (LOx) ORB: -{pintle_orb['num']} ({pintle_orb['tube']:.4f})")
print(f"Fuel ORB: -{fuel_orb['num']}")
print(f"PT ORB: -{pt_orb['num']}")

with open("master.exp", "w") as file:
    file.write(f"""
// Version:  3
[Inch]annular_gap={delta_ann:.4f}
[Degree]annulus_angle=60
[Inch]chamber_id={chamber_diameter:.4f}
[Inch]chamber_od={chamber_od:.4f}
[Inch]cv_max_diam={cv_max_diam:.4f}
[Inch]cv_min_diam={cv_min_diam:.4f}
[Degrees]drill_angle=118
[Inch]faceplate_thickness=0.125
[Inch]foring_depth={foring_gland['depth']:.4f}
[Inch]foring_id={foring['ID']:.4f}
[Inch]foring_spacing={foring_gland['spacing']:.4f}
[Inch]foring_width={foring_gland['width']:.4f}
[Inch]nub_height=faceplate_thickness * 1.5
[Inch]nub_thickness=faceplate_thickness * 0.75
[Inch]pintle_flange_extension=5/16
[Inch]pintle_flange_thickness=0.125
[Inch]pintle_id={Dbore:.4f}
[Inch]pintle_od={Dp:.4f}
[Inch]pintle_oring_depth={pintle_gland['depth']:.4f}
[Inch]pintle_oring_id={pintle_gland['id']:.4f}
[Inch]pintle_oring_width={pintle_gland['width']:.4f}
[Inch]plenum_height={cv_min_diam/2:.4f}
[Inch]sleeve_50_thickness={sleeve_50:.4f}
[Inch]sleeve_extrusion_id={pintle_orb['spotface_diam']:.4f}
[Inch]sleeve_extrusion_height=2.5
[Inch]sleeve_oring_depth={sleeve_gland['depth']:.4f}
[Inch]sleeve_oring_spacing={sleeve_gland['spacing']:.4f}
[Inch]sleeve_oring_width={sleeve_gland['width']:.4f}
[Inch]sleeve_pin_thickness=0.250
[Inch]sleeve_ring_od={sleeve_gland['id'] + sleeve_gland['depth']:.4f}
[Inch]sleeve_thickness={sleeve_thickness:.4f}
[Inch]slot_height={slot_height:.4f}
[Inch]pintle_chamfer_od={pintle_orb['chamfer_od']:.4f}
[Inch]pintle_thread_diam={pintle_orb['thread_diam']:.4f}
(String) pintle_thread_type="{pintle_orb['thread_type']}"
[Inch]pintle_chamfer_length={pintle_orb['chamfer_length']:.4f}
[Degrees]pintle_secondary_chamfer={pintle_orb['secondary_chamfer_angle']:.4f}
[Degrees]pintle_chamfer_angle={pintle_orb['chamfer_angle']:.4f}
[Inch]pintle_tap_length={pintle_orb['thread_length']:.4f}
[Inch]fuel_chamfer_od={fuel_orb['chamfer_od']:.4f}
[Inch]fuel_thread_diam={fuel_orb['thread_diam']:.4f}
(String) fuel_thread_type="{fuel_orb['thread_type']}"
[Inch]fuel_chamfer_length={fuel_orb['chamfer_length']:.4f}
[Degrees]fuel_secondary_chamfer={fuel_orb['secondary_chamfer_angle']:.4f}
[Degrees]fuel_chamfer_angle={fuel_orb['chamfer_angle']:.4f}
[Inch]fuel_tap_length={fuel_orb['thread_length']:.4f}
[Inch]fuel_pipe_diam={fuel_orb['min_id']:.4f}
[Inch]fuel_spotface_diam={fuel_orb['spotface_diam']:.4f}
[Inch]pt_chamfer_od={pt_orb['chamfer_od']:.4f}
[Inch]pt_thread_diam={pt_orb['thread_diam']:.4f}
(String) pt_thread_type="{pt_orb['thread_type']}"
[Inch]pt_chamfer_length={pt_orb['chamfer_length']:.4f}
[Degrees]pt_secondary_chamfer={pt_orb['secondary_chamfer_angle']:.4f}
[Degrees]pt_chamfer_angle={pt_orb['chamfer_angle']:.4f}
[Inch]pt_tap_length={pt_orb['thread_length']:.4f}
[Inch]pt_pipe_diam={pt_orb['min_id']:.4f}
[Inch]slot_width={mill_diam:.4f}
slot_count={slot_count}
bolt_count=4
(String) bolt_type="1/4-20"
[Inch]servo_width=3.4
[Degrees]extrusion_angle=30
(String) pintle_bolt_type="#8-32"
[Inch]pintle_bolt_diam={5/32.0}
pintle_bolt_count=2
[Inch]pin_diam=0.125
""")


fig1, ax1 = plt.subplots(2, 3)
fig2, ax2 = plt.subplots(2, 2)
fig3, ax3 = plt.subplots(1, 1)
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
ax2[1,0].vlines(min_throttle, 0.0, max(np.max(slotheights), np.max(annulargaps)), linestyles='dashed')
ax2[1,0].vlines(cap2_limit, 0.0, max(np.max(slotheights), np.max(annulargaps)), linestyles='dashed', colors=['red'])
ax2[1,0].vlines(cap1_limit, 0.0, max(np.max(slotheights), np.max(annulargaps)), linestyles='dashed', colors=['green'])
ax2[1,0].set_xlabel("Thrust Level (actual/full)")
ax2[1,0].set_ylabel("Size (in)")
ax2[1,0].set_title("Slot Height and Annular Gap vs. Throttle Level")
ax2[1,0].legend(["Slot Height", "Annular Gap", "Reliable Minimum Throttle", "Bottom Cap", "Top Cap"])
ax2[1,0].grid()

ax2[1,1].plot(slotheights, annulargaps)
ax2[1,1].set_xlim(np.min(slotheights), np.max(slotheights))
ax2[1,1].set_ylim(0.0, np.max(annulargaps))
ax2[1,1].vlines(min_slotheight, 0.0, np.max(annulargaps), linestyles='dashed')
ax2[1,1].vlines(cap2_height, 0.0, np.max(annulargaps), linestyles='dashed', colors=['red'])
ax2[1,1].vlines(cap1_height, 0.0, np.max(annulargaps), linestyles='dashed', colors=['green'])
ax2[1,1].set_xlabel("Slot Height (in)")
ax2[1,1].set_ylabel("Annular Gap (in)")
ax2[1,1].set_title("Annular Gap vs. Slot Height")
ax2[1,1].legend(["Annular Gap", "Reliable Minimum Throttle", "Bottom Cap", "Top Cap"])
ax2[1,1].grid()

ax3.plot(throttles, lox_vels)
ax3.plot(throttles, fuel_vels)
ax3.legend(["LOx Velocity", "Fuel Velocity"])
ax3.grid()

plt.show()
