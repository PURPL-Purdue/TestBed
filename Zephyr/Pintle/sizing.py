from math import pi, sqrt, isnan, atan2, cos, sin
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

import get_oring

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


stiffness = 0.10
sleeve_thickness = 0.100
pintle_wall = 0.100
blockage_factor = 0.7

Cdr = 0.61
Cda = 0.98

Pa = 14.695

max_thrust = 1620
Pc = 325
of_ratio = 2
Mdot = 6.28
Mf = Mdot / (1 + of_ratio)
Mo = Mf * of_ratio
rhof = 768.09 * (in_to_m ** 3) / lb_to_kg
rhoo = 1205 / lb_to_kg * (in_to_m ** 3)

expansion_ratio = 4.101
chamber_diameter = (0.107) / in_to_m
exit_area = expansion_ratio * ((0.05337 / 2) / in_to_m)**2 * pi

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

print(f"Pintle Tip Outer Diameter (in): {Dp:.3f}")
print(f"Pintle Tip Inner Diameter (in): {Dbore:.3f}")
print(f"Pintle Slot Height (in): {slot_height:.3f}")
print(f"Pintle Slot Used Circumference (in): {slot_circum:.3f}")
print(f"Annulus Outer Diameter (in): {Dann:.3f}")
print(f"Annulur Gap (in): {(Dann - Dp - 2 * sleeve_thickness)/2:.3f}")
print(f"Fuel Stiffness (%): {100 * Pdf/Pc:.1f}")
print(f"LOx Feed Pressure (psia): {ox_feed: .0f}")
print(f"RP1 Feed Pressure (psia): {fuel_feed: .0f}")
print(f"Mininum LOx Inlet Diameter (in): {2*sqrt(Aslot / pi):.3f}")
print(f"Mininum RP1 Inlet Diameter (in): {2*sqrt(Aann / pi):.3f}")


thread_diam = None
thread_tpi = None
thread_name = None
with open('unf_threads.csv') as thread_csv:
    lines = thread_csv.readlines()
    for line in lines:
        data = line.split(",")
        diam = float(data[0])
        tpi = int(data[1])
        name = data[2].strip()
        if diam >= Dp + 0.100:
            thread_diam = diam
            thread_tpi = tpi
            thread_name = name
            break

if thread_diam is None or thread_tpi is None or thread_name is None:
    print(f"[ERROR]: Pintle Diamter {Dp:0.3f}in is larger than the largest thread size!")
    exit()

print(f"\n")
print(f"Pintle Thread: {thread_name} (OD: {thread_diam:.3f}; TPI: {thread_tpi})")


cea = CEA_Obj(oxName='LOX', fuelName='RP1')

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

print(f"{mass_flows}")

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
alignment_limit = 0.0005 / 0.15
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

sleeve_force = fuel_feed * pi * ((Dp + sleeve_max * 2)**2 - (Dp + sleeve_thickness * 2)**2) / 4.0
print(f"Force on Sleeve (lbf): {sleeve_force:.1f}")
thread_angle = atan2(1/thread_tpi, pi * thread_diam)
servo_torque = (thread_diam / 2) * sleeve_force * (friction * cos(thread_angle) + sin(thread_angle)) / (cos(thread_angle) - friction * sin(thread_angle))

print(f"Torque on Servo (lbf*in): {servo_torque:.1f}")
print(f"Actuation Precision (in): ±{thread_tpi/servo_encoder_count:.3f}")

get_oring.load_orings()

key_height = 0.100
thread_length = 3*gland['spacing']+2*gland['width']+2*slot_height+key_height

with open("pintle_shaft.exp", "w") as expr:
    expr.write(f"""
    // Version:  3
    [Inch]OR={Dp/2:.4f}
    [Inch]drill_tip=0.25
    [Inch]flange_rad=1.5
    [Inch]flange_thickness=0.25
    [Degrees]p0=0
    [Degrees]p1=360
    [Inch]p6=wall
    [Inch]p13=wall*sqrt(2)
    [Inch]p15=slot_circum/slot_count
    [Inch]p21=slot_height
    [Inch]p22=wall
    [Inch]p23=2.5*wall
    p31=slot_count
    [Degrees]p32=15
    [Inch]p33=1
    [Degrees]p34=360
    p35=1
    [Inch]p36=1
    [Inch]p37=0
    [Degrees]p38=0
    [Inch]p39=2*OR+3*slot_height
    [Inch]slot_circum={slot_circum:.3f}
    slot_count=24
    [Inch]slot_height={slot_height:.3f}
    [Inch]thread_rad={thread_diam/2:.4f}
    [Inch]upper_flat=2*slot_height+1.5
    [Inch]wall={pintle_wall:.3f}
    [Inch]thread_length={thread_length:.3f}
    """)

with open("sleeve.exp", "w") as expr:
    gland = get_oring.get_gland_size()

    expr.write(f"""
    // Version:  3
    [Inch]p0={sleeve_thickness:.3f}
    [Inch]p5={Dp/2:.2f}
    [Inch]p47={slot_height+Dp:.3f}
    [Inch]p48={slot_height:.3f}
    [Inch]p45={sleeve_max-sleeve_thickness:.3f}
    [Inch]p117={(thread_diam-Dp)/2-5*sqrt(3)/8/thread_tpi/2:.3f}
    [Inch]p134={thread_length:.3f}
    [Degrees]p7=360
    [Degrees]p102=45
    [Inch]p103={Dp/2+sleeve_max:.3f}
    [Inch]p50=sleeve_key_height
    [Inch]p116=middle_plate_height-lox_manifold_height-sleeve_key_height-gland['width']
    [Inch]sleeve_key_height=0.100
    [Inch]sleeve_key_thickness=0.100
    [Inch]p14=key_length
    [Inch]key_length=0.2
    [Inch]lox_manifold_height=0.279
    """)

shelf_length=1.4
lox_cv_length=2.125
side_thickness=0.5

with open("faceplate.exp", "w") as expr:
    expr.write(f"""
    // Version: 3
    [Inch]annulus_rad={Dann/2:.4f}
    [Inch]chamber_rad={chamber_diameter/2:.4f}
    [Inch]curvature_rad={chamber_diameter*7/9/2:.4f}
    [Degrees]slant_angle={2*180/pi*atan2(sleeve_max-sleeve_thickness,slot_height):.4f}
    [Inch]slant_length=0.5
    [Inch]side_thickness={side_thickness:.3f}
    [Inch]flange_rad=1.201
    [Inch]p1=annulus_rad
    [Degrees]p2=slant_angle
    [Inch]p3=slant_length
    [Inch]p4=chamber_rad
    [Inch]p5=curvature_rad
    [Inch]p6=1.5*annulus_rad
    [Inch]p10=side_thickness
    [Degrees]p7=0
    [Degrees]p8=360
    [Inch]p21=chamber_rad+1.5*annulus_rad
    [Inch]p22=flange_rad
    [Inch]p43=-side_thickness
    [Inch]p47=0.5
    [Inch]p57=chamber_rad+1.5*annulus_rad+flange_rad/2
    p175="5/8-11"
    p176="0.75"
    p177="Inch UNC"
    """)


with open("middleplate.exp", "w") as expr:
    gland = get_oring.get_gland_size()

    radial_rings = get_oring.min_inner(sleeve_max*2 + Dp)

    expr.write(f"""
    // Version: 3
    [Inch]oring_backing={gland['spacing']:.3f}
    [Inch]oring_width={gland['width']:.3f}
    [Inch]oring_depth={gland['depth']:.3f}
    [Inch]slot_height={slot_height:.3f}
    [Inch]key_height=0.100
    [Inch]flange_rad=1.201
    [Inch]middle_plate_thickness={2*slot_height+2*gland['width']+3*gland['spacing']+Dp:.3f}
    [Inch]curvature_rad={chamber_diameter*7/9/2:.4f}
    [Inch]slant_distance={radial_rings['ID']/2:.4f}
    [Inch]side_thickness={side_thickness:.3f}
    [Inch]shelf_length={shelf_length:.4f}
    [Inch]p5=slant_distance
    [Inch]p6=curvature_rad
    [Inch]p7={chamber_diameter/2:.4f}
    [Inch]p9=middle_plate_thickness
    [Inch]lox_cv_length=2.125
    [Inch]p13=shelf_length
    [Inch]lox_cv_height=0.25
    [Inch]p8={Dann/2*1.5:.4f}
    [Inch]p10=side_thickness
    [Inch]lox_cv_length={lox_cv_length:.3f}
    [Degrees]p14=0
    [Degrees]p15=360
    [Inch]p17=lox_cv_length+slant_distance+shelf_length/2
    p175="1/2-13"
    p176="0.75"
    p177="Inch UNC"
    p225={chamber_diameter/2+1.5*Dann/2:.4f}
    p268=flange_rad
    p267=p225+(flange_rad/2)
    p266=0.5
    p435="5/8-11"
    p436="0.75"
    p437="Inch UNC"
    """)

with open("topplate.exp", "w") as expr:
    expr.write(f"""
    [Inch]top_plate_thickness=0.75
    [Inch]shelf_length={shelf_length:.4f}
    [Inch]slant_distance={Dp:.3f}
    [Inch]p5=slant_distance
    [Inch]p7=2.125+shelf_length
    [Inch]p20=lox_cv_length+slant_distance+shelf_length/2
    p69="Inch UNC"
    p71="1/2-13"
    p73="0.75"
    """)



fig1, ax1 = plt.subplots(2, 3)
fig2, ax2 = plt.subplots(2, 2)
ax1[0,0].plot(throttles, chamber_pressures)
ax1[0,0].set_xlim(0.0, 1.0)
ax1[0,0].set_ylim(0.0, Pc)
ax1[0,0].set_xlabel("Throttle Level (ratio)")
ax1[0,0].set_ylabel("Chamber Pressure (psia)")
ax1[0,0].set_title("Chamber Pressure vs. Throttle")

ax1[1,0].plot(throttles, mass_flows)
ax1[1,0].set_xlim(0.0, 1.0)
ax1[1,0].set_ylim(0.0, Mo + Mf)
ax1[1,0].set_xlabel("Throttle Level (ratio)")
ax1[1,0].set_ylabel("Total Mass Flow (lbm / s)")
ax1[1,0].set_title("Total Mass Flow vs. Throttle")

ax1[0,1].plot(chamber_pressures, mass_flows)
ax1[0,1].set_xlim(np.min(chamber_pressures), Pc)
ax1[0,1].set_ylim(0.0, Mo + Mf)
ax1[0,1].set_xlabel("Chamber Pressure (psia)")
ax1[0,1].set_ylabel("Total Mass Flow (lbm / s)")
ax1[0,1].set_title("Total Mass Flow vs. Chamber Pressure")

ax2[0,0].plot(throttles, lox_areas)
ax2[0,0].plot(throttles, fuel_areas)
ax2[0,0].set_xlim(0.0, 1.0)
ax2[0,0].set_ylim(0.0, max(Aann, Aslot))
ax2[0,0].set_xlabel("Throttle Level (ratio)")
ax2[0,0].set_ylabel("Orifice Areas (sq.in)")
ax2[0,0].set_title("Orifice Areas vs. Throttle Level")
ax2[0,0].legend(["LOx", "RP1"])

ax1[1,1].plot(throttles, tmrs)
ax1[1,1].plot(throttles, lmrs)
ax1[1,1].set_xlim(0.0, 1.0)
ax1[1,1].set_ylim(0.0, 1.05 * max(np.max(tmrs), np.max(lmrs)))
ax1[1,1].set_xlabel("Throttle Level (ratio)")
ax1[1,1].set_ylabel("TMR and LMR (ratio)")
ax1[1,1].set_title("TMR and LMR vs. Throttle Level")
ax1[1,1].legend(["TMR", "LMR"])

ax2[0,1].plot(throttles, lox_stiffnesses)
ax2[0,1].plot(throttles, fuel_stiffnesses)
ax2[0,1].set_xlim(0.0, 1.0)
ax2[0,1].set_ylim(0.0, max(np.max(lox_stiffnesses), np.max(fuel_stiffnesses)))
ax2[0,1].set_xlabel("Throttle Level (ratio)")
ax2[0,1].set_ylabel("Stiffness (ratio)")
ax2[0,1].set_title("Stiffness vs. Throttle Level")
ax2[0,1].legend(["LOx", "RP1"])

ax1[0,2].plot(throttles, spray_angles)
ax1[0,2].set_xlim(0.0, 1.0)
ax1[0,2].set_ylim(np.min(spray_angles), np.max(spray_angles))
ax1[0,2].set_xlabel("Throttle Level (ratio)")
ax1[0,2].set_ylabel("Spray Angle (degrees)")
ax1[0,2].set_title("Spray Angle vs. Throttle Level")

ax1[1,2].plot(throttles, ispvacs)
ax1[1,2].plot(throttles, ispambs)
ax1[1,2].vlines(separation_point, 0.0, 1.05 * np.max(ispvacs), linestyles='dashed', label='Flow Separation')
ax1[1,2].set_xlim(0.0, 1.0)
ax1[1,2].set_ylim(0.0, 1.05 * np.max(ispvacs))
ax1[1,2].set_xlabel("Thrust Level (ratio)")
ax1[1,2].set_ylabel("ISP (s)")
ax1[1,2].set_title("Vacuum and Sea Level ISP vs. Throttle Level")
ax1[1,2].legend(["Vacuum", "Sea Level"])

ax2[1,0].plot(throttles, slotheights)
ax2[1,0].plot(throttles, annulargaps)
ax2[1,0].set_xlim(0.0, 1.0)
ax2[1,0].set_ylim(0.0, max(np.max(slotheights), np.max(annulargaps)))
ax2[1,0].vlines(min_throttle, 0.0, max(np.max(slotheights), np.max(annulargaps)), linestyles='dashed', label='Alignment Limit')
ax2[1,0].set_xlabel("Thrust Level (ratio)")
ax2[1,0].set_ylabel("Size (in)")
ax2[1,0].set_title("Slot Height and Annular Gap vs. Throttle Level")
ax2[1,0].legend(["Slot Height", "Annular Gap"])

ax2[1,1].plot(slotheights, annulargaps)
ax2[1,1].set_xlim(np.min(slotheights), np.max(slotheights))
ax2[1,1].set_ylim(0.0, np.max(annulargaps))
ax2[1,1].vlines(min_slotheight, 0.0, np.max(annulargaps), linestyles='dashed', label='Alignment Limit')
ax2[1,1].set_xlabel("Slot Height (in)")
ax2[1,1].set_ylabel("Annular Gap (in)")
ax2[1,1].set_title("Annular Gap vs. Slot Height")

plt.show()
