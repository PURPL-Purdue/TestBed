from math import pi, sqrt, isnan
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

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


stiffness = 0.2
sleeve_thickness = 0.100
pintle_wall = 0.050
blockage_factor = 0.85

Cdr = 0.61
Cda = 0.98

Pa = 14.695

max_thrust = 5000
Pc = 400
Mf = 3.0736 / lb_to_kg
Mo = 6.7616 / lb_to_kg
rhof = 768.09 * (in_to_m ** 3) / lb_to_kg
rhoo = 1205 / lb_to_kg * (in_to_m ** 3)

expansion_ratio = 4.897
chamber_diameter = (0.1665) / in_to_m
exit_area = ((0.1843 / 2) / in_to_m)**2 * pi

Dp = 0.2 * chamber_diameter
Dbore = Dp - 2 * pintle_wall

Aslot = Mo / sqrt(2 * rhoo * stiffness * Pc * lbf_to_lbins) / Cdr
print(f"Slot Area: {Aslot:.6f}")
slot_height = Aslot / (pi * Dbore * blockage_factor)

Rmomentum = Mo ** 2 / rhoo / Aslot

Aann = Mf ** 2 / Rmomentum / rhof
Dann = sqrt((4 * Aann) / pi + Dp ** 2)

Pdf = (Mf ** 2 / (2 * Cda ** 2 * Aann ** 2 * rhof)) / lbf_to_lbins
ox_feed = (1 + stiffness) * Pc
fuel_feed = Pdf + Pc

print(f"Pintle Tip Outer Diameter (in): {Dp:.3f}")
print(f"Pintle Tip Inner Diameter (in): {Dbore:.3f}")
print(f"Pintle Slot Height (in): {slot_height:.3f}")
print(f"Annulus Outer Diameter (in): {Dann:.3f}")
print(f"Annulur Gap (in): {Dann - Dp:.3f}")
print(f"Fuel Stiffness (%): {100 * Pdf/Pc:.1f}")
print(f"LOx Feed Pressure (psia): {ox_feed: .0f}")
print(f"RP1 Feed Pressure (psia): {fuel_feed: .0f}")

cea = CEA_Obj(oxName='LOX', fuelName='RP1')

throttles = np.linspace(0.99, 0.01, 99)
chamber_pressures = np.empty_like(throttles)
for i in range(len(throttles)):
    throttle = throttles[i]

    desired = throttle * max_thrust / exit_area + Pa
    print(f"Starting T={throttle*100:.0f}")
    
    converged = False
    chamber_pressure = Pc * throttle
    while not converged:
        actual = get_throttle_target(chamber_pressure, Mo / Mf, expansion_ratio)
        actual -= desired
        if abs(actual / desired) < 0.001:
            converged = True
        else:
            epsilon = get_throttle_target(chamber_pressure + 0.01, Mo / Mf, expansion_ratio)
            epsilon -= desired

            slope = (epsilon - actual) / 0.01
            chamber_pressure -= actual / slope

    chamber_pressures[i] = chamber_pressure

print(f"Throttles: {throttles}")
print(f"Pressures: {chamber_pressures}")

mass_flows = np.empty_like(throttles)
for i in range(len(throttles)):
    rhoe = cea.get_Densities(chamber_pressures[i], Mo / Mf, expansion_ratio)[2] / (12 ** 3)
    mach = cea.get_MachNumber(chamber_pressures[i], Mo / Mf, expansion_ratio)
    sonic = cea.get_SonicVelocities(chamber_pressures[i], Mo / Mf, expansion_ratio)[2] * 12

    mass_flows[i] = rhoe * mach * sonic * exit_area

print(f"Mass Flows: {mass_flows}")

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

print(f"Lox Areas: {lox_areas}")
print(f"Fuel Areas: {fuel_areas}")

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


fig, ax = plt.subplots(2, 4)
ax[0,0].plot(throttles, chamber_pressures)
ax[0,0].set_xlim(0.0, 1.0)
ax[0,0].set_ylim(0.0, Pc)
ax[0,0].set_xlabel("Throttle Level (ratio)")
ax[0,0].set_ylabel("Chamber Pressure (psia)")
ax[0,0].set_title("Chamber Pressure vs. Throttle")

ax[1,0].plot(throttles, mass_flows)
ax[1,0].set_xlim(0.0, 1.0)
ax[1,0].set_ylim(0.0, Mo + Mf)
ax[1,0].set_xlabel("Throttle Level (ratio)")
ax[1,0].set_ylabel("Total Mass Flow (lbm / s)")
ax[1,0].set_title("Total Mass Flow vs. Throttle")

ax[0,1].plot(chamber_pressures, mass_flows)
ax[0,1].set_xlim(np.min(chamber_pressures), Pc)
ax[0,1].set_ylim(0.0, Mo + Mf)
ax[0,1].set_xlabel("Chamber Pressure (psia)")
ax[0,1].set_ylabel("Total Mass Flow (lbm / s)")
ax[0,1].set_title("Total Mass Flow vs. Chamber Pressure")

ax[1,1].plot(throttles, lox_areas)
ax[1,1].plot(throttles, fuel_areas)
ax[1,1].set_xlim(0.0, 1.0)
ax[1,1].set_ylim(0.0, max(Aann, Aslot))
ax[1,1].set_xlabel("Throttle Level (ratio)")
ax[1,1].set_ylabel("Orifice Areas (sq.in)")
ax[1,1].set_title("Orifice Areas vs. Throttle Level")
ax[1,1].legend(["LOx", "RP1"])

ax[0,2].plot(throttles, tmrs)
ax[0,2].plot(throttles, lmrs)
ax[0,2].set_xlim(0.0, 1.0)
ax[0,2].set_ylim(0.0, 1.05 * max(np.max(tmrs), np.max(lmrs)))
ax[0,2].set_xlabel("Throttle Level (ratio)")
ax[0,2].set_ylabel("TMR and LMR (ratio)")
ax[0,2].set_title("TMR and LMR vs. Throttle Level")
ax[0,2].legend(["TMR", "LMR"])

ax[1,2].plot(throttles, lox_stiffnesses)
ax[1,2].plot(throttles, fuel_stiffnesses)
ax[1,2].set_xlim(0.0, 1.0)
ax[1,2].set_ylim(0.0, max(np.max(lox_stiffnesses), np.max(fuel_stiffnesses)))
ax[1,2].set_xlabel("Throttle Level (ratio)")
ax[1,2].set_ylabel("Stiffness (ratio)")
ax[1,2].set_title("Stiffness vs. Throttle Level")
ax[1,2].legend(["LOx", "RP1"])

ax[0,3].plot(throttles, spray_angles)
ax[0,3].set_xlim(0.0, 1.0)
ax[0,3].set_ylim(np.min(spray_angles), np.max(spray_angles))
ax[0,3].set_xlabel("Throttle Level (ratio)")
ax[0,3].set_ylabel("Spray Angle (degrees)")
ax[0,3].set_title("Spray Angle vs. Throttle Level")

ax[1,3].plot(throttles, ispvacs)
ax[1,3].plot(throttles, ispambs)
ax[1,3].vlines(separation_point, 0.0, 1.05 * np.max(ispvacs), linestyles='dashed', label='Flow Separation')
ax[1,3].set_xlim(0.0, 1.0)
ax[1,3].set_ylim(0.0, 1.05 * np.max(ispvacs))
ax[1,3].set_xlabel("Thrust Level (ratio)")
ax[1,3].set_ylabel("ISP (s)")
ax[1,3].set_title("Vacuum and Sea Level ISP vs. Throttle Level")
ax[1,3].legend(["Vacuum", "Sea Level"])
plt.show()
