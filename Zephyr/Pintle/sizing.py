from math import pi, sqrt
in_to_m = 0.0254
lb_to_kg = 0.4536
N_to_lbf = 0.2248
lbf_to_lbins = 386.09

stiffness = 0.2
sleeve_thickness = 0.100
pintle_wall = 0.050
blockage_factor = 0.85

Cdr = 0.61
Cda = 0.98

max_thrust = 5000
Pc = 400
Mf = 3.0736 / lb_to_kg
Mo = 6.7616 / lb_to_kg
rhof = 768.09 * (in_to_m ** 3) / lb_to_kg
rhoo = 1205 / lb_to_kg * (in_to_m ** 3)

chamber_diameter = 0.1665 / in_to_m
throat_area = pi * ((0.0833 / in_to_m) / 2) ** 2

thrust_coefficient = max_thrust / Pc / throat_area

Dp = 0.2 * chamber_diameter
Dbore = Dp - 2 * pintle_wall

Aslot = Mf / sqrt(2 * rhof * stiffness * Pc * lbf_to_lbins) / Cdr
print(f"Slot Area: {Aslot:.6f}")
slot_height = Aslot / (pi * Dbore * blockage_factor)

Rmomentum = Mf ** 2 / rhof / Aslot

Aann = Mo ** 2 / Rmomentum / rhoo
Dann = sqrt((4 * Aann) / pi + Dp ** 2)

Pdox = (Mo ** 2 / (2 * Cda ** 2 * Aann ** 2 * rhoo)) / lbf_to_lbins
stiffness_ox = Pdox / Pc

print(f"Pintle Tip Outer Diameter (in): {Dp:.3f}")
print(f"Pintle Tip Inner Diameter (in): {Dbore:.3f}")
print(f"Pintle Slot Height (in): {slot_height:.3f}")
print(f"Annulus Outer Diameter (in): {Dann:.3f}")
print(f"Annulur Gap (in): {Dann - Dp:.3f}")
print(f"LOx Pressure Drop (psi): {Pdox:.3f}")
