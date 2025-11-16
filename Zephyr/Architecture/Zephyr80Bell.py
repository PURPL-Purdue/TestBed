import math
import numpy as np

# === User-defined inputs ===
Dt_cm = 6.74       # Throat diameter in cm
De_cm = 15.93      # Exit diameter in cm
L_cm = 17.14  * .8     # Nozzle length in cm
theta_start_deg = 30.0
theta_end_deg = 7.0

# === Convert to meters ===
Dt = Dt_cm / 100.0
De = De_cm / 100.0
L = L_cm / 100.0

# === Radii ===
rt = Dt / 2
re = De / 2

# === Convert angles to slopes ===
slope_start = math.tan(math.radians(theta_start_deg))
slope_end = math.tan(math.radians(theta_end_deg))

# === Fit parabola: r(x) = ax^2 + bx + c ===
# Constraints:
# r(0) = rt → c = rt
# r'(0) = slope_start → b = slope_start
# r(L) = re → a*L^2 + b*L + c = re
# r'(L) = slope_end → 2a*L + b = slope_end

# Solve for a, b, c
a = (slope_end - slope_start) / (2 * L)
b = slope_start
c = rt

# === Midpoint and angles ===
x_mid = L / 2
r_mid = a * x_mid**2 + b * x_mid + c

slope_start_check = b
slope_end_check = 2 * a * L + b
angle_start_check = math.degrees(math.atan(slope_start_check))
angle_end_check = math.degrees(math.atan(slope_end_check))

# === Output in mm ===
print("\n=== Rao-Style Bell Nozzle Geometry (Constrained Parabola) ===")
print(f"Throat Diameter:     {Dt * 1000:.2f} mm")
print(f"Exit Diameter:       {De * 1000:.2f} mm")
print(f"Nozzle Length:       {L * 1000:.2f} mm")
print(f"Start Angle:         {angle_start_check:.2f}°")
print(f"End Angle:           {angle_end_check:.2f}°")
print(f"Midpoint Location:   {x_mid * 1000:.2f} mm")
print(f"Midpoint Radius:     {r_mid * 1000:.2f} mm")