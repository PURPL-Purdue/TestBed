import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# ================= User config =================
m_ox_total = 160.81         # kg (total oxidizer used in 30 s hotfire)
hotfire_time = 30.0         # s
mdot_lox = m_ox_total / hotfire_time   # kg/s (LOX mass flow)

rho_lox = 1141.0            # kg/m^3
f = 0.02
Q_in = mdot_lox / rho_lox   # m^3/s total volumetric LOX flow into manifold

num_outlets = 12
R = 0.17                    # m: inner radius of ring manifold

# Inlet pressure
P_in_psi = 1000.0           # psi
P_in = P_in_psi * 6894.757  # Pa
P_out_target = P_in / num_outlets  # Pa, ideal per-channel outlet

# ================= Outlet positions =================
theta_outlets = np.linspace(0, 360, num_outlets + 1)  # include 360° to close ring
Q_per_outlet = Q_in / num_outlets
segment_length = 2 * np.pi * R / num_outlets  # meters

# Initialize arrays
area_cm2 = np.zeros(num_outlets + 1)
diam_cm = np.zeros(num_outlets + 1)
velocity_m_s = np.zeros(num_outlets + 1)
pressure_Pa = np.zeros(num_outlets + 1)
pressure_Pa[0] = P_in

# ================= Compute areas segment by segment =================
for i in range(num_outlets // 2 + 1):
    Q_downstream = Q_per_outlet * (num_outlets // 2 - i + 1)

    def func(A):
        dP_actual = 0.5 * f * rho_lox * (Q_downstream / A)**2 * segment_length
        # evenly divide total drop along half-ring segments
        deltaP_cum = max(pressure_Pa[0] - P_out_target, 1e-3)
        return dP_actual - (deltaP_cum / (num_outlets // 2))

    # Initial guess avoids NaN
    A_guess = Q_downstream * np.sqrt(f * rho_lox * segment_length / (2 * max(1e-6, pressure_Pa[0] - P_out_target)))
    A_min = max(A_guess * 0.1, 1e-12)
    A_max = A_guess * 10

    try:
        A_sol = brentq(func, A_min, A_max)
    except ValueError:
        A_sol = A_guess  # fallback if solver fails

    area_cm2[i] = A_sol * 1e4
    diam_cm[i] = 2 * np.sqrt(A_sol / np.pi) * 100
    velocity_m_s[i] = Q_per_outlet / A_sol

    dP_actual = 0.5 * f * rho_lox * (Q_downstream / A_sol)**2 * segment_length
    pressure_Pa[i+1] = pressure_Pa[i] - dP_actual

# Symmetric side of the ring
for i in range(num_outlets // 2 + 1, num_outlets + 1):
    mirror_idx = num_outlets - i
    area_cm2[i] = area_cm2[mirror_idx]
    diam_cm[i] = diam_cm[mirror_idx]
    velocity_m_s[i] = velocity_m_s[mirror_idx]
    pressure_Pa[i] = pressure_Pa[mirror_idx]

# Convert pressures to psi for plotting
pressure_psi = pressure_Pa / 6894.757

# ================= Print table =================
print("Outlet Angle | Q (m^3/s) | Area (cm^2) | D (cm) | Pressure (psi) | Velocity (m/s)")
print("-"*90)
for i, angle in enumerate(theta_outlets):
    print(f"{angle:12.1f} | {Q_per_outlet:10.6e} | {area_cm2[i]:10.4f} | {diam_cm[i]:10.4f} | {pressure_psi[i]:10.3f} | {velocity_m_s[i]:10.2f}")

# ================= Visualization =================
plt.figure(figsize=(12, 6))

plt.subplot(3,1,1)
plt.plot(theta_outlets, diam_cm, 'o-', label='Diameter (cm)')
plt.ylabel('Diameter (cm)')
plt.grid(True)
plt.legend()

plt.subplot(3,1,2)
plt.plot(theta_outlets, velocity_m_s, 's-', label='Velocity (m/s)', color='orange')
plt.ylabel('Velocity (m/s)')
plt.grid(True)
plt.legend()

plt.subplot(3,1,3)
plt.plot(theta_outlets, pressure_psi, 'd-', label='Pressure (psi)', color='green')
plt.xlabel('Outlet Angle (deg)')
plt.ylabel('Pressure (psi)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
