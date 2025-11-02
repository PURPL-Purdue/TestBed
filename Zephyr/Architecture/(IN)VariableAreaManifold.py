import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER CONFIGURATION
# ============================================================
num_outlets = 40           # total outlets spaced around 360°
ring_radius = 0.5           # [m]
segment_length = (2 * np.pi * ring_radius) / num_outlets  # segment length [m]

# Fluid properties (for LOX)
rho_lox = 1140.0           # density [kg/m^3]
mu_lox = 0.000197          # viscosity [Pa·s]
f = 0.02                   # Darcy friction factor (approx.)

# Flow and pressure
Q_total = 4.6979e-03       # total inlet flow [m³/s]
P_in = 1000.0              # inlet pressure [psi]
P_out_target = 995.0       # desired outlet pressure [psi]

# Design constraints
area_min = 1.0e-4          # [m²]
area_max = 1.0e-2          # [m²]

# Iteration control
max_iter = 5000
tolerance = 1e-6
relax = 0.5   # slightly stronger relaxation

# --- Constant outlet area
A_outlet = 1.0e-4  # [m²]

# ============================================================
# INITIALIZATION
# ============================================================
psi_to_pa = 6894.76
P_in *= psi_to_pa
P_out_target *= psi_to_pa

half_outlets = num_outlets // 2
Q_half = Q_total / 2
Q_per_outlet = Q_total / num_outlets

# Start with slightly tapered manifold areas to avoid zero imbalance
area = np.linspace(4.2e-4, 3.8e-4, half_outlets)

# ============================================================
# DUAL-FLOW SIMULATION (0° → ±180°)
# ============================================================
pressure_forward = np.full(half_outlets + 1, P_in)
pressure_backward = np.full(half_outlets + 1, P_in)
area_forward = np.copy(area)
area_backward = np.copy(area)

for iteration in range(max_iter):
    # --- Forward half (0° → 180°)
    Q_down = Q_half
    for i in range(half_outlets):
        v = np.clip(Q_down / A_outlet, 1e-6, 1e4)  # velocity using constant outlet area
        dp_fric = 0.5 * f * rho_lox * v**2 * (segment_length / np.sqrt(4 * A_outlet / np.pi))
        pressure_forward[i+1] = pressure_forward[i] - dp_fric

        # manifold area above this outlet
        A_manifold = Q_down / np.clip(v, 1e-6, 1e4)
        area_forward[i] = A_manifold

        Q_down -= Q_per_outlet

    # --- Backward half (0° → -180°)
    Q_down = Q_half
    for i in range(half_outlets):
        v = np.clip(Q_down / A_outlet, 1e-6, 1e4)
        dp_fric = 0.5 * f * rho_lox * v**2 * (segment_length / np.sqrt(4 * A_outlet / np.pi))
        pressure_backward[i+1] = pressure_backward[i] - dp_fric

        A_manifold = Q_down / np.clip(v, 1e-6, 1e4)
        area_backward[i] = A_manifold

        Q_down -= Q_per_outlet

    # Combine both sides for balancing (optional)
    P_local = 0.5 * (pressure_forward[1:] + pressure_backward[1:])
    dp = np.maximum(P_local - P_out_target, 0.0)
    capacity = np.sqrt(dp + 1e-12) * np.sqrt(area_forward + area_min)
    if np.all(capacity <= 0):
        capacity = np.ones_like(capacity)
    Q_out = capacity / np.sum(capacity) * Q_half

    imbalance = (Q_out - Q_per_outlet) / (Q_per_outlet + 1e-15)
    area_forward *= np.clip(1.0 - relax * imbalance, 0.5, 1.5)
    area_backward *= np.clip(1.0 - relax * imbalance, 0.5, 1.5)

    # --- Enforce decreasing taper from inlet to far end
    taper_factor = np.linspace(1.0, 0.6, half_outlets)
    area_forward = np.clip(area_forward * taper_factor, area_min, area_max)
    area_backward = np.clip(area_backward * taper_factor, area_min, area_max)

    if np.max(np.abs(imbalance)) < tolerance:
        converged = True
        break

pressure_half = 0.5 * (pressure_forward[1:] + pressure_backward[1:])
area_half = 0.5 * (area_forward + area_backward)

# ============================================================
# MERGE INTO FULL 360° RING (manifold taper)
# ============================================================
angles_full = np.linspace(0, 360, num_outlets)
# Manifold area above each outlet: taper high at 0°, low at 180°, back high at 360°
taper_factor_full = np.linspace(1.0, 0.6, num_outlets//2)
taper_factor_full = np.concatenate((taper_factor_full, taper_factor_full[::-1]))
area_full = np.clip(taper_factor_full * area_max, area_min, area_max)
D_full = np.sqrt(4 * area_full / np.pi)

pressure_full = np.concatenate((pressure_half, pressure_half[::-1]))

# Outlet velocity uses constant outlet area
velocity_full = Q_per_outlet / A_outlet
pressure_full_psi = pressure_full / psi_to_pa
std_dev = np.std(pressure_full_psi)

# ============================================================
# OUTPUT
# ============================================================
if converged:
    print(f"Converged after {iteration} iterations (max rel imbalance {np.max(np.abs(imbalance)):.3e}).\n")
else:
    print("Warning: did not fully converge after max iterations.\n")

print(f"{'Outlet Angle':>10s} | {'Area (cm^2)':>12s} | {'D (cm)':>8s} | {'Pressure (psi)':>14s} | {'Velocity (m/s)':>14s}")
print("-" * 90)
for i in range(num_outlets):
    print(f"{angles_full[i]:10.1f} | {area_full[i]*1e4:12.4f} | {D_full[i]*100:8.4f} | {pressure_full_psi[i]:14.3f} | {velocity_full:14.3f}")
print(f"\nOutlet pressure standard deviation: {std_dev:.6f} psi")

# ============================================================
# PLOTS (ALL IN ONE FIGURE)
# ============================================================
fig, ax = plt.subplots(3, 1, figsize=(8,12), sharex=True)

ax[0].plot(angles_full, pressure_full_psi, 'b-o', label='Pressure (psi)')
ax[0].set_ylabel("Pressure (psi)")
ax[0].grid(True)
ax[0].legend()

ax[1].plot(angles_full, D_full*100, 'r-o', label='Manifold Diameter (cm)')
ax[1].set_ylabel("Diameter (cm)")
ax[1].grid(True)
ax[1].legend()

ax[2].plot(angles_full, [velocity_full]*num_outlets, 'g-o', label='Outlet Velocity (m/s)')
ax[2].set_xlabel("Outlet angle (deg)")
ax[2].set_ylabel("Velocity (m/s)")
ax[2].grid(True)
ax[2].legend()

plt.suptitle("LOX Ring Outlet Flow Characteristics (Manifold Area Displayed)")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()
