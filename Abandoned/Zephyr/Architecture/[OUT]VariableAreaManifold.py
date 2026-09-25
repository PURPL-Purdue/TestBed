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

# Design constraints
area_min = 1e-4             # [m²]
area_max = 2.0e-2           # [m²]

# Iteration control
max_iter = 5000
meet_tol = 1e3             # [Pa] pressure continuity tolerance
mass_tol = 1e-6            # [m^3/s]
tolerance = 1e-6           # [m^3/s] mass balanced tolerance
alpha = 0.2                # slightly stronger relaxation

# ============================================================
# OUTLET SPECIFIC INPUTS
# ============================================================
Vman = 4                   # [m/s]  The velocity we want throughout the manifold (The target velocity)
Cd_branch = 0.65           #        Discharge efficiency factor, basically the efficiency of the flow from the inlet channels.
A_branch = 1.0e-4          # [m^2]  Area of each of the inlet channels (Same as the one for the inlet manifold)
K_join = 0.0               #        Pressure loss from when each channel joins the main manifold area.
dP_channel_guess_psi = 10  # [psi]  guessed channel -> collector delta p
# ============================================================
# INITIALIZATION
# ============================================================
psi_to_pa    = 6894.76            # Conversion ratio  
P_out_target = 750 * psi_to_pa    # Final pressure (required pressure)
P_chan       = 995 * psi_to_pa # MINIMUM pressure value from the inlet manifold

half_outlets = num_outlets // 2

# Start with slightly tapered manifold areas to avoid zero imbalance
area = np.full(half_outlets, Q_total / Vman)

# ============================================================
# DUAL-FLOW SIMULATION (0° → ±180°)
# ============================================================
pressure_forward = np.full(half_outlets + 1, P_out_target) #Starts from the outlet channel. 
pressure_backward = np.full(half_outlets + 1, P_out_target)
area_forward = np.copy(area)
area_backward = np.copy(area)

converged = False

for iteration in range(max_iter):
    # --- Forward half (0° → 180°)
    Q_running_f = 0.0
    for i in range(half_outlets):
        # Branch inflow into header at node i
        dp_branch = max(P_chan - pressure_forward[i], 0.0)
        q_in = Cd_branch  * A_branch * np.sqrt(2.0 * dp_branch / rho_lox)
        Q_running_f += q_in

        # tee/junction mixing loss
        A_header = np.clip(area_forward[i], area_min, area_max)
        D_header = np.sqrt(4.0 * A_header / np.pi)
        v_header = np.clip(Q_running_f / A_header, 1e-6, 1e4)
        dp_join = K_join * 0.5 * rho_lox * v_header**2
        pressure_forward[i] -= dp_join

        # propogating friction for the next node
        dp_fric = 0.5 * f * rho_lox * v_header**2 * (segment_length / D_header)
        pressure_forward[i+1] = pressure_forward[i] - dp_fric

        # manifold area above this outlet
        A_target = Q_running_f / max(Vman, 1e-6)
        area_forward[i] = np.clip((1 - alpha)*area_forward[i] + alpha*A_target, area_min, area_max)

    

    # --- Backward half (0° → -180°)
    Q_running_b = 0.0
    for i in range(half_outlets):

        dp_branch = max(P_chan - pressure_backward[i], 0.0)
        q_in = Cd_branch * A_branch * np.sqrt(2.0 * dp_branch / rho_lox)
        Q_running_b += q_in

        A_header = np.clip(area_backward[i], area_min, area_max)
        D_header = np.sqrt(4.0 * A_header / np.pi)
        v_header = np.clip(Q_running_b / A_header, 1e-6, 1e4)
        dp_join = K_join * 0.5 * rho_lox * v_header**2
        pressure_backward[i] -= dp_join
        
        dp_fric = 0.5 * f * rho_lox * v_header**2 * (segment_length / D_header)
        pressure_backward[i+1] = pressure_backward[i] - dp_fric

        A_target = Q_running_b / max(Vman, 1e-6)
        area_backward[i] = np.clip((1 - alpha)*area_backward[i] + alpha*A_target, area_min, area_max)

    # Combine both sides for balancing (optional)
    p_meet_f = pressure_forward[-1]
    p_meet_b = pressure_backward[-1]
    meet_err = abs(p_meet_f - p_meet_b)

    Q_collected = Q_running_f + Q_running_b
    if Q_collected > 0:
        scale = np.sqrt(np.clip(Q_total / Q_collected, 0.2, 5.0))
        A_branch *= np.clip(scale, 0.8, 1.25)
    mass_err = abs(Q_collected - Q_total)

    if (meet_err < meet_tol) and (mass_err < mass_tol):
        converged = True
        break

# ============================================================
# MERGE INTO FULL 360° RING (manifold taper)
# ============================================================
pressure_half = 0.5 * (pressure_forward[1:] + pressure_backward[1:])
area_half = 0.5 * (area_forward + area_backward)

# mirror to full ring
pressure_full = np.concatenate((pressure_half, pressure_half[::-1]))
area_full = np.concatenate((area_half, area_half[::-1]))
D_full = np.sqrt(4.0 * area_full / np.pi)

# header flow & velocity profile around the ring (using same branch law)
# build q_in per node for the half, then mirror
q_in_half = Cd_branch * A_branch * np.sqrt(2.0 * np.maximum(P_chan - pressure_half, 0.0) / rho_lox)
Q_profile_half = np.cumsum(q_in_half)
Q_profile_full = np.concatenate((Q_profile_half, Q_profile_half[::-1]))
v_header_full = Q_profile_full / np.clip(area_full, 1e-9, None)

# angles and plotting arrays
angles_full = np.linspace(0, 360, num_outlets, endpoint=False)
pressure_full_psi = pressure_full / psi_to_pa
std_dev = np.std(pressure_full_psi)

# ============================================================
# OUTPUT
# ============================================================
if converged:
    print(f"Converged after {iteration} iterations "
          f"(meet_err={meet_err:.3e} Pa, mass_err={mass_err:.3e} m^3/s).")
else:
    print("Warning: did not fully converge after max iterations "
          f"(meet_err={meet_err:.3e} Pa, mass_err={mass_err:.3e} m^3/s).")

print(f"{'Angle (deg)':>10s} | {'Area (cm^2)':>12s} | {'D (cm)':>8s} | {'Pressure (psi)':>14s} | {'Header v (m/s)':>14s}")
print("-" * 90)
for i in range(num_outlets):
    print(f"{angles_full[i]:10.1f} | {area_full[i]*1e4:12.4f} | {D_full[i]*100:8.4f} | {pressure_full_psi[i]:14.3f} | {v_header_full[i]:14.3f}")
print(f"\nOutlet pressure std dev around ring: {std_dev:.6f} psi")

# ============================================================
# PLOTS
# ============================================================
fig, ax = plt.subplots(3, 1, figsize=(8,12), sharex=True)

ax[0].plot(angles_full, pressure_full_psi, 'b-o', label='Header Pressure (psi)')
ax[0].set_ylabel("Pressure (psi)"); ax[0].grid(True); ax[0].legend()

ax[1].plot(angles_full, D_full*100, 'r-o', label='Manifold Diameter (cm)')
ax[1].set_ylabel("Diameter (cm)"); ax[1].grid(True); ax[1].legend()

ax[2].plot(angles_full, v_header_full, 'g-o', label='Header Velocity (m/s)')
ax[2].set_xlabel("Outlet angle (deg)")
ax[2].set_ylabel("Velocity (m/s)")
ax[2].grid(True); ax[2].legend()

plt.suptitle("LOX Collector Ring — Pressure, Diameter, and Header Velocity")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()