import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# USER CONFIGURATION
# ============================================================
num_outlets = 40          # number of INPUT pipes feeding the ring (branches)
ring_radius = 0.08325      # [m]

# Fluid properties (RP-1)
rho = 810.0                # density [kg/m^3]
f = 0.02                   # Darcy friction factor (assumed constant)

# Total mass flow coming IN from the 40 channels
m_dot_total = 3.0736       # [kg/s]
Q_total = m_dot_total / rho   # [m^3/s]

# Branch (input pipe) properties
Cd_branch = 0.65           # discharge coefficient for each inlet branch
A_branch = 2.0e-5          # [m^2] fixed area of each inlet pipe

# Channel (upstream) pressure feeding the branches
P_chan_psi = 954.43        # [psi]
psi_to_pa = 6894.76
P_chan = P_chan_psi * psi_to_pa   # [Pa]

# Outlet pipe target velocity (this is the single OUTPUT pipe)
V_out = 5.0                # [m/s] desired outlet velocity

SAFETY_FACTOR = 1.5        # factor to oversize header vs outlet area

free_zone_deg = 8.0        # total free zone width (e.g. 8° ⇒ -4°..+4°)
half_free = free_zone_deg / 2.0

# ============================================================
# CONTINUITY: BRANCH FLOWS AND OUTLET AREA
# ============================================================
# Each of the 40 INPUT pipes carries equal flow
Q_branch = Q_total / num_outlets          # [m^3/s]
V_branch = Q_branch / A_branch            # [m/s] velocity in each input pipe

# Outlet pipe area required to get the desired outlet velocity
A_exit = Q_total / V_out                  # [m^2]
D_exit = np.sqrt(4.0 * A_exit / np.pi)    # [m]

# Manifold header area (collector ring)
A_header = SAFETY_FACTOR * A_exit
D_header = np.sqrt(4.0 * A_header / np.pi)

# ============================================================
# BASELINE MANIFOLD PRESSURE FROM ORIFICE EQUATION (NO FRICTION)
# ============================================================

deltaP_branch = 0.5 * rho * (Q_branch / (Cd_branch * A_branch))**2  # [Pa]
P_man = P_chan - deltaP_branch                                      # [Pa]
P_man_psi = P_man / psi_to_pa

# ============================================================
# RING GEOMETRY AND FLOW PROFILE (COLLECTOR)
# ============================================================
half_outlets = num_outlets // 2  # 20 per side

# Each branch occupies an equal angular spacing in the non-free arc
arc_with_branches_deg = 360.0 - free_zone_deg
dtheta_deg = arc_with_branches_deg / num_outlets  # degrees per branch spacing
segment_length = ring_radius * np.deg2rad(dtheta_deg)  # [m] between branches

# Angles for branches:
# Negative side: far (-180+half_free) -> near (-half_free)
angles_neg = np.linspace(-180.0 + half_free, -half_free, half_outlets)
# Positive side: near (+half_free) -> far (+180-half_free)
angles_pos = np.linspace(half_free, 180.0 - half_free, half_outlets)

angles_full = np.concatenate((angles_neg, angles_pos))

# Flow in header on each side (collector behavior):
# From far end to near outlet, flow builds from Q_branch to Q_half:
Q_half = half_outlets * Q_branch  # total per side at the outlet junction

# For friction integration it's convenient to have flows from FAR -> NEAR:
Q_neg_far_to_near = np.linspace(Q_branch, Q_half, half_outlets)  # negative side
Q_pos_far_to_near = np.linspace(Q_branch, Q_half, half_outlets)  # positive side

# For display, we want:
#  - negative side: far -> near (as angles_neg is defined)
#  - positive side: near -> far, so we reverse the "far_to_near" array
Q_neg_display = Q_neg_far_to_near                          # far -> near
Q_pos_display = Q_pos_far_to_near[::-1]                    # near -> far

Q_profile_full = np.concatenate((Q_neg_display, Q_pos_display))  # [m^3/s]

# Header velocity at each branch location (for display)
V_header_full = Q_profile_full / A_header                  # [m/s]

# ============================================================
# FRICTION ESTIMATIONS ALONG THE RING (DARCY–WEISBACH)
# ============================================================

v_neg_seg = Q_neg_far_to_near / A_header  # [m/s]
dp_neg_seg = f * (segment_length / D_header) * 0.5 * rho * v_neg_seg**2  # [Pa]

v_pos_seg = Q_pos_far_to_near / A_header
dp_pos_seg = f * (segment_length / D_header) * 0.5 * rho * v_pos_seg**2  # [Pa]

# Pressures at branch locations:
# Let P_far_neg = P_man, then pressure drops toward the outlet along each side.
P_neg = np.zeros(half_outlets)
P_pos = np.zeros(half_outlets)

# Negative side: far (index 0) -> near (index N-1)
P_neg[0] = P_man
for j in range(1, half_outlets):
    P_neg[j] = P_neg[j-1] - dp_neg_seg[j-1]

# Positive side: far (index 0) -> near (index N-1)
P_pos[0] = P_man
for j in range(1, half_outlets):
    P_pos[j] = P_pos[j-1] - dp_pos_seg[j-1]

pressure_neg_display = P_neg                         # far -> near
pressure_pos_display = P_pos[::-1]                   # near -> far

pressure_full = np.concatenate((pressure_neg_display, pressure_pos_display))
pressure_full_psi = pressure_full / psi_to_pa

# Estimate far-end and near-outlet pressures for text ΔP
P_far_avg = 0.5 * (P_neg[0] + P_pos[0])          # should be ~P_man
P_near_out_avg = 0.5 * (P_neg[-1] + P_pos[-1])   # near outlet on both sides
deltaP_fric = P_far_avg - P_near_out_avg         # [Pa]
deltaP_fric_psi = deltaP_fric / psi_to_pa        # [psi]

# ============================================================
# MANIFOLD FILL TIME
# ============================================================
L_ring = 2 * np.pi * ring_radius               # total ring length [m]
V_manifold = A_header * L_ring                 # volume of the manifold [m^3]
t_fill = V_manifold / Q_total                  # fill time [s]

print("\nManifold Fill Time:")
print(f"  Ring length (L_ring):      {L_ring:.3f} m")
print(f"  Manifold volume:           {V_manifold:.6e} m^3")
print(f"  Fill time (t_fill):        {t_fill*1000:.3f} ms")

# ============================================================
# OUTPUTS
# ============================================================
print("\n=== Collector Manifold with Darcy–Weisbach Friction (40 IN → 1 OUT) ===")
print(f"Upstream channel pressure (P_chan):       {P_chan_psi:.3f} psi")
print(f"Baseline manifold pressure (no friction): {P_man_psi:.3f} psi")
print(f"Estimated friction pressure loss (far → outlet): {deltaP_fric_psi:.6f} psi ({deltaP_fric:.3f} Pa)")
print()
print(f"Total mass flow (m_dot_total):            {m_dot_total:.6f} kg/s")
print(f"Total volume flow (Q_total):              {Q_total:.6e} m^3/s")
print()
print(f"Each input pipe flow (Q_branch):          {Q_branch:.6e} m^3/s")
print(f"Each input pipe velocity:                 {V_branch:.3f} m/s")
print()
print("Outlet pipe sizing from desired velocity:")
print(f"  Desired outlet velocity (V_out):        {V_out:.3f} m/s")
print(f"  Outlet cross-sectional area (A_exit):   {A_exit:.6e} m^2")
print(f"  Outlet diameter (D_exit):               {D_exit*1000:.3f} mm")
print()
print("Header (ring manifold) geometry:")
print(f"  Header cross-sectional area (A_header): {A_header:.6e} m^2")
print(f"  Header diameter (D_header):             {D_header*1000:.3f} mm")
print()
print(f"Regen Channel free zone:                  ±{half_free:.1f}° around 0°")

# ============================================================
# PLOTS OF PRESSURE + HEADER VELOCITY
# ============================================================

# Convert pressure to absolute psi
pressure_full_abs_psi = pressure_full / psi_to_pa  # [psi]

# Compute mean pressure to set a nice axis window
P_mean = np.mean(pressure_full_abs_psi)
y_lo = P_mean - 0.5  # show ±0.5 psi around the mean
y_hi = P_mean + 0.5

fig, axs = plt.subplots(2, 1, figsize=(10, 9))

# ------------------------------------------------------------
# 1) PRESSURE PLOT
# ------------------------------------------------------------
axs[0].plot(angles_full, pressure_full_abs_psi, 'b-o', linewidth=1.6, markersize=4)
axs[0].set_title("Pressure Around Collector Manifold", fontsize=16, fontweight='bold')
axs[0].set_ylabel("Pressure (psi)", fontsize=13)
axs[0].set_xlabel(f"Angle Around Ring (0° = Outlet)", fontsize=13)
axs[0].grid(True, linestyle='--', alpha=0.6)
axs[0].tick_params(labelsize=11)
axs[0].set_ylim(y_lo, y_hi)  # force visible absolute scale

# ------------------------------------------------------------
# 2) HEADER VELOCITY PLOT
# ------------------------------------------------------------
axs[1].plot(angles_full, V_header_full, 'g-o', linewidth=1.6, markersize=4)
axs[1].set_title("Flow Velocity Around Collector Manifold", fontsize=16, fontweight='bold')
axs[1].set_ylabel("Velocity (m/s)", fontsize=13)
axs[1].set_xlabel(f"Angle Around Ring (0° = Outlet)", fontsize=13)
axs[1].grid(True, linestyle='--', alpha=0.6)
axs[1].tick_params(labelsize=11)

plt.tight_layout()
plt.show()