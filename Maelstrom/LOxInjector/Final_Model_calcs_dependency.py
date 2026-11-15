import numpy as np
import matplotlib.pyplot as plt
import FILM_COOLING_PURPL_Testbed_Maelstrom_LOx_Variant_Impinging_Jet_Injector_Calcs as film_cooling


lox_mass_flow_rate = film_cooling.mdot_LOx / film_cooling.num_LOx_orifices
kerosene_mass_flow_rate_per_outer_jet = (
    film_cooling.nonFC_mdot_Kerosene / film_cooling.num_Kerosene_orifices
)

lox_orifice_diameter = film_cooling.LOx_orifice_diameter
kerosene_orifice_diameter = film_cooling.Kerosene_orifice_diameter

lox_density = film_cooling.density_LOx
# Surface tensions are not calculated by the film-cooling script.
lox_surface_tension = 0.013
kerosene_density = film_cooling.density_Kerosene
kerosene_surface_tension = 0.024

# The collision angle is also not calculated by the film-cooling script.
half_collision_angle = np.radians(30.0)

lox_orifice_area = np.pi * lox_orifice_diameter**2 / 4.0
kerosene_orifice_area = np.pi * kerosene_orifice_diameter**2 / 4.0

lox_jet_velocity = lox_mass_flow_rate / (lox_density * lox_orifice_area)
kerosene_jet_velocity = kerosene_mass_flow_rate_per_outer_jet / (kerosene_density * kerosene_orifice_area)

# --- 2. HELPER FUNCTIONS ---
def momentum_rate_jets(orifice_diameter, density, jet_velocity, state, alpha):
    if state == 1:
        return np.pi * (1.0 / 4.0) * (orifice_diameter ** 2) * density * (jet_velocity ** 2)
    elif state == 2:
        return 2.0 * np.pi * (1.0 / 4.0) * (orifice_diameter ** 2) * density * (jet_velocity ** 2) * np.cos(alpha)
    return 0.0

def mass_flow_jets(orifice_diameter, density, jet_velocity, state, alpha):
    if state == 1:
        return np.pi * (1.0 / 4.0) * (orifice_diameter ** 2) * density * jet_velocity
    return 0.0

def mixture_density_calculation(kerose_dens, lox_dens, d_kero, d_lox, v_kero, v_lox):
    Q1 = np.pi * (d_kero ** 2) / 4.0 * v_kero
    Q2 = np.pi * (d_lox ** 2) / 4.0 * v_lox
    return ((2.0 * kerose_dens * Q1) + (lox_dens * Q2)) / (2.0 * Q1 + Q2)

def mixture_sheet_velocity_calculation(kerose_dens, lox_dens, d_kero, d_lox, v_kero, v_lox):
    m1 = mass_flow_jets(d_kero, kerose_dens, v_kero, 1, 0)
    m2 = mass_flow_jets(d_lox, lox_dens, v_lox, 1, 0)
    return ((2.0 * m1 * v_kero) + (m2 * v_lox)) / (2.0 * m1 + m2)

def mixture_surface_tension_calculation(sigma_kero, sigma_lox, d_kero, d_lox, v_kero, v_lox):
    Q1 = np.pi * (d_kero ** 2) / 4.0 * v_kero
    Q2 = np.pi * (d_lox ** 2) / 4.0 * v_lox
    return ((2.0 * sigma_kero * Q1) + (sigma_lox * Q2)) / (2.0 * Q1 + Q2)

# --- 3. MIXTURE PROPERTIES & SHEET COEFFICIENTS ---
rho_s = mixture_density_calculation(kerosene_density, lox_density, kerosene_orifice_diameter, lox_orifice_diameter, kerosene_jet_velocity, lox_jet_velocity)
u_s = mixture_sheet_velocity_calculation(kerosene_density, lox_density, kerosene_orifice_diameter, lox_orifice_diameter, kerosene_jet_velocity, lox_jet_velocity)
sigma_s = mixture_surface_tension_calculation(kerosene_surface_tension, lox_surface_tension, kerosene_orifice_diameter, lox_orifice_diameter, kerosene_jet_velocity, lox_jet_velocity)

m_kero_single = mass_flow_jets(kerosene_orifice_diameter, kerosene_density, kerosene_jet_velocity, 1, 0)
m_lox_single = mass_flow_jets(lox_orifice_diameter, lox_density, lox_jet_velocity, 1, 0)
LHS_Eq1 = (2.0 * m_kero_single) + m_lox_single

p_kero_axial = momentum_rate_jets(kerosene_orifice_diameter, kerosene_density, kerosene_jet_velocity, 2, half_collision_angle)
p_lox_axial = momentum_rate_jets(lox_orifice_diameter, lox_density, lox_jet_velocity, 1, 0)
LHS_Eq2 = p_kero_axial + p_lox_axial

f = LHS_Eq2 / (LHS_Eq1 * u_s)
if not 0.0 <= f < 1.0:
    raise ValueError(f"The momentum ratio f must be in [0, 1); got {f:.6g}.")

d_j = lox_orifice_diameter
We = (rho_s * (u_s**2) * d_j) / sigma_s

def g_factor(f_val):
    return 1.5 * ((1.0 - f_val**2) ** 1.5)

def sheet_thickness(r_tilde, theta):
    """Return sheet thickness in meters for a dimensionless radius r_tilde."""
    r_tilde = np.asarray(r_tilde)
    if np.any(r_tilde <= 0.0):
        raise ValueError("Sheet thickness is singular at the stagnation point (r = 0).")
    return (
        d_j
        * g_factor(f)
        / (4.0 * r_tilde * (1.0 + f * np.cos(theta))**2)
    )

# --- 4. RIGHT-HAND SIDE (RHS) OF THE ODE SYSTEM ---
def rim_derivatives(theta, y):
    u_b, s_b, phi, r_b = y
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    shape_factor = 1.0 + f * np.cos(theta)

    if u_b <= 0.0 or s_b <= 0.0 or r_b <= 0.0 or abs(sin_phi) < 1e-12:
        raise FloatingPointError(
            "The rim solution reached a singular or nonphysical state at "
            f"theta={theta:.6g}: u_b={u_b:.6g}, s_b={s_b:.6g}, "
            f"phi={phi:.6g}, r_b={r_b:.6g}."
        )

    # Eq. (3) expressed with r normalized by d_j, then F = h/r.
    F_val = g_factor(f) / (4.0 * r_b**2 * shape_factor**2)

    # Eqs. (6) and (7), solved for the derivatives of u_b and s_b.
    du_b_dtheta = F_val * (cos_phi - u_b) / (u_b * s_b)
    ds_b_dtheta = F_val * (2.0 * u_b - cos_phi) / (u_b**2)

    # Eq. (8) advances the rim tangent angle; Eq. (9) advances its radius.
    curvature_balance = 2.0 * r_b / We - sin_phi**2 * F_val
    dphi_dtheta = curvature_balance / (u_b**2 * s_b * sin_phi) - 1.0
    dr_b_dtheta = r_b * cos_phi / sin_phi

    return np.array([du_b_dtheta, ds_b_dtheta, dphi_dtheta, dr_b_dtheta])

# --- 5. CUSTOM RUNGE-KUTTA 4th ORDER (RK4) SOLVER ---
theta_start = 0.0
theta_end = np.pi
num_steps = 4000

theta_vals = np.linspace(theta_start, theta_end, num_steps)
dtheta = theta_vals[1] - theta_vals[0]

epsilon = 0.01
r_0 = (We * g_factor(f)) / (8.0 * (1.0 + f))

# The paper's printed u_b(0)=epsilon and phi(0) conditions are inconsistent;
# use its phi condition and enforce the defining relation u_b(0)=cos(phi(0)).
phi_init = np.arcsin(np.sqrt(1.0 - epsilon))
u_b_init = np.cos(phi_init)
r_init = r_0 * (1.0 - epsilon)

# Assume the initial rim area equals one jet's cross-sectional area,
# nondimensionalized by d_j**2 (d_j is the LOX-orifice diameter here).
s_b_init = np.pi / 4.0
if s_b_init <= 0.0:
    raise ValueError("The initial rim area is nonpositive for the selected parameters.")

y = np.array([u_b_init, s_b_init, phi_init, r_init])

trajectory = []
for th in theta_vals:
    trajectory.append(y.copy())
    
    # Standard RK4 algorithm steps
    k1 = rim_derivatives(th, y)
    k2 = rim_derivatives(th + 0.5 * dtheta, y + 0.5 * dtheta * k1)
    k3 = rim_derivatives(th + 0.5 * dtheta, y + 0.5 * dtheta * k2)
    k4 = rim_derivatives(th + dtheta, y + dtheta * k3)
    
    y += (dtheta / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

trajectory = np.array(trajectory)
r_vals = trajectory[:, 3]

# --- 6. PLOTTING THE THICKNESS FIELD AND BOUNDING RIM ---
theta_full = np.concatenate([-theta_vals[::-1], theta_vals])
r_full = np.concatenate([r_vals[::-1], r_vals])

x_coords = r_full * np.cos(theta_full)
y_coords = r_full * np.sin(theta_full)

d_j_mm = d_j * 1e3
x_coords_mm = x_coords * d_j_mm
y_coords_mm = y_coords * d_j_mm

plot_resolution = 601
plot_extent = r_vals.max()
plot_axis = np.linspace(-plot_extent, plot_extent, plot_resolution)
x_axis_mm = plot_axis * d_j_mm
x_grid_mm, y_grid_mm = np.meshgrid(x_axis_mm, x_axis_mm)
x_grid, y_grid = np.meshgrid(plot_axis, plot_axis)
r_grid = np.hypot(x_grid, y_grid)
theta_grid = np.arctan2(np.abs(y_grid), x_grid)
rim_radius_grid = np.interp(theta_grid, theta_vals, r_vals)
inside_sheet = (r_grid > 0.0) & (r_grid <= rim_radius_grid)

thickness_nm = np.full_like(r_grid, np.nan)
thickness_nm[inside_sheet] = sheet_thickness(
    r_grid[inside_sheet], theta_grid[inside_sheet]
) * 1e9
rim_thickness_nm = sheet_thickness(r_vals, theta_vals) * 1e9
rim_thickness_full_nm = np.concatenate(
    [rim_thickness_nm[::-1], rim_thickness_nm]
)
rim_angle_deg = np.linspace(0.0, 180.0, num_steps)
rim_model_theta = np.pi - np.abs(np.radians(rim_angle_deg) - np.pi / 2.0)
rim_radius_by_angle = np.interp(rim_model_theta, theta_vals, r_vals)
rim_thickness_by_angle_nm = (
    sheet_thickness(rim_radius_by_angle, rim_model_theta) * 1e9
)
minimum_thickness_nm = 1e3
thickness_above_cutoff_nm = np.ma.masked_less(
    np.ma.masked_invalid(thickness_nm), minimum_thickness_nm
)
visible_rim_x_mm = np.where(
    rim_thickness_full_nm >= minimum_thickness_nm, x_coords_mm, np.nan
)
visible_rim_y_mm = np.where(
    rim_thickness_full_nm >= minimum_thickness_nm, y_coords_mm, np.nan
)
visible_rim_thickness_nm = np.ma.masked_less(
    rim_thickness_nm, minimum_thickness_nm
)

fig, (ax_sheet, ax_rim) = plt.subplots(
    1, 2, figsize=(13, 6), gridspec_kw={"width_ratios": [1.25, 1]}
)
thickness_plot = ax_sheet.pcolormesh(
    x_grid_mm,
    y_grid_mm,
    thickness_above_cutoff_nm,
    shading="auto",
    cmap="viridis",
    norm=plt.matplotlib.colors.LogNorm(
        vmin=minimum_thickness_nm,
        vmax=float(thickness_above_cutoff_nm.max()),
    ),
)
ax_sheet.plot(
    visible_rim_x_mm,
    visible_rim_y_mm,
    color="white",
    lw=1.2,
    label="Rim (thickness ≥ 1000 nm)",
)
ax_sheet.set_title("Liquid-sheet thickness")
ax_sheet.set_xlabel("x (mm)")
ax_sheet.set_ylabel("y (mm)")
ax_sheet.set_aspect("equal")
ax_sheet.legend(loc="upper right")
fig.colorbar(thickness_plot, ax=ax_sheet, label="sheet thickness (nm, log scale)")

ax_rim.plot(
    rim_angle_deg, rim_thickness_by_angle_nm, color="tab:blue", lw=2
)
ax_rim.set_title("Rim thickness vs. angle :)))) ")
ax_rim.set_xlabel("angel")
ax_rim.set_xlim(0.0, 180.0)
ax_rim.set_ylabel("Sheet thikness (nm)")
ax_rim.grid(True, alpha=0.35)
fig.suptitle("Triple-jet liquid sheet")
fig.tight_layout()
plt.show()