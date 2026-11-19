from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class LiquidStream:
    name: str
    mass_flow_rate_lbm_s: float
    density_kg_m3: float
    surface_tension_n_m: float


@dataclass(frozen=True)
class LiquidDefaults:
    name: str
    density_kg_m3: float
    surface_tension_n_m: float


def prompt_float(label, default: float | None, minimum):
    while True:
        prompt = f"{label} [{default:g}]: " if default is not None else f"{label}: "
        raw_value = input(prompt).strip()
        if not raw_value and default is None:
            print("A value is required.")
            continue
        try:
            value = default if not raw_value else float(raw_value)
        except ValueError:
            print("Enter a number.")
            continue
        if not np.isfinite(value) or value <= minimum:
            print(f"Enter a finite number greater than {minimum:g}.")
            continue
        return value


def prompt_liquid(role, defaults: LiquidDefaults):
    name = input(f"{role} liquid name [{defaults.name}]: ").strip() or defaults.name
    
    return LiquidStream(
        name=name,
        mass_flow_rate_lbm_s=prompt_float(
            "Total mass flow rate (lbm/s)", None, 0.0
        ),
        density_kg_m3=prompt_float(
            "Density (kg/m^3)", defaults.density_kg_m3, 0.0
        ),
        surface_tension_n_m=prompt_float(
            "Surface tension (N/m)", defaults.surface_tension_n_m, 0.0
        ),
    )


# --- 1. CONVERSION FACTORS & INPUT PARAMETERS ---
inch_to_meter = 0.0254
lbm_to_kg = 0.45359237

central_defaults = LiquidDefaults(
    name="Central liquid",
    density_kg_m3=1141.0,
    surface_tension_n_m=0.013,
)
outer_defaults = LiquidDefaults(
    name="Outer liquid",
    density_kg_m3=810.0,
    surface_tension_n_m=0.024,
)

print("Three-jet liquid-sheet fan model")
central_liquid = prompt_liquid("central", central_defaults)
outer_liquid = prompt_liquid("outer", outer_defaults)
outer_jet_count = 2
central_orifice_diameter = (
    prompt_float("Central orifice diameter (in)", 0.046, 0.0) * inch_to_meter
)
outer_orifice_diameters = np.array(
    [
        prompt_float(f"Outer orifice {index + 1} diameter (in)", 0.033, 0.0)
        * inch_to_meter
        for index in range(outer_jet_count)
    ]
)
collision_angle_degrees = prompt_float("Outer-jet collision angle (degrees)", 30.0, 0.0)

if collision_angle_degrees >= 90.0:
    raise ValueError("The collision angle must be less than 90 degrees.")

central_mass_flow_rate = central_liquid.mass_flow_rate_lbm_s * lbm_to_kg
outer_total_mass_flow_rate = outer_liquid.mass_flow_rate_lbm_s * lbm_to_kg
outer_mass_flow_rate_per_jet = outer_total_mass_flow_rate / outer_jet_count

central_density = central_liquid.density_kg_m3
central_surface_tension = central_liquid.surface_tension_n_m

outer_density = outer_liquid.density_kg_m3
outer_surface_tension = outer_liquid.surface_tension_n_m

half_collision_angle = np.radians(collision_angle_degrees)

central_orifice_area = np.pi * central_orifice_diameter**2 / 4.0
outer_orifice_areas = np.pi * outer_orifice_diameters**2 / 4.0

central_jet_velocity = central_mass_flow_rate / (
    central_density * central_orifice_area
)
outer_jet_velocities = outer_mass_flow_rate_per_jet / (
    outer_density * outer_orifice_areas
)

# --- 2. HELPER FUNCTIONS ---
def momentum_rate_jets(orifice_diameters, density, jet_velocities, alpha):
    jet_momenta = (
        np.pi
        / 4.0
        * np.asarray(orifice_diameters) ** 2
        * density
        * np.asarray(jet_velocities) ** 2
    )
    return np.sum(jet_momenta) * np.cos(alpha)


def mass_flow_jets(orifice_diameters, density, jet_velocities):
    return (
        np.pi
        / 4.0
        * np.asarray(orifice_diameters) ** 2
        * density
        * np.asarray(jet_velocities)
    )


def mixture_density_calculation(
    outer_density, central_density, d_outer, d_central, v_outer, v_central,
):
    q_outer = np.sum(np.pi * (d_outer**2) / 4.0 * v_outer)
    q_central = np.pi * (d_central**2) / 4.0 * v_central
    return (
        (outer_density * q_outer + central_density * q_central)
        / (q_outer + q_central)
    )


def mixture_sheet_velocity_calculation(
    outer_density, central_density, d_outer, d_central, v_outer, v_central,
):
    m_outer = mass_flow_jets(d_outer, outer_density, v_outer)
    m_central = mass_flow_jets(d_central, central_density, v_central)
    return (
        (np.sum(m_outer * v_outer) + np.sum(m_central * v_central))
        / (np.sum(m_outer) + np.sum(m_central))
    )


def mixture_surface_tension_calculation(
    sigma_outer, sigma_central, d_outer, d_central, v_outer, v_central,
):
    q_outer = np.sum(np.pi * (d_outer**2) / 4.0 * v_outer)
    q_central = np.pi * (d_central**2) / 4.0 * v_central
    return (
        (sigma_outer * q_outer + sigma_central * q_central)
        / (q_outer + q_central)
    )


# --- 3. MIXTURE PROPERTIES & SHEET COEFFICIENTS ---
rho_s = mixture_density_calculation(
    outer_density,
    central_density,
    outer_orifice_diameters,
    central_orifice_diameter,
    outer_jet_velocities,
    central_jet_velocity,
)
u_s = mixture_sheet_velocity_calculation(
    outer_density,
    central_density,
    outer_orifice_diameters,
    central_orifice_diameter,
    outer_jet_velocities,
    central_jet_velocity,
)
sigma_s = mixture_surface_tension_calculation(
    outer_surface_tension,
    central_surface_tension,
    outer_orifice_diameters,
    central_orifice_diameter,
    outer_jet_velocities,
    central_jet_velocity,
)

m_outer_jets = mass_flow_jets(
    outer_orifice_diameters, outer_density, outer_jet_velocities
)
m_central_single = mass_flow_jets(
    central_orifice_diameter, central_density, central_jet_velocity
)
LHS_Eq1 = np.sum(m_outer_jets) + m_central_single

p_outer_axial = momentum_rate_jets(
    outer_orifice_diameters,
    outer_density,
    outer_jet_velocities,
    half_collision_angle,
)
p_central_axial = momentum_rate_jets(
    central_orifice_diameter, central_density, central_jet_velocity, 0.0
)
LHS_Eq2 = p_outer_axial + p_central_axial

f = LHS_Eq2 / (LHS_Eq1 * u_s)
if not 0.0 <= f < 1.0:
    raise ValueError(f"The momentum ratio f must be in [0, 1); got {f:.6g}.")

d_j = central_orifice_diameter
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
# nondimensionalized by d_j**2 (d_j is the central-orifice diameter here).
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
minimum_radius_for_cutoff = (
    d_j
    * g_factor(f)
    * 1e9
    / (
        4.0
        * minimum_thickness_nm
        * (1.0 + f * np.cos(theta_vals))**2
    )
)
visible_sheet_radius = np.minimum(r_vals, minimum_radius_for_cutoff)
cad_boundary_x_mm = visible_sheet_radius * np.cos(theta_vals) * d_j_mm
cad_boundary_y_mm = visible_sheet_radius * np.sin(theta_vals) * d_j_mm
cad_x_min_mm = float(np.min(cad_boundary_x_mm))
cad_x_max_mm = float(np.max(cad_boundary_x_mm))
cad_station_x_mm = np.linspace(cad_x_min_mm, cad_x_max_mm, 11)
cad_station_rows = []
for station_x_mm in cad_station_x_mm:
    crossing_y_mm = []
    for x0, x1, y0, y1 in zip(
        cad_boundary_x_mm[:-1],
        cad_boundary_x_mm[1:],
        cad_boundary_y_mm[:-1],
        cad_boundary_y_mm[1:],
    ):
        if (x0 - station_x_mm) * (x1 - station_x_mm) <= 0.0:
            if x0 == x1:
                crossing_y_mm.extend((y0, y1))
            else:
                fraction = (station_x_mm - x0) / (x1 - x0)
                crossing_y_mm.append(y0 + fraction * (y1 - y0))
    if not crossing_y_mm:
        raise RuntimeError(
            f"No cutoff-boundary point found at x={station_x_mm:.6f} mm."
        )
    half_width_mm = float(np.max(crossing_y_mm))
    cad_station_rows.append((station_x_mm, half_width_mm))

maximum_width_index = int(np.argmax(cad_boundary_y_mm))
maximum_width_mm = 2.0 * float(cad_boundary_y_mm[maximum_width_index])
maximum_width_x_mm = float(cad_boundary_x_mm[maximum_width_index])

max_x_distance_in = (
    np.max(visible_sheet_radius * np.abs(np.cos(theta_vals))) * d_j_mm
    / 25.4
)
print(
    "Calculated cutoff-outline dimensions "
    f"(threshold={minimum_thickness_nm:.0f} nm): "
    f"x={cad_x_min_mm:.6f} to {cad_x_max_mm:.6f} mm, "
    f"maximum width={maximum_width_mm:.6f} mm "
    f"at x={maximum_width_x_mm:.6f} mm."
)
print("CAD profile stations (all dimensions in mm; symmetric about y=0):")
print("distance_from_min_x, x, lower_y, upper_y, full_width")
for station_x_mm, half_width_mm in cad_station_rows:
    print(
        f"{station_x_mm - cad_x_min_mm:.6f}, {station_x_mm:.6f}, "
        f"{-half_width_mm:.6f}, {half_width_mm:.6f}, "
        f"{2.0 * half_width_mm:.6f}"
    )
print(
    "Maximum x-distance of the sheet from (0, 0) "
    f"(thickness >= {minimum_thickness_nm:.0f} nm): "
    f"{max_x_distance_in:.3f} in"
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
fig.suptitle(
    f"{outer_liquid.name} / {central_liquid.name} liquid sheet "
    f"({outer_jet_count} outer jets)"
)
fig.tight_layout()
plt.show()
