# Pintle Fill Time Analysis Code
# Authors: Chanak Gautam
# First Created: 10/29/2026
# Last Updated: 10/29/2026
# Calculations done in SI unit

import math
import os
import matplotlib.pyplot as plt
import pint
from pyfluids import Fluid, FluidsList, Input
from CoolProp.CoolProp import PropsSI

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


# ---------------------------------------------------------------------------
# Inputs (imperial). All pressures are absolute.
# ---------------------------------------------------------------------------
LIQUID = FluidsList.nDodecane               # propellant (n-dodecane is the usual RP-1 surrogate)
GAS = FluidsList.Nitrogen                   # purge gas

V0_imp = Q_(2.0, "in**3")                   # manifold internal volume
P_f_imp = Q_(500.0, "psi")                  # feed pressure upstream of liquid inlet orifice
P_a_imp = Q_(14.7, "psi")                   # pressure downstream of gas outlet
P_amb_imp = Q_(14.7, "psi")                 # initial manifold pressure (purge gas at rest)
T_gas_imp = Q_(70.0, "degF")                # purge gas temperature (held constant, isothermal)
T_liq_imp = Q_(70.0, "degF")                # propellant temperature
D_in_imp = Q_(0.0625, "in")                 # liquid inlet orifice diameter
D_out_imp = Q_(0.0625, "in")                # gas outlet orifice diameter

Cd_liq = 0.61                               # estimate: sharp-edged orifice, liquid
Cd_gas = 0.80                               # estimate: gas through valve/fitting
dt = 1e-5                                   # s
t_max = 10.0                                # s, abort if the manifold hasn't filled by now


# ---------------------------------------------------------------------------
# Convert inputs to SI
# ---------------------------------------------------------------------------
V0 = V0_imp.m_as("m**3")                    # m^3
P_f = P_f_imp.m_as("Pa")                    # Pa
P_a = P_a_imp.m_as("Pa")                    # Pa
P_amb = P_amb_imp.m_as("Pa")                # Pa
T = T_gas_imp.m_as("K")                     # K
T_liq = T_liq_imp.m_as("K")                 # K
A_in = math.pi / 4 * D_in_imp.m_as("m") ** 2    # m^2
A_out = math.pi / 4 * D_out_imp.m_as("m") ** 2  # m^2

if P_f <= P_amb:
    raise ValueError("Feed pressure must exceed initial manifold pressure for liquid to flow in")


# ---------------------------------------------------------------------------
# Fluid properties (pyfluids takes temperature in degC, pressure in Pa)
# ---------------------------------------------------------------------------
liquid = Fluid(LIQUID).with_state(Input.pressure(P_f), Input.temperature(T_liq - 273.15))
gas = Fluid(GAS).with_state(Input.pressure(P_amb), Input.temperature(T - 273.15))

rho_L = liquid.density                      # kg/m^3, liquid treated as incompressible

R_u = PropsSI("gas_constant", GAS.coolprop_name)    # J/(mol K), universal gas constant
M_g = PropsSI("molar_mass", GAS.coolprop_name)      # kg/mol
R = R_u / M_g                               # J/(kg K), specific gas constant of purge gas
gamma = gas.specific_heat / PropsSI("Cvmass", "T", T, "P", P_amb, GAS.coolprop_name)

# Choked when P_a / P_m is at or below this critical pressure ratio
P_ratio_crit = math.pow(2 / (gamma + 1), gamma / (gamma - 1))


# ---------------------------------------------------------------------------
# Mass flow helpers
# Flows are clamped to zero instead of reversing when the pressure difference flips.
# ---------------------------------------------------------------------------
def liquid_mdot(P_m):
    """Liquid mass flow in through the inlet orifice [kg/s], incompressible orifice equation."""
    if P_m >= P_f:
        return 0.0
    return Cd_liq * A_in * math.sqrt(2 * rho_L * (P_f - P_m))


def gas_mdot(P_m):
    """Purge gas mass flow out through the outlet orifice [kg/s], isentropic ideal-gas nozzle.
    Returns (mdot, choked)."""
    if P_m <= P_a:
        return 0.0, False

    P_ratio = P_a / P_m
    if P_ratio <= P_ratio_crit:
        mdot = Cd_gas * A_out * P_m * math.sqrt(gamma / (R * T)) \
            * math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)))
        return mdot, True

    mdot = Cd_gas * A_out * P_m * math.sqrt(
        2 * gamma / ((gamma - 1) * R * T)
        * (math.pow(P_ratio, 2 / gamma) - math.pow(P_ratio, (gamma + 1) / gamma)))
    return mdot, False


# ---------------------------------------------------------------------------
# Initialize (t = 0, manifold full of purge gas)
# rho_g comes from the real-gas EOS rather than P/RT, so Z carries the real-gas correction
# (the ideal-gas rho_g in the flowchart would make Z = 1 by construction).
# ---------------------------------------------------------------------------
V_g = V0                                    # m^3
V_L = 0.0                                   # m^3
P_m = P_amb                                 # Pa
t_tot = 0.0                                 # s

rho_g = gas.density                         # kg/m^3
rho_g_0 = rho_g                             # kg/m^3, initial gas density
m_g = rho_g * V0                            # kg
Z = P_m * V_g / (m_g * R * T)               # compressibility factor


# Print out initial state
print(f"Liquid: {LIQUID.name}, rho_L = {Q_(rho_L, 'kg/m**3').m_as('lb/ft**3'):.3f} lb/ft^3")
print(f"Gas:    {GAS.name}, R = {R:.3f} J/(kg K), gamma = {gamma:.4f}, Z = {Z:.5f}")
print(f"        rho_g = {Q_(rho_g, 'kg/m**3').m_as('lb/ft**3'):.5f} lb/ft^3, "
      f"m_g = {Q_(m_g, 'kg').m_as('lb'):.4e} lb")
print(f"A_in = {Q_(A_in, 'm**2').m_as('in**2'):.4e} in^2, A_out = {Q_(A_out, 'm**2').m_as('in**2'):.4e} in^2")


# ---------------------------------------------------------------------------
# Time-marching loop
# Gas pressure is taken straight from the real-gas EOS at (rho_g, T) each step;
# Z is still computed for reporting.
# ---------------------------------------------------------------------------
history = {"t": [], "P_m": [], "V_L": [], "mdot_L": [], "mdot_g": [], "choked": [], "Z": []}

while True:
    mdot_L = liquid_mdot(P_m)
    mdot_g, choked = gas_mdot(P_m)

    history["t"].append(t_tot)
    history["P_m"].append(P_m)
    history["V_L"].append(V_L)
    history["mdot_L"].append(mdot_L)
    history["mdot_g"].append(mdot_g)
    history["choked"].append(choked)
    history["Z"].append(Z)

    # Update control volumes over dt (gas can't leave faster than what's left)
    V_g_prev = V_g
    V_L += mdot_L / rho_L * dt
    m_g -= min(mdot_g * dt, m_g)
    V_g = V0 - V_L

    # Manifold full: interpolate within the last step for when V_g hit zero
    if V_g <= 0:
        t_fill = t_tot + dt * V_g_prev / (V_g_prev - V_g)
        break

    t_tot += dt
    if t_tot > t_max:
        raise RuntimeError(f"Manifold not full after {t_max} s (V_g = {V_g:.3e} m^3)")

    # Update gas state
    rho_g = m_g / V_g
    P_m = PropsSI("P", "Dmass", rho_g, "T", T, GAS.coolprop_name)
    Z = P_m / (rho_g * R * T)


# Print out results
P_m_hist = history["P_m"]
i_peak = max(range(len(P_m_hist)), key=P_m_hist.__getitem__)
n_choked = sum(history["choked"])
print(f"\nFill time = {t_fill * 1e3:.3f} ms ({len(P_m_hist)} steps, dt = {dt:.0e} s)")
print(f"Peak P_m  = {Q_(P_m_hist[i_peak], 'Pa').m_as('psi'):.2f} psia at t = {history['t'][i_peak] * 1e3:.3f} ms")
print(f"Gas outlet choked for {n_choked * dt * 1e3:.3f} ms of the fill")
print(f"Purge gas vented: {(1 - m_g / (rho_g_0 * V0)) * 100:.2f} % of initial mass")


# ---------------------------------------------------------------------------
# Plots: manifold pressure and liquid inlet velocity vs time
# Liquid velocity is the mean velocity through the inlet orifice area, mdot_L / (rho_L * A_in).
# ---------------------------------------------------------------------------
t_ms = [t * 1e3 for t in history["t"]]
P_m_psia = [Q_(P, "Pa").m_as("psi") for P in P_m_hist]
v_L_fps = [Q_(mdot / (rho_L * A_in), "m/s").m_as("ft/s") for mdot in history["mdot_L"]]

line_color = "#2a78d6"
ref_color = "#52514e"

fig, (ax_P, ax_v) = plt.subplots(2, 1, sharex=True, figsize=(8, 6), constrained_layout=True)

ax_P.plot(t_ms, P_m_psia, color=line_color, linewidth=2)
ax_P.axhline(P_a_imp.m_as("psi"), color=ref_color, linewidth=1, linestyle="--")
ax_P.annotate("P_a", xy=(t_fill * 1e3, P_a_imp.m_as("psi")), xytext=(-8, 4),
              textcoords="offset points", ha="right", color=ref_color, fontsize=9)
ax_P.set_title("Manifold pressure", loc="left")
ax_P.set_ylabel("P_m [psia]")

ax_v.plot(t_ms, v_L_fps, color=line_color, linewidth=2)
ax_v.set_title("Liquid velocity through inlet orifice", loc="left")
ax_v.set_ylabel("v_L [ft/s]")
ax_v.set_xlabel("Time [ms]")

for ax in (ax_P, ax_v):
    ax.axvline(t_fill * 1e3, color=ref_color, linewidth=1, linestyle=":")
    ax.grid(True, color="#e5e4e0", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
ax_v.annotate(f"t_fill = {t_fill * 1e3:.1f} ms", xy=(t_fill * 1e3, 1), xycoords=("data", "axes fraction"),
              xytext=(-6, -4), textcoords="offset points", ha="right", va="top", color=ref_color, fontsize=9)

fig.suptitle(f"Manifold fill: {LIQUID.name} into {GAS.name}-purged manifold", x=0.01, ha="left")
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fill_time.png"), dpi=150)
plt.show()
