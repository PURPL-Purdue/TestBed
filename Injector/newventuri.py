# ════════════════════════════════════════════════════════════════
#   Maelstrom Cavitating Venturi Sizing + Inlet Pressure Mapping
#   Author: Dominik Sloup
#   Updated: 06/01/2025
#   SI Units Throughout
# ════════════════════════════════════════════════════════════════

import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input

# ------------------------------------------------------------
# UNIT CONSTANTS
# ------------------------------------------------------------
psi_to_pa = 6894.76
pa_to_psi = 1.0 / psi_to_pa
in_to_m = 0.0254

# ------------------------------------------------------------
# FLUID + GEOMETRY PARAMETERS
# ------------------------------------------------------------
rho_rp1 = 810.0     # RP-1 density [kg/m^3]
C_d     = 0.931     # Venturi discharge coefficient
K_rec   = 0.67      # Pressure recovery factor

# Desired mass flow from feed system code (replace later)
m_dot_target = 0.533   # kg/s

# Line diameter for venturi sizing (you used -8 tube)
d_line_in = 0.426
A_line_m2 = np.pi * (d_line_in * in_to_m)**2 / 4

# Surrogate saturation pressure (need actual RP-1 vapor pressure)
p_sat_pa = Fluid(FluidsList.nDodecane).with_state(
    Input.quality(0.0),
    Input.temperature(20)
).pressure
p_sat_psi = p_sat_pa * pa_to_psi
print(p_sat_pa)
print(Fluid(FluidsList.Water).with_state(
    Input.quality(0.0),
    Input.temperature(20)
).pressure)

# ------------------------------------------------------------
#   1. DIRECT CAVITATING VENTURI SIZING
# ------------------------------------------------------------
def compute_venturi_throat_area(m_dot, Cd, rho, p0_pa, p_sat_pa):
    """
    Computes throat area assuming cavitating venturi behavior.
    Uses the modified Bernoulli equation.
    """
    A_star = 1 / np.sqrt(
        (1 / A_line_m2**2)
        - (2 * rho * (p_sat_pa - p0_pa)) / m_dot**2
    )
    return A_star

# Example design point
p1_psi = 750
x_recovery = 0.9               # p0 / p1 ratio
p0_psi = p1_psi / x_recovery
p0_pa  = p0_psi * psi_to_pa

A_star_m2 = compute_venturi_throat_area(
    m_dot_target, C_d, rho_rp1, p0_pa, p_sat_pa
)

d_star_m = 2 * np.sqrt(A_star_m2 / np.pi)
d_star_in = d_star_m / in_to_m

print("\n──────── SIZING RESULTS ────────")
print(f"Desired m_dot             = {m_dot_target:.3f} kg/s")
print(f"Venturi inlet pressure    = {p0_psi:.1f} psi")
print(f"Venturi throat diameter   = {d_star_in:.3f} in")
print(f"Venturi throat area       = {A_star_m2:.6e} m²\n")

# ------------------------------------------------------------
#   2. PRESSURE MAPPING FUNCTION
# ------------------------------------------------------------
def analyze_venturi_from_inlet(
    p0_pa,
    m_dot,
    Cd,
    A_star_m2,
    rho,
    K_rec,
    p_sat_pa
):
    """
    Given inlet pressure p0 and required mass flow m_dot:
      • compute throat pressure p_t
      • compute recovered manifold pressure p1
      • detect cavitation
    """
    if m_dot <= 0:
        return p0_pa, p0_pa, False

    CdA = Cd * A_star_m2

    # standard orifice equation
    Δp = (m_dot / CdA)**2 / (2 * rho)

    p_t = p0_pa - Δp
    p1 = p_t + K_rec * Δp
    cavitating = p_t < p_sat_pa

    return p_t, p1, cavitating

# ------------------------------------------------------------
#   3. INLET PRESSURE SWEEP
# ------------------------------------------------------------
p0_psi_array = np.linspace(10, 900, 200)
p0_pa_array  = p0_psi_array * psi_to_pa

# placeholder engine mass flow – REPLACE with your actual engine solver
m_dot_engine_array = 0.0009 * p0_psi_array

p_t_array = np.zeros_like(p0_pa_array)
p1_array  = np.zeros_like(p0_pa_array)
cav_array = np.zeros_like(p0_pa_array, dtype=bool)

for i in range(len(p0_pa_array)):
    p_t, p1, cav = analyze_venturi_from_inlet(
        p0_pa_array[i],
        m_dot_engine_array[i],
        C_d,
        A_star_m2,
        rho_rp1,
        K_rec,
        p_sat_pa
    )
    p_t_array[i] = p_t
    p1_array[i]  = p1
    cav_array[i] = cav

p_t_psi = p_t_array * pa_to_psi
p1_psi = p1_array * pa_to_psi

# ------------------------------------------------------------
#   4. PLOT
# ------------------------------------------------------------
plt.figure(figsize=(10, 6))

plt.plot(m_dot_engine_array, p0_psi_array, label="Inlet Pressure $p_0$", linestyle="--")
plt.plot(m_dot_engine_array, p_t_psi,        label="Throat Pressure $p_t$")
plt.plot(m_dot_engine_array, p1_psi,        label="Recovered Pressure $p_1$")

plt.axhline(p_sat_psi, color='r', linestyle=':', label="Saturation Pressure (proxy)")

plt.xlabel("Mass Flow Rate [kg/s]")
plt.ylabel("Pressure [psi]")
plt.title("Venturi Pressure Mapping vs Engine Mass Flow")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ------------------------------------------------------------
# Print sample rows
# ------------------------------------------------------------
print("──────── PRESSURE MAP SAMPLE ────────")
for i in [0, 50, 100, 150]:
    print(
        f"p0={p0_psi_array[i]:6.1f} psi | "
        f"pt={p_t_psi[i]:7.2f} psi | "
        f"p1={p1_psi[i]:7.2f} psi | "
        f"cav={cav_array[i]}"
    )