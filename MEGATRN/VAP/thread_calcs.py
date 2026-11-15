# Pintle Thread Calc Code
# Authors: Vasundhra Arulazi
# First Created: 11/1/2026
# Last Updated: 11/15/2026
"""
VAP pintle tip <-> shaft thread calcs. Formulas only -- every number comes from the YAML.

Copper tip = MALE thread, 316 SS pintle shaft = FEMALE thread.
Tip is the COMPRESSION member, shaft is the TENSION member.

  1. Separating force   F_p = (P_manifold - P_atm) * A_eff
  2. Stiffness split     k = E*A/L,  C = k_s / (k_s + k_p)   (1 - C goes to the tip)
  3. Preload             F_i - losses >= (1 - C) * F_p,  losses = thermal + embed + relax
                         upper bound = min(Cu stem, Cu strip, shoulder, SS sleeve, SS strip)
                         T = K d F_i  (full Shigley K, swept over friction)
  4. Threads             standard UNF + UNC; Cu wall grown (bore necked) so each one fits.
                         RoyMech min engagement, FED-STD-H28 strip check at F_max = F_i + C F_p
  5. Von Mises           Cu ring: axial + hoop + torsion (general plane stress)

Refs: Shigley 10th ed. Ch. 8, RoyMech thread engagement, FED-STD-H28, ASME B1.1, VDI 2230.
Run:  python MEGATRN/VAP/thread_calcs.py [inputs.yaml]     (default: thread_calcs.yaml next to this file)

Units: in, lbf, psi, deg F.
"""
import os
import sys
from math import pi, sqrt, tan, cos, atan, radians, ceil

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
YAML_PATH = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "thread_calcs.yaml")

THREAD_HALF_ANGLE = radians(30)     # constant: 60 deg UN thread form


def to_number(value):
    """YAML reads things like 17.0e6 as text -- turn number-looking text back into numbers."""
    if isinstance(value, list):
        return [to_number(x) for x in value]
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value
    return value


with open(YAML_PATH) as f:
    inputs = {key: to_number(value) for key, value in yaml.safe_load(f).items()}

# ---- derived inputs (null in the YAML -> default from the engine numbers) ----
tip_od = inputs["diam_ratio"] * inputs["chamber_id"]
room_temp = inputs["room_temp"]
manifold_pressure = inputs["manifold_pressure"] or (1 + inputs["lox_stiffness"]) * inputs["chamber_pressure"]
pressure_dia = inputs["effective_dia"] or tip_od
shaft_od = inputs["shaft_od"] or tip_od
lox_mass_flow = inputs["mass_flow"] * inputs["of_ratio"] / (1 + inputs["of_ratio"])


# ============================ YIELD vs TEMPERATURE ============================
cu_yield_fit = np.polyfit(inputs["cu_yield_temps"], inputs["cu_yield_ratio"], inputs["cu_yield_fit_order"])

#HELPER FUNCTIONS
def cu_yield_at(temp):
    """Cu yield at temp: RT yield x curve fit of empirical % of RT. No credit below room temp."""
    temp = min(max(temp, room_temp), max(inputs["cu_yield_temps"]))
    ratio = np.polyval(cu_yield_fit, temp) / np.polyval(cu_yield_fit, room_temp)
    return inputs["cu_yield"] * min(1.0, float(ratio))


def ss_yield_at(temp):
    """SS yield at temp: RT yield x linear interp of empirical % of RT. No credit below room temp."""
    temp = min(max(temp, room_temp), max(inputs["ss_yield_temps"]))
    return inputs["ss_yield"] * float(np.interp(temp, inputs["ss_yield_temps"], inputs["ss_yield_ratio"]))


# ============================ THREAD GEOMETRY ============================
def load_threads():
    """Read every thread table listed in the YAML -> [(major_dia, tpi, name), ...]."""
    threads = []
    for filename in inputs["thread_tables"]:
        series = os.path.splitext(filename)[0].split("_")[0].upper()          # unf_threads.csv -> UNF
        with open(os.path.join(HERE, filename)) as f:
            for line in f:
                if line.strip():
                    major_dia, tpi, name = line.strip().split(",", 2)
                    threads.append((float(major_dia), int(tpi), f"{name.strip()} {series}"))
    return sorted(threads)


def thread_geometry(major_dia, tpi):
    """UN thread dims + class 2A (external) / 2B (internal) worst-case limits, ASME B1.1."""
    pitch = 1.0 / tpi
    pitch_dia = major_dia - 0.649519 * pitch                 # basic pitch dia
    root_dia = major_dia - 1.299038 * pitch                  # external (Cu) root dia
    int_minor_dia = major_dia - 1.082532 * pitch             # basic internal (SS) minor dia

    # tolerances
    pitch_dia_tol = 0.0015 * major_dia**(1 / 3) + 0.0015 * sqrt(9 * pitch) + 0.015 * pitch**(2 / 3)
    allowance = 0.300 * pitch_dia_tol
    major_dia_tol = 0.060 * pitch**(2 / 3)
    int_minor_dia_tol = 0.25 * pitch - 0.40 * pitch**2

    # worst-case limits
    ext_major_min = major_dia - allowance - major_dia_tol
    ext_pitch_min = pitch_dia - allowance - pitch_dia_tol
    int_pitch_max = pitch_dia + 1.30 * pitch_dia_tol
    int_minor_max = int_minor_dia + int_minor_dia_tol

    # FED-STD-H28 shear area per inch of engagement
    cu_shear_area_per_in = pi * tpi * int_minor_max * (1 / (2 * tpi) + 0.57735 * (ext_pitch_min - int_minor_max))
    ss_shear_area_per_in = pi * tpi * ext_major_min * (1 / (2 * tpi) + 0.57735 * (ext_major_min - int_pitch_max))

    return dict(pitch=pitch, major_dia=major_dia, pitch_dia=pitch_dia, root_dia=root_dia,
                cu_shear_area_per_in=cu_shear_area_per_in, ss_shear_area_per_in=ss_shear_area_per_in,
                tensile_stress_area=pi / 4 * (major_dia - 0.9743 * pitch)**2)   # solid, standard


# ============================ 1. SEPARATING FORCE ============================
pressure_area = pi / 4 * pressure_dia**2
separating_force = (manifold_pressure - inputs["atm_pressure"]) * pressure_area


# ============================ 2-5 ============================
def shoulder_width(major_dia):
    """Shoulder bears on the SS end face, so it can't be wider than the SS wall over the threads."""
    return min(inputs["shoulder_width"], (shaft_od - major_dia) / 2)


def thread_torque_arm(friction, thread):
    """Thread torque per unit preload, T_thread = F_i * arm (Shigley power-screw term)."""
    lead_angle = atan(thread["pitch"] / (pi * thread["pitch_dia"]))
    flank_sec = 1 / cos(THREAD_HALF_ANGLE)
    return ((thread["pitch_dia"] / 2) * (tan(lead_angle) + friction * flank_sec)
            / (1 - friction * tan(lead_angle) * flank_sec))


def torque_coefficient(friction, thread):
    """Full Shigley K, T = K d F. Collar term = Cu shoulder friction radius, not a hex head."""
    shoulder_outer_radius = shaft_od / 2
    shoulder_inner_radius = shoulder_outer_radius - shoulder_width(thread["major_dia"])
    shoulder_friction_radius = ((2 / 3) * (shoulder_outer_radius**3 - shoulder_inner_radius**3)
                                / (shoulder_outer_radius**2 - shoulder_inner_radius**2))
    return (thread_torque_arm(friction, thread) + friction * shoulder_friction_radius) / thread["major_dia"]

#ACTUAL CALCULATION CODE IS BELOW
def analyze(major_dia, tpi, name, bore, cu_temp=None, ss_temp=None):
    #setting up
    cu_temp = inputs["cu_temp"] if cu_temp is None else cu_temp
    ss_temp = inputs["ss_temp"] if ss_temp is None else ss_temp
    thread = thread_geometry(major_dia, tpi)
    grip = inputs["grip_length"]
    cu_yield, ss_yield = cu_yield_at(cu_temp), ss_yield_at(ss_temp)
    cu_wall = (thread["root_dia"] - bore) / 2

    # ---- 2. stiffness split ----
    # this works out how the pressure load gets divided between Cu and SS
    cu_ring_area = pi / 4 * (thread["root_dia"]**2 - bore**2)     # thread root to bore (compression member)
    ss_sleeve_area = pi / 4 * (shaft_od**2 - major_dia**2)       # SS over the threads (tension member)
    cu_stiffness = inputs["cu_modulus"] * cu_ring_area / grip
    ss_stiffness = inputs["ss_modulus"] * ss_sleeve_area / grip
    shaft_load_share = ss_stiffness / (ss_stiffness + cu_stiffness)            # C

    # ---- 3. preload: F_i(1 - relax) - embed - thermal >= (1 - C) F_p ----
    #works out how much preload you need
    joint_compliance = 1 / ss_stiffness + 1 / cu_stiffness
    embed_loss = inputs["embedding"] / joint_compliance
    cold_strain_mismatch = max(0.0, inputs["cu_cold_contraction"] - inputs["ss_cold_contraction"])
    thermal_loss = cold_strain_mismatch * grip / joint_compliance if inputs["include_thermal"] else 0.0
    preload_required = (((1 - shaft_load_share) * separating_force + embed_loss + thermal_loss)
                        / (1 - inputs["relax_fraction"]))
    relax_loss = inputs["relax_fraction"] * preload_required

    # torque set at the highest friction (worst case still reaches preload_required);
    # lowest friction then over-tightens -> highest preload the copper sees
    friction_values = inputs["friction_sweep"]
    torque_coeff = {mu: torque_coefficient(mu, thread) for mu in friction_values}
    torque_spec = torque_coeff[max(friction_values)] * major_dia * preload_required
    preload_at_friction = {mu: torque_spec / (torque_coeff[mu] * major_dia) for mu in friction_values}
    preload_max = max(preload_at_friction.values())

    # ---- 4. engagement (RoyMech); strength ratio J > 1 -> SS (internal) is the weak side ----
    engagement_roymech = 2 * thread["tensile_stress_area"] / (0.5 * pi * thread["pitch_dia"])
    strength_ratio_J = ((thread["cu_shear_area_per_in"] * inputs["cu_ultimate"])
                        / (thread["ss_shear_area_per_in"] * inputs["ss_ultimate"]))
    if strength_ratio_J > 1:
        engagement_roymech *= strength_ratio_J
    threads_engaged = ceil(engagement_roymech * tpi)
    engagement_length = threads_engaged * thread["pitch"]

    # force limits (at temperature)
    shoulder_area = pi / 4 * (shaft_od**2 - (shaft_od - 2 * shoulder_width(major_dia))**2)
    force_limits = {
        "Cu stem": cu_yield * cu_ring_area,
        "Cu strip": 0.577 * cu_yield * thread["cu_shear_area_per_in"] * engagement_length,
        "Shoulder": cu_yield * shoulder_area,
        "SS sleeve": ss_yield * ss_sleeve_area,
        "SS strip": 0.577 * ss_yield * thread["ss_shear_area_per_in"] * engagement_length,
    }
    thread_force_max = preload_max + shaft_load_share * separating_force        # lowest friction
    thread_force_min = preload_required + shaft_load_share * separating_force   # highest friction
    fos_worst = {k: v / thread_force_max for k, v in force_limits.items()}
    fos_best = {k: v / thread_force_min for k, v in force_limits.items()}

    # ---- 5. von Mises on the Cu ring (11 axial, 22 hoop, 12 torsion), worst over friction ----
    cu_ring_polar_moment = pi / 32 * (thread["root_dia"]**4 - bore**4)
    hoop_stress_operating = manifold_pressure * bore / (2 * cu_wall)            # thin wall, LOx inside
    von_mises = {}
    for mu in friction_values:
        axial_stress = -preload_at_friction[mu] / cu_ring_area                  # compressive
        thread_torque = preload_at_friction[mu] * thread_torque_arm(mu, thread)
        torsion_shear = thread_torque * (thread["root_dia"] / 2) / cu_ring_polar_moment
        von_mises[mu] = {
            case: sqrt(axial_stress**2 - axial_stress * hoop_stress + hoop_stress**2 + 3 * torsion_shear**2)
            for case, hoop_stress in (("tighten", 0.0), ("operate", hoop_stress_operating))
        }
    for case in ("tighten", "operate"):
        fos_worst[f"Cu vM {case}"] = cu_yield / max(von_mises[mu][case] for mu in friction_values)
        fos_best[f"Cu vM {case}"] = cu_yield / von_mises[max(friction_values)][case]

    lox_velocity = lox_mass_flow / (inputs["lox_density"] * pi / 4 * bore**2) / 12    # ft/s

    return dict(name=name, major_dia=major_dia, tpi=tpi, bore=bore, cu_yield=cu_yield,
                cu_wall=cu_wall, ss_wall=(shaft_od - major_dia) / 2, lox_velocity=lox_velocity,
                cu_stiffness=cu_stiffness, ss_stiffness=ss_stiffness, shaft_load_share=shaft_load_share,
                embed_loss=embed_loss, thermal_loss=thermal_loss, relax_loss=relax_loss,
                preload_required=preload_required, torque_coeff=torque_coeff, torque_spec=torque_spec,
                preload_at_friction=preload_at_friction, engagement_roymech=engagement_roymech,
                strength_ratio_J=strength_ratio_J, threads_engaged=threads_engaged,
                engagement_length=engagement_length, force_limits=force_limits,
                min_force_limit=min(force_limits.values()), von_mises=von_mises, fos_worst=fos_worst,
                fos_low=min(fos_worst.values()), fos_high=min(fos_best.values()),
                governs=min(fos_worst, key=fos_worst.get))


def design(major_dia, tpi, name):
    """Grow the Cu wall (neck the bore) from min_wall until fos_target is met or the bore hits bore_min."""
    if (shaft_od - major_dia) / 2 < inputs["min_wall"]:
        return None
    root_dia = major_dia - 1.299038 / tpi
    cu_wall, best = inputs["min_wall"], None
    while root_dia - 2 * cu_wall >= inputs["bore_min"] - 1e-9:
        result = analyze(major_dia, tpi, name, root_dia - 2 * cu_wall)
        if best is None or result["fos_low"] > best["fos_low"]:
            best = result
        if result["fos_low"] >= inputs["fos_target"]:
            return result
        cu_wall += inputs["wall_step"]
    return best


# ============================ OUTPUT ============================
if __name__ == "__main__":
    threads = load_threads()
    results = [r for r in (design(*t) for t in threads) if r]
    if not results:
        sys.exit(f"No thread fits: check shaft_od / bore_min / min_wall in {YAML_PATH}")
    results.sort(key=lambda r: (r["fos_low"] >= inputs["fos_target"], r["fos_low"]), reverse=True)

    atm = inputs["atm_pressure"]
    print(f"Inputs: {YAML_PATH}")
    print("1. Separating force")
    print(f"   dP = {manifold_pressure:.1f} - {atm} = {manifold_pressure - atm:.1f} psi,"
          f"  A_eff = {pressure_area:.3f} in^2  ->  F_p = {separating_force:.0f} lbf")
    print(f"   Temps: Cu {inputs['cu_temp']} F (Sy {cu_yield_at(inputs['cu_temp'])/1e3:.1f} ksi),"
          f" SS {inputs['ss_temp']} F (Sy {ss_yield_at(inputs['ss_temp'])/1e3:.1f} ksi)")
    print(f"\n{len(threads)} threads checked, {len(results)} fit (OD {shaft_od:.3f}, bore >= {inputs['bore_min']},"
          f" walls >= {inputs['min_wall']}). FOS range = lowest .. highest friction\n")

    print(f"{'Thread':<18}{'Cu wall':>8}{'bore':>7}{'LOx v':>7}{'SS wall':>8}{'F_req':>7}{'T_spec':>8}"
          f"{'threads':>8}{'L_e':>7}{'governs':>15}{'FOS range':>13}  ok")
    for r in results:
        ok = "Y" if r["fos_low"] >= inputs["fos_target"] else "x"
        print(f"{r['name']:<18}{r['cu_wall']:>8.3f}{r['bore']:>7.3f}{r['lox_velocity']:>6.1f}'{r['ss_wall']:>8.3f}"
              f"{r['preload_required']:>7.0f}{r['torque_spec']:>8.0f}{r['threads_engaged']:>8d}"
              f"{r['engagement_length']:>7.3f}{r['governs']:>15}{r['fos_low']:>7.2f}-{r['fos_high']:<5.2f}  {ok}")

    best = results[0]
    C = best["shaft_load_share"]
    print(f"\n=== Best: {best['name']}  (Cu wall {best['cu_wall']:.3f} in, bore {best['bore']:.3f} in) ===")
    print(f"2. Stiffness   k_s = {best['ss_stiffness']/1e6:.2f}e6, k_p = {best['cu_stiffness']/1e6:.2f}e6 lbf/in"
          f"  ->  C = {C:.2f}, tip takes {1 - C:.0%}")
    print(f"3. Losses      embed {best['embed_loss']:.0f} + thermal {best['thermal_loss']:.0f}"
          f" + relax {best['relax_loss']:.0f} lbf  ->  F_req = {best['preload_required']:.0f} lbf")
    print(f"   Upper bound {best['min_force_limit']:.0f} lbf  ["
          + ", ".join(f"{k} {v:.0f}" for k, v in best["force_limits"].items()) + "]")
    print(f"   Torque spec {best['torque_spec']:.0f} in-lbf:")
    for mu in inputs["friction_sweep"]:
        preload = best["preload_at_friction"][mu]
        print(f"     mu {mu:.2f}:  K = {best['torque_coeff'][mu]:.3f}   F_i = {preload:5.0f}"
              f"   F_max = {preload + C * separating_force:5.0f} lbf")
    print(f"4. Engagement  RoyMech {best['engagement_roymech']:.3f} in (J = {best['strength_ratio_J']:.2f})"
          f"  ->  {best['threads_engaged']} threads = {best['engagement_length']:.3f} in")
    print("   FOS at max preload: " + "   ".join(f"{k} {v:.2f}" for k, v in best["fos_worst"].items()))

    print(f"\nCu temperature sensitivity ({best['name']}):")
    for temp in inputs["cu_temp_sweep"]:
        r = analyze(best["major_dia"], best["tpi"], best["name"], best["bore"], cu_temp=temp)
        print(f"   {temp:>5} F:  Sy {r['cu_yield']/1e3:5.1f} ksi   FOS {r['fos_low']:.2f}-{r['fos_high']:.2f}"
              f"   governs {r['governs']}")