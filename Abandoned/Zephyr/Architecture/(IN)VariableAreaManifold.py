import math
import matplotlib.pyplot as plt

# =============================
# USER CONFIGURATION
# =============================

CONFIG = {
    "num_outlets_total": 50,         # total outlets in full ring
    "A_outlet_cm2": 0.150,           # outlet area (cm²) - not used for plenum A but kept for reference
    "rho": 810.0,                    # fluid density (kg/m³)
    "mass_flow_kg_s": 3.0736,        # total inlet mass flow (kg/s)
    # plenum design velocity (m/s)
    "v_plenum_m_s": 3.0,

    # fraction of the half-ring mass flow that enters at the start (1.0 for single-end feed)
    "m_start_fraction": 1.0,

    # global scale on computed areas (1.0 => no scaling)
    "plenum_scale_factor": 1.0,

    "verbose": True
}

# =============================
# 📏 UNIT CONVERSIONS
# =============================

CM2_TO_M2 = 1e-4
def cm2_to_m2(a): return a * CM2_TO_M2
def m2_to_cm2(a): return a / CM2_TO_M2

# additional conversions for plotting in inches
CM2_TO_MM2 = 100.0
MM2_TO_IN2 = 0.00155

# helper: area -> equivalent circular diameter (mm)
def area_m2_to_diameter_mm(A_m2):
    if A_m2 <= 0:
        return 0.0
    d_m = math.sqrt(4.0 * A_m2 / math.pi)
    return d_m * 1000.0

# =============================
# 💨 HALF-RING MANIFOLD (constant velocity design)
# =============================

def ring_manifold_half(cfg):
    n_total = cfg["num_outlets_total"]
    n_half = n_total // 2

    m_total = cfg["mass_flow_kg_s"]
    m_half = m_total / 2.0                # mass flow allocated to one half-ring

    rho = cfg["rho"]
    v_plenum = cfg["v_plenum_m_s"]

    # mass that actually enters the start of this half-ring
    m_start = m_half * cfg["m_start_fraction"]

    # minimum flow at the last outlet (only one last outlet draws flow)
    m_per_outlet = m_start / n_half

    if cfg["verbose"]:
        print(f"Total m: {m_total:.6f} kg/s | Half: {m_half:.6f} kg/s")
        print(f"m_start_fraction: {cfg['m_start_fraction']:.3f} -> m_start = {m_start:.6f} kg/s")
        print(f"Plenum velocity (design): {v_plenum:.3f} m/s")
        print("-" * 60)

    # compute inlet area (largest section)
    A_start_m2 = m_start / (rho * v_plenum)
    A_start_cm2 = m2_to_cm2(A_start_m2)
    d_start_mm = area_m2_to_diameter_mm(A_start_m2)

    if cfg["verbose"]:
        print(f"A_start = {A_start_m2:.6e} m² = {A_start_cm2:.4f} cm² | equiv. dia = {d_start_mm:.1f} mm")
        print("-" * 60)

    # -----------------------------
    # LINEAR TAPER (constant-velocity manifold)
    # -----------------------------
    A_plenum_cm2 = []
    for i in range(n_half):

        # normalized index 0..1
        i_frac = i / (n_half - 1) if n_half > 1 else 0.0

        m_remain = m_start - i_frac * (m_start - m_per_outlet)

        # area from continuity: A = ṁ / (ρ v)
        A_i_m2 = m_remain / (rho * v_plenum)
        A_plenum_cm2.append(m2_to_cm2(A_i_m2))

        if cfg["verbose"]:
            d_mm = area_m2_to_diameter_mm(A_i_m2)
            print(f"Outlet {i+1:02d}: m_remain={m_remain:.6f} kg/s | "
                  f"A_raw={m2_to_cm2(A_i_m2):.4f} cm² | dia ~ {d_mm:.1f} mm")

    
    # apply global scale if desired
    scale = cfg["plenum_scale_factor"]
    A_plenum_cm2 = [a * scale for a in A_plenum_cm2]

    if cfg["verbose"]:
        DeadCSA =  400    # Dead CSA added to protect from friction
        print("-" * 60)
        print(f"Applied global scale factor: {scale:.3f}×")
        print(f"Largest area (scaled + DeadCSA): {(A_plenum_cm2[0]*100 + DeadCSA):.4f} mm² | {(A_plenum_cm2[0]*100 + DeadCSA)*MM2_TO_IN2:.4f} in²")
        mid_idx = len(A_plenum_cm2)//2
        print(f"Mid area (scaled + DeadCSA): {(A_plenum_cm2[mid_idx]*100 + DeadCSA):.4f} mm² | {(A_plenum_cm2[mid_idx]*100 + DeadCSA)*MM2_TO_IN2:.4f} in²")
        print(f"Smallest area (end of half-ring, scaled + DeadCSA): {(A_plenum_cm2[-1]*100 + DeadCSA):.4f} mm² | {(A_plenum_cm2[-1]*100 + DeadCSA)*MM2_TO_IN2:.4f} in²")
        print("-" * 60)

    return A_plenum_cm2

# =============================
# PLOTTING (in inches, +400 mm²)
# =============================

def plot_half_ring(A_plenum_cm2):
    n_half = len(A_plenum_cm2)

    # Convert plenum area (cm²) → mm² → add 400 → in²
    A_in2 = [ (a*100 + 400) * MM2_TO_IN2 for a in A_plenum_cm2 ]

    # mirrored for ring shape
    degrees_left = [-180 + (180 * i / n_half) for i in range(n_half)]
    areas_left = A_in2[::-1]

    degrees_right = [0 + (180 * i / n_half) for i in range(n_half)]
    areas_right = A_in2

    degrees = degrees_left + degrees_right
    areas = areas_left + areas_right

    plt.figure(figsize=(10,5))
    plt.plot(degrees, areas, marker='o', linewidth = 2)
    plt.xlabel("Degrees Around Manifold [-180° to 180°]", fontsize = 14)
    plt.ylabel("Plenum Cross-Sectional Area [in²]", fontsize = 14)
    plt.title("Ring Manifold Plenum Area vs Angular Distance from Input", fontsize = 18)
    plt.grid(True)
    plt.show()

# =============================
# ▶️ MAIN
# =============================

if __name__ == "__main__":
    A_plenum = ring_manifold_half(CONFIG)
    print("\nSimulation complete ✅\n")
    plot_half_ring(A_plenum)
