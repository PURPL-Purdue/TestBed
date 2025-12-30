import math
import matplotlib.pyplot as plt

# =============================
# 🧰 USER CONFIGURATION
# =============================

CONFIG = {
    "num_outlets_total": 40,       # total number of outlet holes (full ring)
    "P_res_psi": 1000.0,           # fixed pump pressure, PSI
    "P_chamber_psi": 600.0,        # downstream chamber pressure, PSI
    "A_outlet_cm2": .150,          # outlet area (cm²)
    "rho": 810.0,                  # fluid density (kg/m³)
    "mass_flow_kg_s": 4.204,       # total inlet mass flow (kg/s)
    "verbose": True
}

# =============================
# 📏 UNIT CONVERSIONS
# =============================

PSI_TO_PA = 6894.75729
CM2_TO_M2 = 1e-4

def psi_to_pa(psi): return psi * PSI_TO_PA
def m2_to_cm2(a): return a / CM2_TO_M2
def cm2_to_m2(a): return a * CM2_TO_M2

# =============================
# 💨 HALF-RING CONSTANT VELOCITY MANIFOLD
# =============================

def ring_manifold_half(cfg):
    n_total = cfg["num_outlets_total"]
    n_half = n_total // 2
    m_total = cfg["mass_flow_kg_s"]
    m_half = m_total / 2
    rho = cfg["rho"]
    A_outlet_m2 = cm2_to_m2(cfg["A_outlet_cm2"])

    # outlet velocity (constant bleed)
    v_target = m_half / (rho * n_half * A_outlet_m2)

    A_plenum = []
    outlet_velocities = []
    outlet_massflows = []

    for i in range(n_half):
        m_remain = m_half - i * (m_half / n_half)
        A_i_m2 = m_remain / (rho * v_target)
        A_plenum.append(m2_to_cm2(A_i_m2))

        outlet_velocities.append(v_target)
        outlet_massflows.append(m_half / n_half)

        if cfg["verbose"]:
            print(f"Outlet {i+1:02d}: "
                  f"{v_target:.3f} m/s | "
                  f"{m_half/n_half:.3f} kg/s | "
                  f"Plenum {m2_to_cm2(A_i_m2):.3f} cm²")

    # ================================================
    # 🚀 NEW: APPLY HEADER AREA = 1.5 × OUTLET AREA
    # ================================================

    # desired minimum area at the end of manifold:
    # A_header = 1.5 * A_outlet
    A_outlet_cm2 = cfg["A_outlet_cm2"]
    A_required_min_cm2 = 1.5 * A_outlet_cm2

    # find current smallest plenum area
    A_min_current_cm2 = A_plenum[-1]

    # scaling factor for ALL plenum areas
    scale_factor = A_required_min_cm2 / A_min_current_cm2

    # scaled plenum areas
    A_plenum_scaled = [a * scale_factor for a in A_plenum]

    print("\n--- Scaled Plenum Geometry (1.5× header rule) ---")
    print(f"Required min header area : {A_required_min_cm2:.4f} cm²")
    print(f"Current min plenum area  : {A_min_current_cm2:.4f} cm²")
    print(f"Scale factor applied     : {scale_factor:.6f}")
    print(f"Largest area (scaled)    : {A_plenum_scaled[0]:.4f} cm²")
    print(f"Smallest area (scaled)   : {A_plenum_scaled[-1]:.4f} cm²")
    print("--------------------------------------------------\n")

    return outlet_velocities, outlet_massflows, A_plenum_scaled

# =============================
# 📊 PLOTTING FUNCTION
# =============================

def plot_half_ring(A_plenum):
    n_half = len(A_plenum)

    degrees_left = [-180 + (180 * i / n_half) for i in range(n_half)]
    areas_left = A_plenum[::-1]

    degrees_right = [0 + (180 * i / n_half) for i in range(n_half)]
    areas_right = A_plenum

    degrees = degrees_left + degrees_right
    areas = areas_left + areas_right

    plt.figure(figsize=(10,5))
    plt.plot(degrees, areas, marker='o', color='purple')
    plt.xlabel("Degrees around ring (-180° to 180°)")
    plt.ylabel("Plenum cross-sectional area (cm²)")
    plt.title("Half-Ring Manifold Plenum Area (mirrored)")
    plt.grid(True)
    plt.show()

# =============================
# ▶️ MAIN EXECUTION
# =============================

if __name__ == "__main__":
    v_out, m_out, A_plenum = ring_manifold_half(CONFIG)
    print("\nSimulation complete ✅\n")
    plot_half_ring(A_plenum)
