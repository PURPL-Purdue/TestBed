import math
import matplotlib.pyplot as plt

# =============================
# 🧰 USER CONFIGURATION
# =============================

CONFIG = {
    "num_outlets_total": 40,       # total number of outlet holes (full ring)
    "P_res_psi": 1000.0,           # fixed pump pressure, PSI
    "P_chamber_psi": 600.0,        # downstream chamber pressure, PSI
    "A_outlet_cm2": .150,           # outlet area (cm²)
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
    n_half = n_total // 2                # number of outlets in one half
    m_total = cfg["mass_flow_kg_s"]
    m_half = m_total / 2                 # mass flow in one half
    rho = cfg["rho"]
    A_outlet = cm2_to_m2(cfg["A_outlet_cm2"])
    P_res = psi_to_pa(cfg["P_res_psi"])
    P_chamber = psi_to_pa(cfg["P_chamber_psi"])

    # Target velocity for each outlet
    v_target = m_half / (rho * n_half * A_outlet)

    A_plenum = []
    outlet_velocities = []
    outlet_massflows = []

    for i in range(n_half):
        # Remaining mass flow in this half
        m_remain = m_half - i * (m_half / n_half)
        # Tapered plenum area
        A_i_m2 = m_remain / (rho * v_target)
        A_i_cm2 = m2_to_cm2(A_i_m2)
        A_plenum.append(A_i_cm2)

        outlet_velocities.append(v_target)
        outlet_massflows.append(m_half / n_half)

        if cfg["verbose"]:
            print(f"Outlet {i+1:02d}: "
                  f"{v_target:.3f} m/s | "
                  f"{m_half/n_half:.3f} kg/s | "
                  f"Plenum {A_i_cm2:.3f} cm²")

    return outlet_velocities, outlet_massflows, A_plenum

# =============================
# 📊 PLOTTING FUNCTION
# =============================

def plot_half_ring(A_plenum):
    n_half = len(A_plenum)
    # Degrees for negative side (-180 to 0), mirrored
    degrees_left = [-180 + (180 * i / n_half) for i in range(n_half)]
    areas_left = A_plenum[::-1]  # reverse for correct taper
    # Degrees for positive side (0 to 180)
    degrees_right = [0 + (180 * i / n_half) for i in range(n_half)]
    areas_right = A_plenum

    # Combine for full ring
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
    print("Note: this is one half of the ring manifold. The other half is symmetrical.")

    # Plot plenum area around the ring
    plot_half_ring(A_plenum)
