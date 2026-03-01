import math
import matplotlib.pyplot as plt

in_to_m = 0.0254
lb_to_kg = 0.4536
N_to_lbf = 0.2248
lbf_to_lbins = 386.09
CM2_TO_IN2 = 1 / 6.4516

# accepts metric but outputs inches
def calc_cv(mdot, inlet_rad, rho, orifice_num, deadCSA):
    inlet_ID = inlet_rad * 2
    inlet_area = math.pi * inlet_ID**2 / 4

    plenum_velocity = mdot / (rho*inlet_area)
    #distributed_mdot = mdot / lox_orifices_total

    N = orifice_num
    N_half = N // 2 

    mdot_orifice = mdot / N    
    mdot_in_half = mdot / 2

    A_half_cm2 = [0.0] * (N_half + 1)

    for k in range(N_half + 1):
        mdot_k = mdot_in_half - k * mdot_orifice
        if mdot_k < 0:
            mdot_k = 0.0
        A_m2 = mdot_k / (lox_rho * plenum_velocity)
        A_half_cm2[k] = (A_m2 * 10000) + 4  # m^2 -> cm^2
    
    A_half_cm2 = [a + deadCSA for a in A_half_cm2]

    A_half_in2 = [a * CM2_TO_IN2 for a in A_half_cm2]

    return A_half_in2

# plots CSA vs ° around manifold
def plot_ring_half(A_half_in2, signal, filename):
    n = len(A_half_in2)
    if n < 2:
        raise ValueError("Need at least 2 points in A_half_cm2 to plot.")

    # -180 to 0 to 180 
    deg_left  = [-180 + 180 * i / (n - 1) for i in range(n)]   # -180 .. 0
    deg_right = [   0 + 180 * i / (n - 1) for i in range(n)]   # 0 .. 180

    area_left  = A_half_in2[::-1]   # mirror
    area_right = A_half_in2

    # avoid duplicating 0 
    deg = deg_left[:-1] + deg_right
    area = area_left[:-1] + area_right

    plt.figure(figsize=(10, 5))
    plt.plot(deg, area, marker='o', linewidth=2)
    plt.xlabel("Degrees Around Manifold [-180° to 180°]", fontsize=12)
    plt.ylabel("Manifold Cross-Sectional Area [in^2]", fontsize=12)
    if signal == 1:
        plt.title("CV LOx Ring Manifold Area vs Angle", fontsize=14)
    else:
        plt.title("CV RP-1 Ring Manifold Area vs Angle", fontsize=14)
    
    plt.grid(True)
    plt.tight_layout()
    #plt.savefig(filename, dpi=200)
    plt.savefig(filename, dpi=200)
    plt.close()
    return

lox_orifices_total = 12
rp1_orifies_total = 6 # chosen value to model

lox_rho = 1141
rp1_rho = 810

lox_deadCSA = 4.0 # 4cm minimmum CSA
rp1_deadCSA = lox_deadCSA # change?

Mf = 3.0736 * lb_to_kg
Mo = 6.7616 * lb_to_kg

lox_inlet_rad = 0.422 * in_to_m
rp1_inlet_rad = 0.305 * in_to_m 

A_half_ox = calc_cv(Mo, lox_inlet_rad, lox_rho, lox_orifices_total, lox_deadCSA)
A_half_rp1 = calc_cv(Mf, rp1_inlet_rad, rp1_rho, rp1_orifies_total, rp1_deadCSA)

print(A_half_ox, A_half_rp1)

plot_ring_half(A_half_ox, 1, "lox_ring_area.png")
plot_ring_half(A_half_rp1, 2, "rp1_ring_area.png")

print(A_half_ox, A_half_rp1)

#print(A_half_ox, A_half_rp1)

# 24 orifices on pintle
# 12 rotating orifices for fuel inlet

# 0.844 - LOx hole diameter
# 1.3125 bore dia

# 0.61 - RP-1 hole diameter
# 1.0625 bore dia