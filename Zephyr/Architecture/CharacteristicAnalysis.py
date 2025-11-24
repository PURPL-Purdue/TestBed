import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.widgets import Slider
from rocketcea.cea_obj import CEA_Obj

# =====================================================
# === USER CONFIGURATION ===============================
# =====================================================

################################################################
# INSTRUCTIONS: Run script to see performance and sizing, once you find an acceptable O/F ratio
# and chamber pressure, either change the reference point to that value, or enter those values into ZephyrCharacteristicGenerator.py
################################################################

# --- COMPARISON REFERENCE POINT (PLEASE CHANGE) ---
comparisonPc = 600
comparisonOF = 2.2

OXIDIZER = "LOX"
FUEL = "RP1"

TARGET_THRUST_LBF = 5000.0
CHAMBER_PRESSURES = [600]
OF_TARGET = 2.2
OF_RANGE = np.arange(1, 2.5 + 0.01, 0.01)




AMBIENT_P_PSI = 14.7
EFFICIENCY_FACTOR = .85
HOTFIRE_SECONDS = 30

NOZZLE_HALF_ANGLE_DEG = 15.0

# --- L* ---
LSTAR = 1.1  # m

# =====================================================
# === CONSTANTS =======================================
# =====================================================

G0 = 9.80665
PSI_TO_PA = 6894.76
LBF_TO_N = 4.44822
R_UNIVERSAL = 8314.4621
KG_TO_LB = 2.20462
GAL_PER_LITER = 0.264172
L_PER_M3 = 1000.0
M3_TO_GAL = GAL_PER_LITER * L_PER_M3

DENSITY_RP1 = 810.0
DENSITY_LOX = 1140.0
OFSTEP = 0.1

# =====================================================
# === HELPER FUNCTIONS ================================
# =====================================================

def cstar_from_T_gamma_R(Tc, gamma, R_specific):
    term = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    return math.sqrt(R_specific * Tc) / (gamma ** 0.5 * term)

def pe_over_pc_from_M(M, gamma):
    return (1.0 + 0.5 * (gamma - 1.0) * M**2) ** (-gamma / (gamma - 1.0))

def area_ratio_from_M(M, gamma):
    term = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * M**2)
    exponent = (gamma + 1.0) / (2.0 * (gamma - 1.0))
    return (1.0 / M) * (term ** exponent)

def Mach_from_pe_pc(pe_pc, gamma):
    if pe_pc <= 0 or pe_pc >= 1:
        raise ValueError("pe_pc must be between 0 and 1.")
    term = (pe_pc) ** (-(gamma - 1.0) / gamma)
    M = math.sqrt((2.0 / (gamma - 1.0)) * (term - 1.0))
    return M

def compute_Cf_ideal(gamma, Pe_Pc, Pa_Pc, Ae_over_At):
    term1 = (2.0 * gamma**2) / (gamma - 1.0)
    exponent = (gamma + 1.0) / (gamma - 1.0)
    bracket = (2.0 / (gamma + 1.0)) ** exponent
    momentum = math.sqrt(term1 * bracket * (1.0 - Pe_Pc ** ((gamma - 1.0) / gamma)))
    pressure_term = (Pe_Pc - Pa_Pc) * Ae_over_At
    return momentum + pressure_term

# =====================================================
# === CORE CEA WRAPPER ================================
# =====================================================

def get_numbers_extended(OF_ratio, pc_psi, ox, fuel):
    cea = CEA_Obj(oxName=ox, fuelName=fuel)
    temps_rankine = cea.get_Temperatures(Pc=pc_psi, MR=OF_ratio, eps=1)
    Tc = temps_rankine[0] * 5.0 / 9.0
    mw_gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=OF_ratio)
    MolWt, gamma = mw_gamma
    R_specific = R_UNIVERSAL / MolWt
    cstar = cstar_from_T_gamma_R(Tc, gamma, R_specific)
    return {'cstar': cstar, 'Tc': Tc, 'gamma': gamma, 'MolWt': MolWt, 'R_specific': R_specific}

# =====================================================
# === GEOMETRY & PERFORMANCE ==========================
# =====================================================

def compute_sizing(F_newtons, pc_psi, of):
    data = get_numbers_extended(of, pc_psi, OXIDIZER, FUEL)
    gamma, cstar, Tc, R_specific = data['gamma'], data['cstar'], data['Tc'], data['R_specific']
    Pc_Pa = pc_psi * PSI_TO_PA
    Pa_Pa = AMBIENT_P_PSI * PSI_TO_PA
    Pa_over_Pc = Pa_Pa / Pc_Pa

    if Pa_over_Pc < 1.0:
        Me = Mach_from_pe_pc(Pa_over_Pc, gamma)
        eps = area_ratio_from_M(Me, gamma)
        Pe_Pc = pe_over_pc_from_M(Me, gamma)
    else:
        Me, eps, Pe_Pc = 1.0, 1.0, 1.0

    Cf = compute_Cf_ideal(gamma, Pe_Pc, Pa_over_Pc, eps)
    At = F_newtons / (Pc_Pa * Cf)
    Dt = math.sqrt(4.0 * At / math.pi)

    cstar_eff = cstar * EFFICIENCY_FACTOR
    mdot = Pc_Pa * At / cstar_eff

    Ae = At * eps
    De = math.sqrt(4.0 * Ae / math.pi)
    Isp = (Cf * cstar) / G0
    Ve = Isp * G0 * EFFICIENCY_FACTOR

    D_chamber = 2.0 * Dt
    V_chamber = LSTAR * At
    L_chamber = V_chamber / (math.pi/4 * D_chamber**2)

    L_converge = (D_chamber - Dt) / (2 * math.tan(math.radians(45)))
    V_converge = (1/3) * math.pi * L_converge * ((D_chamber/2)**2 + (D_chamber/2)*(Dt/2) + (Dt/2)**2)

    V_total = V_chamber + V_converge
    Lstar_with_converge = V_total / At

    R_throat = 0.382 * Dt
    L_throat = R_throat + 0.5 * Dt
    L_nozzle = (De - Dt) / (2 * math.tan(math.radians(NOZZLE_HALF_ANGLE_DEG)))

    A_cyl_wall = math.pi * D_chamber * L_chamber
    r1 = D_chamber / 2.0
    r2 = Dt / 2.0
    slant_converge = math.sqrt((r1 - r2)**2 + L_converge**2)
    A_converge_wall = math.pi * (r1 + r2) * slant_converge

    A_total_wall = A_cyl_wall + A_converge_wall

    return {
        'pc_psi': pc_psi, 'of': of, 'gamma': gamma, 'Tc': Tc,
        'cstar': cstar, 'Isp': Isp, 'Ve': Ve, 'Cf': Cf, 'eps': eps,
        'mdot': mdot, 'At': At, 'Dt': Dt, 'Ae': Ae, 'De': De,
        'D_chamber': D_chamber, 'L_chamber': L_chamber*100.0,
        'L_chamber_cyl': L_chamber, 'L_converge': L_converge,
        'L_chamber_full': L_chamber + L_converge, 'V_chamber': V_chamber,
        'Lstar': LSTAR, 'V_total': V_total, 'Lstar_with_converge': Lstar_with_converge,
        'L_throat': L_throat*100.0, 'L_nozzle': L_nozzle*100.0,
        'A_cyl_wall': A_cyl_wall * 10000,
        'A_converge_wall': A_converge_wall * 10000,
        'A_total_wall': A_total_wall * 10000
    }

# =====================================================
# === PRETTY PRINT ====================================
# =====================================================

def pretty_print(r):
    print("="*60)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" gamma = {r['gamma']:.4f} | eps = {r['eps']:.3f}")
    print(f" Chamber Temp = {r['Tc']:.1f} K")
    print(f" Isp = {r['Isp']:.2f} s | Ve = {r['Ve']:.1f} m/s | c* = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" Dia Throat = {r['Dt']*100:.2f} cm | Dia Exit = {r['De']*100:.2f} cm")
    print(f" Chamber D = {r['D_chamber']*100:.2f} cm | Chamber L = {r['L_chamber']:.2f} cm")
    print(f" Converging length = {r['L_converge']*100:.2f} cm")
    print(f" Full chamber length= {r['L_chamber_full']*100:.2f} cm")
    print(f" Chamber volume (with converging section) = {r['V_chamber']*1e6:.2f} cm³")
    print(f" L* (just cylinder) = {r['Lstar']:.3f} m | L* (with converge) = {r['Lstar_with_converge']:.3f} m")
    print(f" Throat length = {r['L_throat']:.2f} cm | Nozzle length = {r['L_nozzle']:.2f} cm")
    print(f" SA Cyl = {r['A_cyl_wall']:.2f} cm^2 | SA Cyl+converge = {r['A_total_wall']:.2f} cm^2")

    of, mdot, t_fire = r['of'], r['mdot'], HOTFIRE_SECONDS
    mdot_fuel = mdot / (1.0 + of)
    mdot_ox = mdot - mdot_fuel
    m_fuel_total = mdot_fuel * t_fire
    m_ox_total = mdot_ox * t_fire

    V_fuel_gal = (m_fuel_total / DENSITY_RP1) * M3_TO_GAL
    V_ox_gal = (m_ox_total / DENSITY_LOX) * M3_TO_GAL

    print("\n")
    print(f"--- Hotfire {HOTFIRE_SECONDS:.1f}s ---")
    print(f" Fuel: {m_fuel_total:.2f} kg ({V_fuel_gal:.2f} gal)")
    print(f" Oxidizer: {m_ox_total:.2f} kg ({V_ox_gal:.2f} gal)")
    print(f" Total: {(m_fuel_total+m_ox_total):.2f} kg ({V_fuel_gal+V_ox_gal:.2f} gal)")
    print("="*60 + "\n")

# =====================================================
# === INTERACTIVE SLIDER (FIXED RED DOTS) =============
# =====================================================

def interactive_OF_plot(pc_min=300, pc_max=700, pc_step=25, target_thrust_N=None):

    pc_values = np.arange(pc_min, pc_max + 1e-9, pc_step)

    # fixed reference
    ref = compute_sizing(target_thrust_N, comparisonPc, comparisonOF)
    Tc_ref = ref['Tc']
    Dch_ref = ref['D_chamber'] * 100.0
    Isp_ref = ref['Isp']
    mdot_ref = ref['mdot']
    mdot_fuel_ref = ref['mdot'] / (1.0 + comparisonOF)
    A_cyl_ref = ref['A_cyl_wall']
    A_total_ref = ref['A_total_wall']

    # Precompute curves
    data = {}
    for of in OF_RANGE:
        Tc_list = []
        Dch_list = []
        Isp_list = []
        mdot_list = []
        mdot_fuel_list = []
        A_cyl_list = []
        A_total_list = []

        for pc in pc_values:
            r = compute_sizing(target_thrust_N, pc, of)
            Tc_list.append(r['Tc'])
            Dch_list.append(r['D_chamber'] * 100.0)
            Isp_list.append(r['Isp'])
            mdot_list.append(r['mdot'])
            mdot_fuel_list.append(r['mdot'] / (1 + of))
            A_cyl_list.append(r['A_cyl_wall'])
            A_total_list.append(r['A_total_wall'])

        data[of] = {
            "Tc": Tc_list,
            "Dch": Dch_list,
            "Isp": Isp_list,
            "mdot": mdot_list,
            "mdot_fuel": mdot_fuel_list,
            "A_cyl": A_cyl_list,
            "A_total": A_total_list,
        }

    # Figure setup
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))
    plt.subplots_adjust(bottom=0.18)

    of0 = OF_TARGET
    of_key = min(OF_RANGE, key=lambda x: abs(x - of0))
    current = data[of_key]

    plot_keys = ["Tc", "Dch", "Isp", "A_cyl", "mdot", "mdot_fuel"]
    ylabels = [
        "Tc (K)", "Chamber Dia (cm)", "Isp (s)",
        "Surface Area (cm²)", "mdot (kg/s)", "fuel mdot (kg/s)"
    ]

    axes_list = axs.ravel()
    lines = []

    multi_curve_index = 3

    for i, (ax, key, ylabel) in enumerate(zip(axes_list, plot_keys, ylabels)):

        if i == multi_curve_index:
            ln1, = ax.plot(pc_values, current["A_cyl"], "o-", label="Cyl")
            ln2, = ax.plot(pc_values, current["A_total"], "s-", label="Cyl + Converge")
            lines.append((ln1, ln2))

            ax.plot(comparisonPc, A_cyl_ref, "ro", markersize=7)
            ax.plot(comparisonPc, A_total_ref, "ro", markersize=7)

            ax.legend()

        else:
            ln, = ax.plot(pc_values, current[key], "o-")
            lines.append(ln)

            # fixed dot
            if key == "Tc":
                ax.plot(comparisonPc, Tc_ref, "ro", markersize=7)
            elif key == "Dch":
                ax.plot(comparisonPc, Dch_ref, "ro", markersize=7)
            elif key == "Isp":
                ax.plot(comparisonPc, Isp_ref, "ro", markersize=7)
            elif key == "mdot":
                ax.plot(comparisonPc, mdot_ref, "ro", markersize=7)
            elif key == "mdot_fuel":
                ax.plot(comparisonPc, mdot_fuel_ref, "ro", markersize=7)

        ax.set_xlabel("Pc (psi)")
        ax.set_ylabel(ylabel)
        ax.grid(True)

        # === FIXED AXIS RANGES ===
        if key == "Tc":
            ax.set_ylim(1500, 4000)
        elif key == "Isp":
            ax.set_ylim(100, 300)
        elif key == "mdot":
            ax.set_ylim(8, 12)
        elif key == "mdot_fuel":
            ax.set_ylim(2, 8)

    # Slider
    ax_slider = plt.axes([0.15, 0.05, 0.7, 0.04])
    slider = Slider(ax_slider, "O/F", OF_RANGE.min(), OF_RANGE.max(), valinit=of0)


    def update(val):
        of_val = slider.val
        of_key = min(OF_RANGE, key=lambda x: abs(x - of_val))
        new = data[of_key]
        fig.suptitle(f"O/F = {of_key:.3f}", fontsize=16)

        for i, key in enumerate(plot_keys):

            if i == multi_curve_index:
                ln1, ln2 = lines[i]
                ln1.set_ydata(new["A_cyl"])
                ln2.set_ydata(new["A_total"])

            else:
                lines[i].set_ydata(new[key])



    slider.on_changed(update)
    plt.show()

# =====================================================
# === NEW: PERCENT-CHANGE FIGURE (REFERENCE O/F ONLY) ==
# =====================================================

def make_percent_change_figure(pc_min, pc_max, pc_step, target_thrust_N):

    # Build the Pc range to compare against
    pc_values = np.arange(pc_min, pc_max + 1e-9, pc_step)

    # --- Reference point ---
    ref_result = compute_sizing(target_thrust_N, comparisonPc, comparisonOF)

    ref_vals = {
        "Tc": ref_result["Tc"],
        "Dch": ref_result["D_chamber"] * 100.0,
        "Isp": ref_result["Isp"],
        "mdot": ref_result["mdot"],
        "mdot_fuel": ref_result["mdot"] / (1.0 + comparisonOF),
        "A_cyl": ref_result["A_cyl_wall"]
    }

    # --- Storage for percent-change curves ---
    percent_data = {
        "Tc": [],
        "Dch": [],
        "Isp": [],
        "mdot": [],
        "mdot_fuel": [],
        "A_cyl": []
    }

    # Compute percent changes vs Pc at REFERENCE OF ONLY
    for pc in pc_values:
        r = compute_sizing(target_thrust_N, pc, comparisonOF)

        current_vals = {
            "Tc": r["Tc"],
            "Dch": r["D_chamber"] * 100.0,
            "Isp": r["Isp"],
            "mdot": r["mdot"],
            "mdot_fuel": r["mdot"] / (1.0 + comparisonOF),
            "A_cyl": r["A_cyl_wall"]
        }

        # percent change formula
        for key in percent_data:
            pct = 100.0 * (current_vals[key] - ref_vals[key]) / ref_vals[key]
            percent_data[key].append(pct)

    # ------------------------------
    # Create the percent-change figure
    # ------------------------------
    fig2, axs2 = plt.subplots(2, 3, figsize=(16, 9))
    axes = axs2.ravel()

    keys = ["Tc", "Dch", "Isp", "mdot", "mdot_fuel", "A_cyl"]
    labels = [
        "Tc (%)",
        "Chamber Diameter (%)",
        "Isp (%)",
        "mdot (%)",
        "Fuel mdot (%)",
        "Surface Area Cyl (%)"
    ]

    # Determine GLOBAL y-limits to share across all
    all_values = np.array([v for arr in percent_data.values() for v in arr])
    ymin = np.min(all_values) * 1.05
    ymax = np.max(all_values) * 1.05

    for ax, key, label in zip(axes, keys, labels):
        ax.plot(pc_values, percent_data[key], "o-")
        ax.axhline(0, color="gray", linewidth=1)
        ax.set_ylabel(label)
        ax.set_xlabel("Pc (psi)")
        ax.grid(True)
        ax.set_ylim(ymin, ymax)

    fig2.suptitle(
        f"Percent Change vs Pc at Reference O/F = {comparisonOF}, Pc_ref = {comparisonPc} psi",
        fontsize=16
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# =====================================================
# === NEW: PERCENT-CHANGE TABLE =======================
# =====================================================

def print_percent_change_table(pc_min, pc_max, pc_step, target_thrust_N):

    pc_values = np.arange(pc_min, pc_max + 1e-9, pc_step)

    # Reference at comparisonPc, comparisonOF
    ref_result = compute_sizing(target_thrust_N, comparisonPc, comparisonOF)

    ref_vals = {
        "Tc": ref_result["Tc"],
        "Dch": ref_result["D_chamber"] * 100.0,
        "Isp": ref_result["Isp"],
        "mdot": ref_result["mdot"],
        "mdot_fuel": ref_result["mdot"] / (1.0 + comparisonOF),
        "A_cyl": ref_result["A_cyl_wall"]
    }

    # Header
    print("\n==================== PERCENT CHANGE TABLE ====================")
    print(f"Reference Point: Pc = {comparisonPc} psi,  O/F = {comparisonOF}")
    print("All values reported as % change relative to reference.\n")

    header = (
        f"{'Pc (psi)':>8} | "
        f"{'Tc (%)':>10} | {'Dch (%)':>10} | {'Isp (%)':>10} | "
        f"{'mdot (%)':>10} | {'fuel mdot (%)':>15} | {'A_cyl (%)':>12}"
    )
    print(header)
    print("-" * len(header))

    # Compute and print each row
    for pc in pc_values:
        r = compute_sizing(target_thrust_N, pc, comparisonOF)

        current_vals = {
            "Tc": r["Tc"],
            "Dch": r["D_chamber"] * 100.0,
            "Isp": r["Isp"],
            "mdot": r["mdot"],
            "mdot_fuel": r["mdot"] / (1.0 + comparisonOF),
            "A_cyl": r["A_cyl_wall"],
        }

        pct = {}
        for key in ref_vals:
            pct[key] = 100.0 * (current_vals[key] - ref_vals[key]) / ref_vals[key]

        print(
            f"{pc:8.1f} | "
            f"{pct['Tc']:10.3f} | {pct['Dch']:10.3f} | {pct['Isp']:10.3f} | "
            f"{pct['mdot']:10.3f} | {pct['mdot_fuel']:15.3f} | {pct['A_cyl']:12.3f}"
        )

    print("=" * len(header))
    print()


# =====================================================
# === UPDATED MAIN ====================================
# =====================================================

if __name__ == "__main__":
    print("=== Rocket Engine Analysis (auto-optimal chamber geometry) ===\n")

    target_thrust_N = TARGET_THRUST_LBF * LBF_TO_N

    # Print sizing for selected chamber pressures
    for pc in CHAMBER_PRESSURES:
        r = compute_sizing(target_thrust_N, pc, OF_TARGET)
        pretty_print(r)

    # --- Original interactive figure (unchanged) ---
    interactive_OF_plot(target_thrust_N=target_thrust_N)

    # --- NEW percent-change figure ---
    # Uses same Pc range as interactive plot (default 300–700 psi)
    make_percent_change_figure(
        pc_min=300,
        pc_max=600,
        pc_step=5,
        target_thrust_N=target_thrust_N
    )
        # --- NEW percent-change table ---
    print_percent_change_table(
        pc_min=300,
        pc_max=600,
        pc_step=5,
        target_thrust_N=target_thrust_N
    )


