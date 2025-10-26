# nozzle_diagnostics_corrected.py
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt

# ---------- Constants & conversions ----------
g0 = 9.80665             # m/s^2
psi_to_Pa = 6894.76      # Pa per psi
ft_to_m = 0.3048
lbf_to_N = 4.44822

# ---------- User-tunable inputs ----------

#target_thrust_lbf = 5000.0          # lbf
while True:
    target_thrust_lbf = input("Enter target thrust in lbf: ")
    try:
        target_thrust_lbf = float(target_thrust_lbf)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

target_thrust_N = target_thrust_lbf * lbf_to_N
ambient_p_psi = 14.7                # psi (sea level)
pe_psi = ambient_p_psi              # exit pressure (psi)

#OF_min = 0.1                       # minimum O/F ratio for sweep
#OF_max = 3.0                       # maximum O/F ratio for sweep
#OF_points = 80                    # number of values for sweep
while True:
    OF_min = input("Enter minimum O/F ratio for sweep: ")
    try:
        OF_min = float(OF_min)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    OF_max = input("Enter maximum O/F ratio for sweep: ")
    try:
        OF_max = float(OF_max)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

while True:
    OF_points = input("Enter number of values for sweep: ")
    try:
        OF_points = int(OF_points)  # convert to int
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

#OF_target = 1                        # single point O/F to print
while True:
    OF_target = input("Enter single point O/F to print: ")
    try:
        OF_target = float(OF_target)  # convert to float
        break  # exit loop if successful
    except ValueError:
        print("Invalid input, please enter a number.")

#pressures_to_check = [500, 600, 1000]    # chamber pressures [psi] to evaluate
while True:
    p_c_inputs = input("Enter 3 chamber pressures separated by commas, in psi. Ex: 500, 600, 1000: ")

    try:
        # Try converting all entries to floats
        pressures_to_check = [float(p.strip()) for p in p_c_inputs.split(",")]
        break  # success → exit loop
    except ValueError:
        print("Invalid input. Please enter only numbers, separated by commas.")


OF_range = np.linspace(OF_min, OF_max, OF_points)     # O/F sweep for plots

# ---------- Initialize CEA ----------
cea = CEA_Obj(oxName="LOX", fuelName="RP-1")

# ---------- CEA wrappers ----------
def get_gamma(pc_psi, of):
    _, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psi, MR=of)
    return gamma

def get_cstar_m_s(pc_psi, of):
    return cea.get_Cstar(Pc=pc_psi, MR=of) * ft_to_m

def get_Tc(pc_psi, of):
    temps_rankine = cea.get_Tcomb(Pc=pc_psi, MR=of)  # returns in Rankine
    temps_kelvin = (5/9) * temps_rankine
    return temps_kelvin
    #return cea.get_Tcomb(Pc=pc_psi, MR=of)  # in Rankine

def get_eps_at_PcOvPe(pc_psi, of, PcOvPe):
    return cea.get_eps_at_PcOvPe(Pc=pc_psi, MR=of, PcOvPe=PcOvPe)

# ---------- Core solver using exact CEA nozzle performance ----------
def compute_sizing(F_newtons, pc_psi, of, pe_psi=ambient_p_psi):
    gamma = get_gamma(pc_psi, of)
    cstar = get_cstar_m_s(pc_psi, of)
    Tc = get_Tc(pc_psi, of)

    # Calculate expansion ratio
    eps = get_eps_at_PcOvPe(pc_psi, of, pc_psi / pe_psi)

    # Calculate vacuum-specific impulse
    Cf, IspVac, IspSL = cea.get_PambCf(Pamb = 14.7, Pc = pc_psi, MR=of, eps =eps)
    Isp = (cstar * Cf) / 9.81

    # Calculate exhaust velocity (APPLY EFFICIENCY HERE)
    Ve = Isp * g0

    Pc_Pa = pc_psi * psi_to_Pa
    Pe_Pa = pe_psi * psi_to_Pa
    Pa_Pa = ambient_p_psi * psi_to_Pa

    # Mass flow from thrust equation (with pressure correction)
    denom = Ve + (Pe_Pa - Pa_Pa) * (cstar * eps) / Pc_Pa
    mdot = F_newtons / denom if denom > 0 else np.nan

    # Throat area directly from c* relation
    At = mdot * cstar / Pc_Pa
    Dt = np.sqrt(4.0 * At / np.pi)

    # Exit area and diameter
    Ae = At * eps
    De = np.sqrt(4.0 * Ae / np.pi)

    return {
        'pc_psi': pc_psi,
        'of': of,
        'gamma': gamma,
        'Tc': Tc,
        'cstar': cstar,
        'Isp': Isp,
        'Ve': Ve,
        'eps': eps,
        'Pe_psi': pe_psi,
        'mdot': mdot,
        'At': At,
        'Dt': Dt,
        'Ae': Ae,
        'De': De
    }

# ---------- Pretty print ----------
def pretty_print(r):
    print("="*70)
    print(f"Pc = {r['pc_psi']} psi | O/F = {r['of']:.3f}")
    print(f" gamma = {r['gamma']:.4f} | eps (Ae/At) = {r['eps']:.3f}")
    print(f" Tc = {r['Tc']:.1f} K")
    print(f" Isp = {r['Isp']:.2f} s | Ve = {r['Ve']:.1f} m/s | c* = {r['cstar']:.1f} m/s")
    print(f" mdot = {r['mdot']:.3f} kg/s")
    print(f" At = {r['At']:.6f} m^2 | Dt = {r['Dt']*100:.2f} cm ({r['Dt']*39.37:.2f} in)")
    print(f" Ae = {r['Ae']:.6f} m^2 | De = {r['De']*100:.2f} cm ({r['De']*39.37:.2f} in)")
    print("="*70 + "\n")

# ---------- Run single-point diagnostics ----------
if __name__ == "__main__":
    print("=== Single-point diagnostics ===\n")
    for pc in pressures_to_check:
        r = compute_sizing(target_thrust_N, pc, OF_target, pe_psi)
        pretty_print(r)

    # ---------- Prepare arrays for plots ----------
    fig, axs = plt.subplots(4, 2, figsize=(14, 14))
    axs = axs.flatten()

    for pc in pressures_to_check:
        Ve_list, mdot_list, cstar_list, gamma_list, eps_list = [], [], [], [], []
        At_list, Ae_list, Dt_list, De_list, Tc_list = [], [], [], [], []

        for of in OF_range:
            r = compute_sizing(target_thrust_N, pc, of, pe_psi)
            Ve_list.append(r['Ve'])
            mdot_list.append(r['mdot'])
            cstar_list.append(r['cstar'])
            gamma_list.append(r['gamma'])
            eps_list.append(r['eps'])
            At_list.append(r['At'])
            Ae_list.append(r['Ae'])
            Dt_list.append(r['Dt']*100)  # cm
            De_list.append(r['De']*100)  # cm
            Tc_list.append(r['Tc'])

        axs[0].plot(OF_range, Ve_list, label=f"Pc={pc} psi")
        axs[0].set_title("Exhaust Velocity Ve (m/s)")

        axs[1].plot(OF_range, mdot_list, label=f"Pc={pc} psi")
        axs[1].set_title("Mass Flow mdot (kg/s)")

        axs[2].plot(OF_range, cstar_list, label=f"Pc={pc} psi")
        axs[2].set_title("c* (m/s)")

        axs[3].plot(OF_range, gamma_list, label=f"Pc={pc} psi")
        axs[3].set_title("Gamma (γ)")

        axs[4].plot(OF_range, np.array(At_list)*1e4, label=f"At (cm²), Pc={pc}")
        axs[4].plot(OF_range, np.array(Ae_list)*1e4, '--', label=f"Ae (cm²), Pc={pc}")
        axs[4].set_title("Areas (At & Ae) vs O/F (cm²)")

        axs[5].plot(OF_range, Dt_list, label=f"Dt (cm), Pc={pc}")
        axs[5].plot(OF_range, De_list, '--', label=f"De (cm), Pc={pc}")
        axs[5].set_title("Diameters (Dt & De) vs O/F (cm)")

        axs[6].plot(OF_range, eps_list, label=f"Pc={pc} psi")
        axs[6].set_title("Expansion Ratio eps (Ae/At) vs O/F")

        axs[7].plot(OF_range, Tc_list, label=f"Pc={pc} psi")
        axs[7].set_title("Chamber Temperature Tc (K) vs O/F")

    for ax in axs:
        ax.set_xlabel("O/F")
        ax.grid(True)
        ax.legend()

    plt.tight_layout()
    plt.show()
