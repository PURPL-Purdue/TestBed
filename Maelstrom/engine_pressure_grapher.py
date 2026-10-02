# Maelstrom Torch Testing Grapher
# Authors: Dominik Sloup
# First Created: 04/22/2023
# Last Updated: 04/25/2023
# Calculations done in SI units

import math
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.widgets import TextBox # type: ignore
import CEA_Wrap as CEA # type: ignore

# cstar = CEA.RocketProblem(materials=[CEA.Fuel("C2H5OH"), CEA.Oxidizer("O2")], pressure = 350, o_f = 0.28, analysis_type = "frozen").run_cea().cstar
# import matlab.engine as m

# try:
#    m = m.connect_matlab()
#except m.EngineError:
#    print("Matlab engine not running. Starting a new session...")
#    m = m.start_matlab()

#m.addpath(r'/Users/dominiksloup/Documents/MATLAB/CEA/')
    
# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────

gamma_ox = 1.4
gamma_rp = 1.22
R = 8.314                               # J /(mol·K)
mm_ox, mm_eth = 32.00, 46.07            # g/mol
R_ox = R / (mm_ox/1000)                # J/(kg·K)
R_eth = R / (mm_eth/1000)              # J/(kg·K)
C_star_raw = 1048.6                  # 0.8 * CEA @ 0.1 OF
C_star_eff = 0.9 # C Star efficiency
C_star = C_star_raw * C_star_eff
crit_ratio_ox = (2/(gamma_ox+1))**(gamma_ox/(gamma_ox-1))
fuel = CEA.Fuel("Jet-A(L)")
ox = CEA.Oxidizer("O2")
rho_fu = 810 #kg/m^3

relaxation = 0.05
psi_to_pa  = 6894.76
max_stiff = 0

cont_pc = None
cont_of = None
hat = None
hat2 = None
pc_labels = []
of_labels = []
num_orifices_ox = 12
num_orifices_fu = 24
global T0
T0 = 298 #

# feed-pressure grid
Pox_psi = np.linspace(1, 800, num=200)
Prp_psi = np.linspace(1, 800, num=200)
P_ox_pa, P_rp_pa = Pox_psi*psi_to_pa, Prp_psi*psi_to_pa
X, Y = np.meshgrid(Pox_psi, Prp_psi)

# choking coefficients
Cd_rp = 0.8   # estimate
Cd_ox = 0.8   #estimate
#K_rp = math.sqrt(gamma_rp/(R_rp*T0)) * (2/(gamma_rp+1))**((gamma_rp+1)/(2*(gamma_rp-1)))
K_o = math.sqrt(gamma_ox/(R_ox*T0)) * (2/(gamma_ox+1))**((gamma_ox+1)/(2*(gamma_ox-1)))

# ──────────────────────────────────────────────────────────────
#  CORE CALCULATIONS
# ──────────────────────────────────────────────────────────────
def Mdot_rp(Prp, Pc, Afu):
    """Mass flow rate through liquid orifice."""
    mdot = np.zeros_like(Pc)
    mdot = np.where(
        (Prp > Pc) & (Prp - Pc != 0),
        Cd_rp * Afu * np.sqrt(2 * rho_fu * (Prp - Pc)),
        0
    )
    return np.maximum(mdot, 0)

def Mdot_ox(Pox, Pc, Ao):
    """Mass flow rate through oxidizer orifice."""
    mdot = np.zeros_like(Pc)
    rho_ox = (Pox / (R_ox * T0))
    mdot = np.where((Pc / Pox > crit_ratio_ox), 
            Cd_ox * Ao * np.sqrt(2 * rho_ox *(Pox) * (gamma_ox/(gamma_ox-1)) * ((Pc/Pox)**(2/gamma_ox) - (Pc/ Pox)**((gamma_ox+1)/gamma_ox))),
            Cd_ox * Ao * Pox * K_o
        )
        
    return np.maximum(mdot, 0)

def compute_fields(Afu, Ao, Ae):
    error = []

    # Prepare arrays
    Pc_pa = np.full((len(Prp_psi), len(Pox_psi)), 1 * psi_to_pa)
    mdot_ox = Cd_ox * Ao * P_ox_pa[None, :] * K_o

    for i in range(200):  # Iterate for convergence
        mdot_rp = Mdot_rp(P_rp_pa[:, None], Pc_pa, Afu)
        mdot_ox = Mdot_ox(P_ox_pa[None, :], Pc_pa, Ao)
                                  
        Pc_pa_new = (mdot_ox + mdot_rp) * (C_star / Ae)
        # Under-relaxation to help convergence
        Pc_pa = (1-relaxation) * Pc_pa + relaxation * Pc_pa_new
        error.append((np.abs(Pc_pa[100,160] - Pc_pa_new[100,160])))

    mask_error = np.abs(Pc_pa - Pc_pa_new) > 0.5
    Pc = Pc_pa / psi_to_pa
    OF = np.where(mdot_rp != 0, mdot_ox / mdot_rp,100)
    
    vapor_p_rp = 2000 / psi_to_pa #psi
    rp_p_loss = (0.5*(rho_fu) * (mdot_rp/(rho_fu*Afu))**2)  / psi_to_pa   #ethanol dynamic pressure loss
    mask = (((X - Pc) / Pc) <= max_stiff) | (Y - rp_p_loss <= vapor_p_rp)  # unchoke threshold
    return Pc, OF, mask, mask_error
    
def pick_levels(arr):
    return np.linspace(np.nanmin(arr), np.nanmax(arr), 5)[1:4]      # 25–75 %

# ──────────────────────────────────────────────────────────────
#  INITIAL VALUES
# ──────────────────────────────────────────────────────────────
d_fuel0, d_ox0, d_exit0 = 0.047, 0.150, 0.93       # in
inch = 0.0254
Ae = math.pi * (d_exit0 * inch / 2) ** 2
Afu = num_orifices_fu * (math.pi * (d_fuel0 * inch / 2) ** 2)
Ao = num_orifices_ox * (math.pi * ((d_ox0 * inch / 2) ** 2))
Pc0, OF0, mask0, mask_error = compute_fields(Afu, Ao, Ae)


# ──────────────────────────────────────────────────────────────
#  FIGURE & MAIN AXES  (≈1000 px wide, *not* maximised)
# ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 10), dpi=70)               # 10 in × 100 dpi ≈ 1000 px
main_ax = fig.add_axes([0.05, 0.42, 0.88, 0.50])          # [L, B, W, H]

pcm = main_ax.pcolormesh(X, Y, Pc0, shading="auto", vmin=0, vmax=500)  # example limits
cbar = fig.colorbar(pcm, ax=main_ax, label="P₍c,eff₎ (psi)")

cont_pc   = main_ax.contour(X, Y, Pc0, levels = [150, 250, 350],
                            colors='k', linewidths=1)
pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=14)

cont_of   = main_ax.contour(X, Y, OF0, levels = [1, 1.5, 2, 2.5],
                            colors='k', linewidths=1)
of_labels = main_ax.clabel(cont_of, fmt='OF: %.2f', inline=True, fontsize=14)

hat1 = main_ax.contourf(X, Y, mask0, levels=[0.5,1.0], colors=['none'],
                       hatches=['\\\\'], zorder=2, edgecolor='red')

hat2 = main_ax.contourf(X, Y, mask_error, levels=[0.5,1.0], colors=['none'], 
                          hatches=['\\\\'], zorder=3)
                    
#print("Unique values in mask0:", np.unique(mask0))
#print(hat.collections)

if hasattr(hat, 'collections'):
    for coll in hat.collections:
       coll.set(edgecolor='red', linewidth=0.5)

main_ax.set(xlim=(0,800), ylim=(0,800), aspect='equal',
            xlabel='Oxidizer Feed Pressure (psi)',
            ylabel='Fuel Feed Pressure (psi)',
            title='Injector Choking Map  &  Chamber-Pressure Contours')

point_handle = None     # placeholder for marker

# ──────────────────────────────────────────────────────────────
#  TEXT-BOX LAYOUT
# ──────────────────────────────────────────────────────────────
BOX_W, BOX_H = 0.18, 0.06
H_SP         = 0.07
ROW1_Y, ROW2_Y = 0.25, 0.10            # ← widened vertical gap

# horizontal centring under plot+colour-bar
fig.canvas.draw()
left  = min(main_ax.get_position().x0, cbar.ax.get_position().x0)
right = max(main_ax.get_position().x1, cbar.ax.get_position().x1)
X_START = (left + right - (3*BOX_W + 2*H_SP)) / 2

def make_box(x, y, label, init=''):
    ax_box = fig.add_axes([x, y, BOX_W, BOX_H])
    ax_box.text(0.5, 1.2, label, transform=ax_box.transAxes,
                ha='center', va='bottom')
    return TextBox(ax_box, '', initial=init)

# row 1 – diameters
tb_h = make_box(X_START,                  ROW1_Y, 'Fuel Orifice diameter (in)',  str(d_fuel0))
tb_o = make_box(X_START+BOX_W+H_SP,       ROW1_Y, 'Ox Orifice Diameter (in)',    str(d_ox0))
tb_e = make_box(X_START+2*(BOX_W+H_SP),   ROW1_Y, 'Throat Diameter (in)',  str(d_exit0))

# row 2 – pressures
pf_box = make_box(X_START,                  ROW2_Y, 'Fuel Line Pressure (psi)', '')
po_box = make_box(X_START+BOX_W+H_SP,       ROW2_Y, 'Oxidizer Line Pressure (psi)', '')
pc_box = make_box(X_START+2*(BOX_W+H_SP),   ROW2_Y, 'Chamber Pressure (psi)',  '')

# OF read-out (bottom-most)
of_text = fig.text(X_START+BOX_W,0.04,'OF = ', va='center', fontsize=12, weight='bold')

# ──────────────────────────────────────────────────────────────
#  CALLBACKS
# ──────────────────────────────────────────────────────────────
def compute_inlet(_):
    global point_handle, Pc0, OF0

    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    # Pc0, OF0, mask0, mask_error = compute_fields(Afu, Ao, Ae)  # Recompute everything

    pf = float(pf_box.text) if pf_box.text else None
    po = float(po_box.text) if po_box.text else None

    if pf is not None and po is not None:
        pc = Pc0[int(pf/4), int(po/4)]
        pc_box.set_val(f'{pc:.2f}')

        if point_handle:
            point_handle.remove()

        point_handle = main_ax.plot(po, pf, 'o', ms=10,
                                     color='white', mec='black', zorder=5)[0]

        of_val = OF0[int(pf/4), int(po/4)]

        thresh_ox = po * crit_ratio_ox
        of_text.set_text(
            f'OF = {of_val:.2f}    '
            f'Red-line Pc = {thresh_ox:.1f} psi'
        )
        plt.draw()


def update(_):
    """Re-draw contours when any diameter changes."""
    global cont_pc, pc_labels, cont_of, of_labels, hat, hat2
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    Pc, OF, mask, mask_error = compute_fields(Afu, Ao, Ae)

    # pcm.set_array(Pc.ravel())
    # pcm.set_clim(Pc.min(), Pc.max())
    # cbar.update_normal(pcm)

    for grp in (cont_pc, cont_of, hat, hat2):
        if grp and hasattr(grp, 'collections'):
            for coll in grp.collections:
                coll.remove()

    for txt in (*pc_labels, *of_labels):
        txt.remove()

    cont_pc   = main_ax.contour(X, Y, Pc, levels=[150,250,350],
                                colors='k', linewidths=1)
    pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=10)

    cont_of   = main_ax.contour(X, Y, OF, levels=[1, 1.5, 2, 2.5],
                                colors='k', linewidths=1)
    of_labels = main_ax.clabel(cont_of, fmt='OF: %.2f', inline=True, fontsize=10)

    hat = main_ax.contourf(X, Y, mask, levels=[0.5,1],
                           colors=['none'], hatches=['\\\\'], zorder=2)
    
    hat2 = main_ax.contourf(X, Y, mask_error, levels=[0.5,1],
                            colors=['none'], hatches=['\\\\'], zorder=3)

    for grp in (hat, hat2):
        if hasattr(grp, 'collections'):
            for coll in grp.collections:
                coll.set(edgecolor='red', linewidth=0.5)

    plt.draw()

# hook up callbacks
for box in (tb_h, tb_o, tb_e):
    box.on_submit(update)
for box in (pf_box, po_box, pc_box):
    box.on_submit(compute_inlet)
    #print("flag5")


update(None)       # first draw
plt.show()