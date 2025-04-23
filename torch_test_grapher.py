# Maelstrom Torch Testing Grapher
# Authors: Dominik Sloup
# First Created: 04/22/2023
# Last Updated: 04/23/2023
# Calculations done in SI units

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox

# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────
gamma = 1.4
R_ox, R_h2 = 260.0, 4124.0             # J/(kg·K)
C_star     = 2018.2                    # m/s  (CEA)
T0         = 300.0                     # K
crit_ratio = (2/(gamma+1))**(gamma/(gamma-1))
psi_to_pa  = 6894.76

# feed-pressure grid
Pox_psi = np.linspace(0, 500, 200)
Ph2_psi = np.linspace(0, 500, 200)
P_ox_pa, P_h2_pa = Pox_psi*psi_to_pa, Ph2_psi*psi_to_pa
X, Y = np.meshgrid(Pox_psi, Ph2_psi)

# choking coefficients
K_h = math.sqrt(gamma/(R_h2*T0)) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
K_o = math.sqrt(gamma/(R_ox*T0)) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))

# ──────────────────────────────────────────────────────────────
#  CORE CALCULATIONS
# ──────────────────────────────────────────────────────────────
def compute_fields(d_fuel, d_ox, d_exit):
    inch = 0.0254
    Ah = math.pi*(d_fuel*inch/2)**2
    Ao = math.pi*(d_ox  *inch/2)**2
    Ae = math.pi*(d_exit*inch/2)**2

    Pc_pa = (Ao*P_ox_pa[None,:]*K_o + Ah*P_h2_pa[:,None]*K_h) / (Ae/C_star)
    Pc    = Pc_pa / psi_to_pa
    OF    = (P_ox_pa[None,:]*Ao/np.sqrt(R_ox)) / (P_h2_pa[:,None]*Ah/np.sqrt(R_h2))
    mask  = (Pc/X >= crit_ratio) | (Pc/Y >= crit_ratio)
    return Pc, OF, mask

def pick_levels(arr):
    return np.linspace(np.nanmin(arr), np.nanmax(arr), 5)[1:4]      # 25–75 %

# ──────────────────────────────────────────────────────────────
#  INITIAL VALUES
# ──────────────────────────────────────────────────────────────
d_fuel0, d_ox0, d_exit0 = 0.028, 0.032, 0.118
Pc0, OF0, mask0 = compute_fields(d_fuel0, d_ox0, d_exit0)

# ──────────────────────────────────────────────────────────────
#  FIGURE & MAIN AXES  (≈1000 px wide, *not* maximised)
# ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 10), dpi=70)               # 10 in × 100 dpi ≈ 1000 px
main_ax = fig.add_axes([0.05, 0.42, 0.88, 0.50])          # [L, B, W, H]

pcm  = main_ax.pcolormesh(X, Y, Pc0, shading="auto")
cbar = fig.colorbar(pcm, ax=main_ax, label="P₍c,eff₎ (psi)")

cont_pc   = main_ax.contour(X, Y, Pc0, levels=pick_levels(Pc0),
                            colors='k', linewidths=1)
pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=8)

cont_of   = main_ax.contour(X, Y, OF0, levels=[3,4,5],
                            colors='k', linewidths=1)
of_labels = main_ax.clabel(cont_of, fmt='OF: %d', inline=True, fontsize=8)

hat = main_ax.contourf(X, Y, mask0, levels=[0.5,1], colors=['none'],
                       hatches=['////'], zorder=2)
for coll in hat.collections:
    coll.set(edgecolor='red', linewidth=0.5)

main_ax.set(xlim=(0,500), ylim=(0,500), aspect='equal',
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
pf_box = make_box(X_START,                  ROW2_Y, 'Fuel Line Pressure (psi)')
po_box = make_box(X_START+BOX_W+H_SP,       ROW2_Y, 'Oxidizer Line Pressure (psi)')
pc_box = make_box(X_START+2*(BOX_W+H_SP),   ROW2_Y, 'Chamber Pressure (psi)')

# OF read-out (bottom-most)
of_text = fig.text(X_START+BOX_W,0.04,'OF = ', va='center', fontsize=12, weight='bold')

# ──────────────────────────────────────────────────────────────
#  CALLBACKS
# ──────────────────────────────────────────────────────────────
def compute_inlet(_):
    """Solve third pressure, update marker/OF, and display unchoke threshold."""
    global point_handle
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    inch = 0.0254
    Ah, Ao, Ae = [math.pi*(d*inch/2)**2 for d in (dh, dox, de)]

    pf = float(pf_box.text) if pf_box.text else None
    po = float(po_box.text) if po_box.text else None
    pc = float(pc_box.text) if pc_box.text else None

    # solve for the missing pressure
    if pf is not None and po is not None:
        pc = (Ao*K_o*po + Ah*K_h*pf) * (C_star/Ae)
        pc_box.set_val(f'{pc:.2f}')
    elif pf is not None and pc is not None:
        po = (pc*(Ae/C_star) - Ah*K_h*pf) / (Ao*K_o)
        po_box.set_val(f'{po:.2f}')
    elif po is not None and pc is not None:
        pf = (pc*(Ae/C_star) - Ao*K_o*po) / (Ah*K_h)
        pf_box.set_val(f'{pf:.2f}')

    # plot the operating point
    if pf is not None and po is not None:
        if point_handle:
            point_handle.remove()
        point_handle, = main_ax.plot(po, pf, 'o', ms=10,
                                     color='white', mec='black', zorder=5)

        # compute OF ratio
        of_val = (po*Ao/np.sqrt(R_ox)) / (pf*Ah/np.sqrt(R_h2))

        # compute unchoke thresholds (Pc at which flow JUST unchokes)
        #	choked if P_up/Pc >= crit_ratio  →  unchokes when Pc > P_up/crit_ratio
        thresh_fuel = pf * crit_ratio
        thresh_ox   = po * crit_ratio
        thresh_pc   = min(thresh_fuel, thresh_ox)   # whichever orifice unchokes first

        # update text: OF and red-line Pc (psi)
        of_text.set_text(
            f'OF = {of_val:.2f}    '
            f'Red-line Pc = {thresh_pc:.1f} psi'
        )

        plt.draw()

def update(_):
    """Re-draw contours when any diameter changes."""
    global cont_pc, pc_labels, cont_of, of_labels, hat
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    Pc, OF, mask = compute_fields(dh, dox, de)

    pcm.set_array(Pc.ravel())
    pcm.set_clim(Pc.min(), Pc.max())
    cbar.update_normal(pcm)

    for grp in (cont_pc.collections, cont_of.collections, hat.collections):
        for coll in grp:
            coll.remove()
    for txt in (*pc_labels, *of_labels):
        txt.remove()

    cont_pc   = main_ax.contour(X, Y, Pc, levels=pick_levels(Pc),
                                colors='k', linewidths=1)
    pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=8)

    cont_of   = main_ax.contour(X, Y, OF, levels=[3,4,5],
                                colors='k', linewidths=1)
    of_labels = main_ax.clabel(cont_of, fmt='OF: %d', inline=True, fontsize=8)

    hat = main_ax.contourf(X, Y, mask, levels=[0.5,1],
                           colors=['none'], hatches=['////'], zorder=2)
    for coll in hat.collections:
        coll.set(edgecolor='red', linewidth=0.5)

    plt.draw()

# hook up callbacks
for box in (tb_h, tb_o, tb_e):
    box.on_submit(update)
for box in (pf_box, po_box, pc_box):
    box.on_submit(compute_inlet)

update(None)       # first draw
plt.show()