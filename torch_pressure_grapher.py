# Maelstrom Torch Testing Grapher
# Authors: Dominik Sloup
# First Created: 04/22/2023
# Last Updated: 12/02/2025
# Calculations done in SI units
import math
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from matplotlib.widgets import TextBox
# ──────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────
engine = CEA_Obj(oxName='GOX', fuelName='H2(g)')
gamma = 1.4
R_ox, R_h2 = 260.0, 4124.0 # J/(kg·K)
C_star = 2500 # m/s (CEA)
C_star_eff = 0.75
C_star_actual = 1896 # Adjusted for imperfect mixing
T0 = 293.15 # K
crit_ratio = (2/(gamma+1))**(gamma/(gamma-1))
psi_to_pa = 6894.76
max_feed = 850
# feed-pressure grid
Pox_psi = np.linspace(1, max_feed, 25)
Ph2_psi = np.linspace(1, max_feed, 25)
P_ox_pa, P_h2_pa = Pox_psi*psi_to_pa, Ph2_psi*psi_to_pa
X, Y = np.meshgrid(Pox_psi, Ph2_psi)
# choking coefficients
K_h = math.sqrt(gamma/(R_h2*T0)) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
K_o = math.sqrt(gamma/(R_ox*T0)) * (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
# ──────────────────────────────────────────────────────────────
# CORE CALCULATIONS
# ──────────────────────────────────────────────────────────────
def compute_fields(d_fuel, d_ox, d_exit):
    inch = 0.0254
    Ah = math.pi * (d_fuel * inch / 2)**2
    Ao = math.pi * (d_ox * inch / 2)**2
    Ae = math.pi * (d_exit * inch / 2)**2
    OF = (P_ox_pa[None, :] * Ao / np.sqrt(R_ox)) / (P_h2_pa[:, None] * Ah / np.sqrt(R_h2))
    C_star_arr = np.zeros_like(OF)
    for i in range(OF.shape[0]):
        for j in range(OF.shape[1]):
            C_star_arr[i, j] = engine.get_Cstar(Pc = 350, MR = OF[i, j])
    C_star_actual_arr = C_star_arr * C_star_eff
    Pc_pa = (Ao * P_ox_pa[None, :] * K_o + Ah * P_h2_pa[:, None] * K_h) / (Ae / C_star_actual_arr)
    Pc = Pc_pa / psi_to_pa
    mask = (Pc / X >= crit_ratio) | (Pc / Y >= crit_ratio)
    return Pc, OF, mask
def pick_levels(arr):
    return np.linspace(np.nanmin(arr), np.nanmax(arr), 5)[1:4] # 25–75 %
# ──────────────────────────────────────────────────────────────
# INITIAL VALUES
# ──────────────────────────────────────────────────────────────
d_fuel0, d_ox0, d_exit0 = 0.032, 0.28, 0.147
Pc0, OF0, mask0 = compute_fields(d_fuel0, d_ox0, d_exit0)
# ──────────────────────────────────────────────────────────────
# FIGURE & MAIN AXES
# ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 10), dpi=70)
main_ax = fig.add_axes([0.05, 0.42, 0.88, 0.50])
pcm = main_ax.pcolormesh(X, Y, Pc0, shading="auto")
cbar = fig.colorbar(pcm, ax=main_ax, label="P₍c,eff₎ (psi)")
cont_pc = main_ax.contour(X, Y, Pc0, levels=pick_levels(Pc0), colors='k', linewidths=1)
pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=8)
cont_of = main_ax.contour(X, Y, OF0, levels=[1,2,3,4,5], colors='k', linewidths=1)
of_labels = main_ax.clabel(cont_of, fmt='OF: %d', inline=True, fontsize=8)
point_handle = None
# ──────────────────────────────────────────────────────────────
# TEXT-BOX LAYOUT
# ──────────────────────────────────────────────────────────────
BOX_W, BOX_H = 0.1, 0.06
H_SP = 0.07
ROW1_Y, ROW2_Y = 0.25, 0.10
fig.canvas.draw()
left = min(main_ax.get_position().x0, cbar.ax.get_position().x0)
right = max(main_ax.get_position().x1, cbar.ax.get_position().x1)
X_START = (left + right - (3*BOX_W + 2*H_SP)) / 2 - 0.1
def make_box(x, y, label, init=''):
    ax_box = fig.add_axes([x, y, BOX_W, BOX_H])
    ax_box.text(0.5, 1.2, label, transform=ax_box.transAxes, ha='center', va='bottom')
    return TextBox(ax_box, '', initial=init)
# Row 1 – diameters
tb_h = make_box(X_START, ROW1_Y, 'Fuel Orifice diameter (in)', str(d_fuel0))
tb_o = make_box(X_START+BOX_W+H_SP, ROW1_Y, 'Ox Orifice Diameter (in)', str(d_ox0))
tb_e = make_box(X_START+2*(BOX_W+H_SP), ROW1_Y, 'Throat Diameter (in)', str(d_exit0))
# Row 2 – pressures
pf_box = make_box(X_START, ROW2_Y, 'Fuel Line Pressure (psi)')
po_box = make_box(X_START+BOX_W+H_SP, ROW2_Y, 'Oxidizer Line Pressure (psi)')
pc_box = make_box(X_START+2*(BOX_W+H_SP), ROW2_Y, 'Chamber Pressure (psi)')
of_box = make_box(X_START+3*(BOX_W+H_SP), ROW2_Y, 'OF Ratio (N/A)')
of_text = fig.text(X_START+BOX_W,0.04,'OF = ', va='center', fontsize=12, weight='bold')
# ──────────────────────────────────────────────────────────────
# CALLBACKS
# ──────────────────────────────────────────────────────────────
def compute_inlet(_):
    global point_handle
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    inch = 0.0254
    Ah, Ao, Ae = [math.pi*(d*inch/2)**2 for d in (dh, dox, de)]
    pf = float(pf_box.text) if pf_box.text else None
    po = float(po_box.text) if po_box.text else None
    pc = float(pc_box.text) if pc_box.text else None
    of = float(of_box.text) if of_box.text else None
    if pf is not None and po is not None:
        pc = (Ao*K_o*po + Ah*K_h*pf) * (C_star_actual/Ae)
        of_val = (po*Ao/np.sqrt(R_ox)) / (pf*Ah/np.sqrt(R_h2))
        pc_box.set_val(f'{pc:.2f}')
        of_box.set_val(f'{of_val:.2f}')
    elif pf is not None and pc is not None:
        po = (pc*(Ae/C_star_actual) - Ah*K_h*pf) / (Ao*K_o)
        of_val = (po*Ao/np.sqrt(R_ox)) / (pf*Ah/np.sqrt(R_h2))
        of_box.set_val(f'{of_val:.2f}')
        po_box.set_val(f'{po:.2f}')
    elif po is not None and pc is not None:
        pf = (pc*(Ae/C_star_actual) - Ao*K_o*po) / (Ah*K_h)
        of_val = (po*Ao/np.sqrt(R_ox)) / (pf*Ah/np.sqrt(R_h2))
        pf_box.set_val(f'{pf:.2f}')
        of_box.set_val(f'{of_val:.2f}')
    elif of is not None and pc is not None:
        po = (pc * (Ae / C_star_actual)) / (Ao * K_o * ((of + 1) / of))
        pf = (Ao * K_o * po) / (Ah * K_h * of)
        pf_box.set_val(f'{pf:.2f}')
        po_box.set_val(f'{po:.2f}')
    # plot operating point
    if pf is not None and po is not None:
        if point_handle:
            point_handle.remove()
        point_handle, = main_ax.plot(po, pf, 'o', ms=10, color='white', mec='black', zorder=5)
        # compute OF ratio and unchoke thresholds
        of_val = (po*Ao/np.sqrt(R_ox)) / (pf*Ah/np.sqrt(R_h2))
        thresh_fuel = pf * crit_ratio
        thresh_ox = po * crit_ratio
        thresh_pc = min(thresh_fuel, thresh_ox)
        of_text.set_text(f'OF = {of_val:.2f} Red-line Pc = {thresh_pc:.1f} psi')
        plt.draw()
hat = None  # make sure this exists globally

def apply_hatching(mask):
    global hat
    # Remove previous hatching
    if hat is not None and hasattr(hat, 'collections'):
        for coll in hat.collections:
            try:
                coll.remove()
            except:
                pass
    
    # Draw new hatching
    hat = main_ax.contourf(X, Y, 
                           mask.astype(float), 
                           levels=[0.5, 1.5], 
                           colors='none', 
                           hatches=['////'], 
                           zorder=2)
    
    # Style the red border lines
    if hasattr(hat, 'collections'):
        for coll in hat.collections:
            coll.set_edgecolor('red')
            coll.set_linewidth(0.5)
def update(_):
    global cont_pc, pc_labels, cont_of, of_labels, hat
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    Pc, OF, mask = compute_fields(dh, dox, de)
    pcm.set_array(Pc.ravel())
    pcm.set_clim(Pc.min(), Pc.max())
    cbar.update_normal(pcm)
    # Remove old contours safely
    for coll in getattr(cont_pc, 'collections', []):
        coll.remove()
    for coll in getattr(cont_of, 'collections', []):
        coll.remove()
    for txt in (*pc_labels, *of_labels):
        txt.remove()
    # Redraw contours
    cont_pc = main_ax.contour(X, Y, Pc, levels=pick_levels(Pc), colors='k', linewidths=1)
    pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=8)
    cont_of = main_ax.contour(X, Y, OF, levels=[1,2,3,4,5], colors='k', linewidths=1)
    of_labels = main_ax.clabel(cont_of, fmt='OF: %d', inline=True, fontsize=8)
    apply_hatching(mask)
    plt.draw()
# hook up callbacks
for box in (tb_h, tb_o, tb_e):
    box.on_submit(update)
for box in (pf_box, po_box, pc_box, of_box):
    box.on_submit(compute_inlet)
update(None) # initial draw
apply_hatching(mask0) # initial hatching
plt.show()