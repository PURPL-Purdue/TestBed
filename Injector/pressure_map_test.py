import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import TextBox, Button
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
from ruamel.yaml import YAML

# ── optional import ──────────────────────────────────────────────────────────
try:
    from pressure_map_gasgas import pressure_map
    HAS_PRESSURE_MAP = True
except ImportError:
    HAS_PRESSURE_MAP = False
    print("WARNING: pressure_map_gasgas not found – using synthetic demo data.")

# ── YAML loader ──────────────────────────────────────────────────────────────
def find_yaml(filename="Maelstrom_torch.yaml", start_dir=None):
    start_dir = Path(start_dir or Path.cwd())
    for path in start_dir.rglob(filename):
        return path
    return None

# ── constants ────────────────────────────────────────────────────────────────
ft_to_m  = 0.3048
m_to_in  = 39.3701
psi_to_pa = 6894.76

# ── helpers ──────────────────────────────────────────────────────────────────
def finite_minmax(a):
    a = np.asarray(a, dtype=float)
    m = np.isfinite(a)
    if not np.any(m):
        raise ValueError("Array has no finite values.")
    return float(np.min(a[m])), float(np.max(a[m]))


# ── interpolator builder ─────────────────────────────────────────────────────
def build_interpolators(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo):
    """Return (fwd_pc, fwd_of, inv_pfuel, inv_pox) interpolators."""
    pc_ax = np.asarray(pc_scale, float)
    of_ax = np.asarray(OF_scale, float)

    # Forward: (Pc, OF) → feed pressures
    fwd_pfuel = RegularGridInterpolator(
        (pc_ax, of_ax), fu_p_map, method="linear", bounds_error=False, fill_value=np.nan)
    fwd_pox = RegularGridInterpolator(
        (pc_ax, of_ax), ox_p_map, method="linear", bounds_error=False, fill_value=np.nan)

    # Inverse: (p_fuel, p_ox) → Pc / OF
    pf_min, pf_max = finite_minmax(fu_p_map)
    po_min, po_max = finite_minmax(ox_p_map)
    n_pf, n_po = PC_pfpo.shape
    pf_ax = np.linspace(pf_min, pf_max, n_pf)
    po_ax = np.linspace(po_min, po_max, n_po)

    inv_pc = RegularGridInterpolator(
        (pf_ax, po_ax), PC_pfpo, method="linear", bounds_error=False, fill_value=np.nan)
    inv_of = RegularGridInterpolator(
        (pf_ax, po_ax), OF_pfpo, method="linear", bounds_error=False, fill_value=np.nan)

    return fwd_pfuel, fwd_pox, inv_pc, inv_of, (pf_ax, po_ax, pc_ax, of_ax)


# ── main plot + interactive app ──────────────────────────────────────────────
def launch_app(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo, max_pc=250):

    interps = build_interpolators(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo)
    fwd_pfuel, fwd_pox, inv_pc, inv_of, axes = interps
    pf_ax, po_ax, pc_ax, of_ax = axes

    pf_min, pf_max = pf_ax[0], pf_ax[-1]
    po_min, po_max = po_ax[0], po_ax[-1]

    # ── figure layout ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(13, 11))

    # GridSpec: top = heatmap, bottom = control panel
    gs = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[3, 1],
        hspace=0.15,
    )
    ax_map  = fig.add_subplot(gs[0])
    ax_ctrl = fig.add_subplot(gs[1])
    ax_ctrl.set_facecolor("#f0f0f0")
    ax_ctrl.set_xticks([])
    ax_ctrl.set_yticks([])
    for spine in ax_ctrl.spines.values():
        spine.set_edgecolor("#cccccc")

    # ── heatmap ──────────────────────────────────────────────────────────────
    im = ax_map.imshow(
        PC_pfpo,
        origin="lower", aspect="auto",
        extent=[po_min, po_max, pf_min, pf_max],
        cmap="viridis", vmin=0, vmax=max_pc,
    )
    cb = fig.colorbar(im, ax=ax_map, pad=0.01)
    cb.set_label("Pc (psi)", fontsize=11, labelpad=8)

    # OF contours
    PO, PF = np.meshgrid(po_ax, pf_ax)
    of_finite = OF_pfpo[np.isfinite(OF_pfpo)]
    of_levels = np.linspace(np.percentile(of_finite, 5), np.percentile(of_finite, 95), 7) \
        if of_finite.size else np.linspace(1.5, 5.5, 7)
    cs = ax_map.contour(PO, PF, OF_pfpo, levels=of_levels, colors="k", linewidths=0.8)
    ax_map.clabel(cs, fmt=lambda v: f"OF={v:.2f}", inline=True, fontsize=8)

    # Pc constant-pressure contours (dashed grey)
    pc_levels = np.arange(25, max_pc + 1, 25)
    pc_cs = ax_map.contour(PO, PF, PC_pfpo, levels=pc_levels,
                           colors="grey", linewidths=0.6, linestyles="dashed", alpha=1)
    ax_map.clabel(pc_cs, fmt=lambda v: f"{v:.0f} psi", inline=True, fontsize=7, colors="k")

    # S26 Maelstrom operating envelope: OF in [2.15, 2.75], Pc in [150, 175] psi
    OF_lo, OF_hi = 2.15, 2.75
    PC_lo, PC_hi = 150.0, 175.0
    envelope_mask = (
        (OF_pfpo >= OF_lo) & (OF_pfpo <= OF_hi) &
        (PC_pfpo >= PC_lo) & (PC_pfpo <= PC_hi)
    ).astype(float)
    ax_map.contourf(PO, PF, envelope_mask, levels=[0.5, 1.5],
                    colors=["#8b0000"], alpha=0.0, zorder=12)
    hatch_cs = ax_map.contourf(PO, PF, envelope_mask, levels=[0.5, 1.5],
                               hatches=["/////"], colors="none",
                               edgecolors="#8b0000", linewidths=0.5, zorder=12)
    # also draw the boundary
    from matplotlib.lines import Line2D
    envelope_handle = Line2D([0], [0], color="#8b0000", linewidth=2.2,
                             linestyle="dashed", label="S26 Maelstrom")

    ax_map.set_xlabel("p_ox (psi)", fontsize=11)
    ax_map.set_ylabel("p_fu (psi)", fontsize=11)
    ax_map.set_title("Torch Igniter Operating Envelope", fontsize=12)

    # ── hardcoded campaign points ────────────────────────────────────────────
    # format: (p_fuel, p_ox)  →  plotted as (p_ox, p_fuel) on the map
    old_campaign = [
        (470, 390),
        (440, 360),
        (390, 380),
        (460, 275),
        (460, 245),
    ]
    new_campaign = [
        (686, 260),
        (590, 260),
        (480, 260),
        (430, 240),
        (350, 290),
        (400, 270),
    ]

    old_ox, old_fu = zip(*[(po, pf) for pf, po in old_campaign])
    new_ox, new_fu = zip(*[(po, pf) for pf, po in new_campaign])

    ax_map.plot(old_ox, old_fu, "s", ms=9, zorder=11,
                markerfacecolor="none", markeredgecolor="black", markeredgewidth=1.8,
                label="S25 Torch + GG")
    ax_map.plot(new_ox, new_fu, "o", ms=9, zorder=11,
                markerfacecolor="black", markeredgecolor="black", markeredgewidth=1.0,
                label="S26 Torch Only")
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D as _L2D
    sq_handle = _L2D([0],[0], marker='s', color='w', markerfacecolor='none',
                     markeredgecolor='black', markeredgewidth=1.8, markersize=9, label='Old campaign')
    ci_handle = _L2D([0],[0], marker='o', color='w', markerfacecolor='black',
                     markeredgecolor='black', markersize=9, label='New campaign')
    ax_map.legend(handles=[sq_handle, ci_handle, envelope_handle],
                  loc='upper left', fontsize=9, framealpha=0.85)

    # accumulating points — blue for feeds→output, orange for output→feeds
    pts_a, = ax_map.plot([], [], "o", color="#4c8bf5", ms=7, zorder=10,
                         markeredgecolor="white", markeredgewidth=0.8)
    pts_b, = ax_map.plot([], [], "o", color="#e8622a", ms=7, zorder=10,
                         markeredgecolor="white", markeredgewidth=0.8)
    pts_a_data = []   # list of (p_ox, p_fuel)
    pts_b_data = []

    # ── control panel layout ─────────────────────────────────────────────────
    #
    #   MODE A  [p_fuel] [p_ox]  [LOOKUP →]   |   MODE B  [Pc] [OF]  [← LOOKUP]
    #
    PANEL_BG   = "#f0f0f0"
    LABEL_CLR  = "#333333"
    INPUT_BG   = "white"
    INPUT_FG   = "#111111"
    BTN_A_CLR  = "#4c8bf5"
    BTN_B_CLR  = "#e8622a"
    RESULT_CLR = "#111111"
    MONO       = "monospace"

    def make_textbox(rect, init=""):
        """TextBox with no built-in label (we draw labels manually)."""
        ax_tb = fig.add_axes(rect, facecolor=INPUT_BG)
        tb = TextBox(ax_tb, "", initial=init,
                     color=INPUT_BG, hovercolor="#e8e8e8",
                     label_pad=0)
        tb.text_disp.set_color(INPUT_FG)
        tb.text_disp.set_fontfamily(MONO)
        tb.text_disp.set_fontsize(11)
        return tb

    def make_button(rect, label, color):
        ax_bt = fig.add_axes(rect, facecolor=color)
        bt = Button(ax_bt, label, color=color, hovercolor="#555555")
        bt.label.set_color("white")
        bt.label.set_fontfamily(MONO)
        bt.label.set_fontsize(10)
        bt.label.set_fontweight("bold")
        return bt

    # ── layout math ──────────────────────────────────────────────────────────
    ctrlb = ax_ctrl.get_position()
    L, B, W, H = ctrlb.x0, ctrlb.y0, ctrlb.width, ctrlb.height

    def fax(xl, xr, yt, yb):
        """ctrl-panel relative [0,1] → figure fraction."""
        return [L + xl*W, B + yb*H, (xr-xl)*W, (yt-yb)*H]

    def ctrl_label(x_rel, y_rel, txt, **kwargs):
        """Place a label in ctrl-panel relative coords via ax_ctrl.transAxes."""
        ax_ctrl.text(x_rel, y_rel, txt, transform=ax_ctrl.transAxes,
                     fontfamily=MONO, fontsize=9, color=LABEL_CLR,
                     ha="center", va="bottom", **kwargs)

    # ── divider & section headers ─────────────────────────────────────────────
    ax_ctrl.axvline(0.5, color="#bbbbbb", linewidth=1.5)

    for x, txt, clr in [
        (0.245, "FEEDS → OUTPUT", BTN_A_CLR),
        (0.745, "OUTPUT → FEEDS", BTN_B_CLR),
    ]:
        ax_ctrl.text(x, 0.97, txt, transform=ax_ctrl.transAxes,
                     color=clr, fontsize=9, fontfamily=MONO, fontweight="bold",
                     ha="center", va="top")

    # ── result text objects ───────────────────────────────────────────────────
    res_a = ax_ctrl.text(0.245, 0.04, "", transform=ax_ctrl.transAxes,
                         color=BTN_A_CLR, fontsize=10, fontfamily=MONO, ha="center", va="bottom")
    res_b = ax_ctrl.text(0.745, 0.04, "", transform=ax_ctrl.transAxes,
                         color=BTN_B_CLR, fontsize=10, fontfamily=MONO, ha="center", va="bottom")

    # ── Mode A: feeds → Pc, OF ───────────────────────────────────────────────
    #   boxes sit in the vertical middle band: yb=0.30, yt=0.68
    tb_pfuel = make_textbox(fax(0.02, 0.17, 0.68, 0.32))
    tb_pox   = make_textbox(fax(0.20, 0.35, 0.68, 0.32))
    btn_a    = make_button( fax(0.37, 0.48, 0.68, 0.32), "LOOKUP →", BTN_A_CLR)

    # labels above each box
    ctrl_label(0.095, 0.69, "p_fuel (psi)")
    ctrl_label(0.275, 0.69, "p_ox (psi)")

    # ── Mode B: Pc, OF → feeds ───────────────────────────────────────────────
    tb_pc = make_textbox(fax(0.52, 0.67, 0.68, 0.32))
    tb_of = make_textbox(fax(0.70, 0.85, 0.68, 0.32))
    btn_b = make_button( fax(0.87, 0.98, 0.68, 0.32), "← LOOKUP", BTN_B_CLR)

    ctrl_label(0.595, 0.69, "Pc (psi)")
    ctrl_label(0.775, 0.69, "OF (–)")

    # clear button — centred at the bottom of the divider
    btn_clear = make_button(fax(0.44, 0.56, 0.38, 0.08), "CLEAR ALL", "#888888")

    # ── callbacks ─────────────────────────────────────────────────────────────
    def add_point_a(p_ox_val, p_fuel_val):
        pts_a_data.append((p_ox_val, p_fuel_val))
        xs, ys = zip(*pts_a_data)
        pts_a.set_data(xs, ys)
        fig.canvas.draw_idle()

    def add_point_b(p_ox_val, p_fuel_val):
        pts_b_data.append((p_ox_val, p_fuel_val))
        xs, ys = zip(*pts_b_data)
        pts_b.set_data(xs, ys)
        fig.canvas.draw_idle()

    def on_lookup_a(_):
        res_a.set_text("")
        try:
            pf = float(tb_pfuel.text)
            po = float(tb_pox.text)
        except ValueError:
            res_a.set_color("#ff4455")
            res_a.set_text("⚠ enter numeric values")
            fig.canvas.draw_idle()
            return

        pc_val = float(inv_pc([[pf, po]]))
        of_val = float(inv_of([[pf, po]]))

        if not (np.isfinite(pc_val) and np.isfinite(of_val)):
            res_a.set_color("#ff4455")
            res_a.set_text("⚠ outside map range")
        else:
            res_a.set_color("#4c8bf5")
            res_a.set_text(f"Pc = {pc_val:.1f} psi   |   OF = {of_val:.3f}")
            add_point_a(po, pf)

        fig.canvas.draw_idle()

    def on_lookup_b(_):
        res_b.set_text("")
        try:
            pc  = float(tb_pc.text)
            of  = float(tb_of.text)
        except ValueError:
            res_b.set_color("#ff4455")
            res_b.set_text("⚠ enter numeric values")
            fig.canvas.draw_idle()
            return

        pf_val = float(fwd_pfuel([[pc, of]]))
        po_val = float(fwd_pox([[pc, of]]))

        if not (np.isfinite(pf_val) and np.isfinite(po_val)):
            res_b.set_color("#ff4455")
            res_b.set_text("⚠ outside map range")
        else:
            res_b.set_color("#e8622a")
            res_b.set_text(f"p_fuel = {pf_val:.1f} psi   |   p_ox = {po_val:.1f} psi")
            add_point_b(po_val, pf_val)

        fig.canvas.draw_idle()

    def on_clear(_):
        pts_a_data.clear()
        pts_b_data.clear()
        pts_a.set_data([], [])
        pts_b.set_data([], [])
        res_a.set_text("")
        res_b.set_text("")
        fig.canvas.draw_idle()

    btn_a.on_clicked(on_lookup_a)
    btn_b.on_clicked(on_lookup_b)
    btn_clear.on_clicked(on_clear)

    plt.show()


# ── entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    min_OF     = 1.5
    max_OF     = 5.5
    max_pc     = 250      # psi
    resolution = 120
    cstar_eff  = 0.95

    if HAS_PRESSURE_MAP:
        yaml_path = find_yaml()
        if yaml_path:
            yaml = YAML()
            yaml.preserve_quotes = True
            with open(yaml_path) as f:
                data = yaml.load(f)
            cstar_eff   = data["cstar_eff"]
            throat_area = data["throat_area"] / (m_to_in ** 2)
            fuel_CdA    = data["fuel_CdA"]    / (m_to_in ** 2)
            ox_CdA      = data["gox_CdA"]     / (m_to_in ** 2)
            fuel_choice = data["fuel_choice"]
            ox_choice   = data["ox_choice"]
        else:
            print("YAML not found – using defaults.")
            throat_area = 0.5 / (m_to_in ** 2)
            fuel_CdA    = 0.1 / (m_to_in ** 2)
            ox_CdA      = 0.1 / (m_to_in ** 2)
            fuel_choice = "ethanol"
            ox_choice   = "gox"

        pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo = pressure_map(
            minOF=min_OF, maxOF=max_OF, maxPc=max_pc, resolution=resolution,
            cstar_eff=cstar_eff, throat_area=throat_area,
            fuel_CdA=fuel_CdA, ox_CdA=ox_CdA,
            fuel_choice=fuel_choice, ox_choice=ox_choice,
        )

    else:
        # ── synthetic demo data so the UI still works ──────────────────────
        pc_scale = np.linspace(50, max_pc, resolution)
        OF_scale = np.linspace(min_OF, max_OF, resolution)
        PC, OF   = np.meshgrid(pc_scale, OF_scale, indexing="ij")

        fu_p_map = PC * 1.6 + OF * 10   # fake feed pressures
        ox_p_map = PC * 1.4 + OF * 8

        pf_ax = np.linspace(fu_p_map.min(), fu_p_map.max(), resolution)
        po_ax = np.linspace(ox_p_map.min(), ox_p_map.max(), resolution)
        PF, PO = np.meshgrid(pf_ax, po_ax, indexing="ij")

        # crude inverse: solve linear system analytically for demo
        # p_fu = 1.6*Pc + 10*OF  →  Pc ≈ (p_fu - 10*OF) / 1.6
        # p_ox = 1.4*Pc +  8*OF
        # Solve: [[1.6, 10],[1.4, 8]] * [Pc, OF] = [pf, po]
        A = np.array([[1.6, 10.0], [1.4, 8.0]])
        Ainv = np.linalg.inv(A)
        sol = Ainv @ np.array([PF.ravel(), PO.ravel()])
        PC_pfpo = sol[0].reshape(PF.shape)
        OF_pfpo = sol[1].reshape(PF.shape)

    launch_app(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo, max_pc=max_pc)