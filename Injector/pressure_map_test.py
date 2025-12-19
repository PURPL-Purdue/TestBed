import numpy as np
import matplotlib.pyplot as plt
from pressure_map import pressure_map

# -------------------------
# Constants / conversions
# -------------------------
ft_to_m = 0.3048
m_to_in = 39.3701

# -------------------------
# Helpers
# -------------------------
def finite_minmax(a):
    a = np.asarray(a, dtype=float)
    m = np.isfinite(a)
    if not np.any(m):
        raise ValueError("Array has no finite values.")
    return float(np.min(a[m])), float(np.max(a[m]))

def plot_pc_heatmap_with_of_contours(fu_p_map, ox_p_map, PC_pfpo, OF_pfpo,
                                     of_levels=None, pc_vmin=None, pc_vmax=None):
    """
    Plots a separate figure in (p_fu, p_ox) space:
      - heatmap: Pc (psi) from PC_pfpo
      - contours: OF from OF_pfpo

    Notes:
      - p_fu and p_ox axes are inferred from min/max of forward pressure maps (Pa)
      - NaNs are left as-is (white / background), no masking.
    """

    # infer pressure axes range (Pa) from the forward maps
    pf_min, pf_max = np.nanmin(fu_p_map), np.nanmax(fu_p_map)
    po_min, po_max = np.nanmin(ox_p_map), np.nanmax(ox_p_map)

    # If your maps can contain inf, protect against it
    if not np.isfinite([pf_min, pf_max, po_min, po_max]).all():
        pf = fu_p_map[np.isfinite(fu_p_map)]
        po = ox_p_map[np.isfinite(ox_p_map)]
        pf_min, pf_max = float(np.min(pf)), float(np.max(pf))
        po_min, po_max = float(np.min(po)), float(np.max(po))

    # These are the axes for the inverse maps (constructed via linspace in your solver)
    p_fu_axis = np.linspace(pf_min, pf_max, PC_pfpo.shape[0])
    p_ox_axis = np.linspace(po_min, po_max, PC_pfpo.shape[1])
    PO, PF = np.meshgrid(p_ox_axis, p_fu_axis)  # x=p_ox, y=p_fu

    # default OF contour levels if not provided
    if of_levels is None:
        of_finite = OF_pfpo[np.isfinite(OF_pfpo)]
        if of_finite.size > 0:
            lo, hi = np.percentile(of_finite, [10, 90])
            of_levels = np.linspace(lo, hi, 6)
        else:
            of_levels = np.linspace(min_OF, max_OF, 5)

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)

    # Heatmap of Pc (psi)
    im = ax.imshow(
        PC_pfpo,
        origin="lower",
        aspect="auto",
        extent=[po_min, po_max, pf_min, pf_max],  # x=p_ox, y=p_fu
        vmin=pc_vmin,
        vmax=pc_vmax,
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Pc (psi)")

    # OF equipotential contours on top
    cs = ax.contour(
        PO, PF, OF_pfpo,
        levels=of_levels,
        colors="k",
        linewidths=1.0
    )
    ax.clabel(cs, fmt=lambda v: f"OF={v:.2f}", inline=True, fontsize=9)

    ax.set_xlabel("p_ox (psi)")
    ax.set_ylabel("p_fu (psi)")
    ax.set_title("Pc heatmap with OF contours in (p_fu, p_ox) space")

    plt.show()

def plot_maps(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo):
    # Forward extents: x=OF, y=pc
    extent_forward = [OF_scale.min(), OF_scale.max(), pc_scale.min(), pc_scale.max()]

    # Inverse extents: x=p_ox, y=p_fu (both in Pa)
    pf_min, pf_max = finite_minmax(fu_p_map)
    po_min, po_max = finite_minmax(ox_p_map)
    extent_inv = [po_min, po_max, pf_min, pf_max]

    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    im0 = axs[0, 0].imshow(np.ma.masked_invalid(fu_p_map), origin="lower",
                           aspect="auto", extent=extent_forward)

    axs[0, 0].set_title("Fuel feed pressure (psi)")
    axs[0, 0].set_xlabel("OF (N/A)")
    axs[0, 0].set_ylabel("Chamber pressure (psi)")
    fig.colorbar(im0, ax=axs[0, 0])

    im1 = axs[0, 1].imshow(np.ma.masked_invalid(ox_p_map), origin="lower",
                           aspect="auto", extent=extent_forward)
    axs[0, 1].set_title("Ox feed pressure (psi)")
    axs[0, 1].set_xlabel("OF (N/A)")
    axs[0, 1].set_ylabel("Chamber pressure (psi)")
    fig.colorbar(im1, ax=axs[0, 1])

    im2 = axs[1, 0].imshow(np.ma.masked_invalid(PC_pfpo), origin="lower",
                           aspect="auto", extent=extent_inv)
    axs[1, 0].set_title("Chamber pressure (psi)")
    axs[1, 0].set_xlabel("Ox feed pressure (psi)")
    axs[1, 0].set_ylabel("Fuel feed pressure (psi)")
    fig.colorbar(im2, ax=axs[1, 0])

    im3 = axs[1, 1].imshow(np.ma.masked_invalid(OF_pfpo), origin="lower",
                           aspect="auto", extent=extent_inv)
    axs[1, 1].set_title("OF Ratio (N/A)")
    axs[1, 1].set_xlabel("Ox feed pressure (psi)")
    axs[1, 1].set_ylabel("Fuel feed pressure (psi)")
    fig.colorbar(im3, ax=axs[1, 1])

    plt.show()

# -------------------------
# Main: build forward + inverse maps
# -------------------------

# -------------------------
# Run demo
# -------------------------
if __name__ == "__main__":
    # Example inputs (edit to your engine)
    min_OF = 0.5
    max_OF = 5
    max_pc = 300 # psi
    resolution = 200
    cstar_eff = 0.7

    throat_area = 1.372 / (m_to_in ** 2)

    fuel_CdA = 0.028891 / (m_to_in ** 2)
    ox_CdA = 0.099651 / (m_to_in ** 2)
    rho_fuel = 810.0 # Jet-A-ish kg/m^3

    pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo = pressure_map(
        minOF=min_OF,
        maxOF=max_OF,
        maxPc=max_pc,
        resolution=resolution,
        cstar_eff=cstar_eff,
        throat_area=throat_area,
        fuel_CdA=fuel_CdA,
        ox_CdA=ox_CdA,
        rho_fuel=rho_fuel)

    plot_pc_heatmap_with_of_contours(
    fu_p_map=fu_p_map,
    ox_p_map=ox_p_map,
    PC_pfpo=PC_pfpo,
    OF_pfpo=OF_pfpo,
    of_levels=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
    pc_vmin=0,
    pc_vmax=max_pc)

    plot_maps(pc_scale, OF_scale, fu_p_map, ox_p_map, PC_pfpo, OF_pfpo)