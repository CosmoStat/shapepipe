# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib", "seaborn"]
# ///
"""Column figure for the resolution x noise breakdown grid (v2).

Reads the JSON emitted by ``run_breakdown_grid.py`` (v2 schema) and draws three
vertically stacked panels sharing the x-axis -- multiplicative bias ``m1`` and
``m2``, and additive bias ``|c1|`` vs the resolution ratio ``r`` -- with the
three S/N bands as hue. **Every point carries its bootstrap error bar**; that is
the whole point of v2. Cells whose response was consistent with zero (the
degeneracy guard fired, ``m`` null) are drawn as open markers at 0 with a hatch
so a reader never mistakes "we refused to divide" for "m = 0". The ``m`` panels
keep the ``+/-5e-3`` tolerance band drawn as a sliver at 0. Seaborn is not in the
ShapePipe sif, so this runs on the host:

    uv run scripts/python/plot_breakdown_grid.py \\
        --input <run>/results/baseline/breakdown_grid/breakdown_grid_v2.json \\
        --output <figures>/breakdown_grid.png
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

BANDS = [("high", r"S/N $\approx$ 1000"), ("moderate", r"S/N $\approx$ 50"),
         ("low", r"S/N $\approx$ 15")]


def series(cells, band, key):
    """(r, value, err, degenerate-mask) for one band and one estimator key.

    ``value``/``err`` are NaN where the cell is degenerate (m null); the mask
    flags those so the plotter can render them as open markers at 0."""
    rows = sorted((c for c in cells if c["s2n_band"] == band),
                  key=lambda c: c["res_ratio"])
    r = np.array([c["res_ratio"] for c in rows])
    deg_key = "c_degenerate" if key.startswith("c") else f"{key}_degenerate"
    deg = np.array([bool(c.get(deg_key, False)) or c[key] is None for c in rows],
                   dtype=bool)
    val = np.array([np.nan if c[key] is None else c[key] for c in rows], dtype=float)
    err = np.array([np.nan if c.get(f"{key}_err") is None else c[f"{key}_err"]
                    for c in rows], dtype=float)
    return r, val, err, deg


def draw_panel(ax, cells, key, pal, transform=lambda v: v):
    """Plot one estimator across all bands with error bars; degenerate cells as
    open hatched markers at 0."""
    for band, _ in BANDS:
        r, val, err, deg = series(cells, band, key)
        good = ~deg
        ax.errorbar(r[good], transform(val[good]),
                    yerr=None if err is None else err[good],
                    color=pal[band], marker="o", ms=7, lw=2, capsize=3,
                    elinewidth=1.4, label=dict(BANDS)[band])
        if deg.any():  # open markers at 0 mark "response ~ 0, refused to divide"
            ax.scatter(r[deg], np.zeros(deg.sum()), s=90, facecolors="none",
                       edgecolors=pal[band], linewidths=1.8, hatch="////",
                       zorder=5)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--input", required=True, help="breakdown_grid JSON path (v2).")
    p.add_argument("--output", required=True, help="output PNG path.")
    a = p.parse_args()

    d = json.load(open(a.input))
    cells = d["cells"]
    m_tol = d.get("meta", {}).get("m_tol_resolved", 5e-3)

    sns.set_theme(style="ticks", context="talk")
    pal = dict(zip([b for b, _ in BANDS], sns.color_palette("husl", len(BANDS))))

    fig, (ax_m1, ax_m2, ax_c) = plt.subplots(3, 1, figsize=(7.6, 12.5), sharex=True)

    for ax, key, title, ylab in (
        (ax_m1, "m1", "Multiplicative bias $m_1$", r"multiplicative bias  $m_1$"),
        (ax_m2, "m2", "Multiplicative bias $m_2$", r"multiplicative bias  $m_2$"),
    ):
        draw_panel(ax, cells, key, pal)
        ax.axhspan(-m_tol, m_tol, color="0.75", alpha=0.4, zorder=0,
                   label=r"$|m| < 5\times10^{-3}$ tol.")
        ax.axhline(0, color="0.5", lw=0.8, zorder=0)
        ax.set_ylabel(ylab)
        ax.set_title(title)

    # |c1| on a symlog axis so the tolerance-scale values and any blow-up both
    # stay legible while error bars that cross zero still render.
    draw_panel(ax_c, cells, "c1", pal, transform=np.abs)
    ax_c.set_yscale("symlog", linthresh=1e-4)
    ax_c.set_ylabel(r"additive bias  $|c_1|$")
    ax_c.set_title("Additive bias")
    # "intrinsic (pre-seeing)" is load-bearing: r is the simulation-input truth
    # knob (r -> 0 is the star limit), NOT the convolved observed size ratio
    # (which is floored at 1) nor metacal's deconvolved T_gal/T_psf.
    ax_c.set_xlabel("intrinsic (pre-seeing) resolution ratio  "
                    "$r = \\mathrm{gal\\ hlr}\\,/\\,\\mathrm{psf\\ fwhm}$")

    for ax in (ax_m1, ax_m2, ax_c):
        ax.margins(x=0.04)

    handles = [Line2D([0], [0], color=pal[b], marker="o", lw=2, label=lab)
               for b, lab in BANDS]
    handles.append(Line2D([0], [0], marker="o", ls="none", mfc="none",
                          mec="0.3", mew=1.8, label="degenerate ($R\\approx0$)"))
    ax_m1.legend(handles=handles, frameon=False, fontsize=12, loc="upper left")
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(a.output, dpi=140, bbox_inches="tight")
    print("wrote", a.output)


if __name__ == "__main__":
    main()
