# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib", "seaborn"]
# ///
"""Column figure for the resolution x noise breakdown grid.

Reads the JSON emitted by ``run_breakdown_grid.py`` and draws three vertically
stacked panels sharing the x-axis -- multiplicative bias ``m``, metacal response
``R11``, and additive bias ``|c1|`` vs the resolution ratio ``r`` -- with the
three S/N bands as hue. The ``m`` panel is linear with a widened y-range so the
moderate/low-S/N climb stays on-plot while the +/-5e-3 tolerance band is still
drawn. Seaborn is not in the ShapePipe sif, so this runs on the host:

    uv run scripts/python/plot_breakdown_grid.py \\
        --input <run>/results/baseline/breakdown_grid/breakdown_grid.json \\
        --output <figures>/breakdown_grid.png
"""
import argparse
import json

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BANDS = [("high", r"S/N $\approx$ 1000"), ("moderate", r"S/N $\approx$ 50"),
         ("low", r"S/N $\approx$ 15")]


def series(cells, band):
    rows = sorted((c for c in cells if c["s2n_band"] == band),
                  key=lambda c: c["res_ratio"])
    return (np.array([c["res_ratio"] for c in rows]),
            np.array([c["m1"] for c in rows]),
            np.array([c["R11"] for c in rows]),
            np.array([abs(c["c1"]) for c in rows]))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--input", required=True, help="breakdown_grid.json path.")
    p.add_argument("--output", required=True, help="output PNG path.")
    a = p.parse_args()

    d = json.load(open(a.input))
    cells = d["cells"]
    m_tol = d.get("meta", {}).get("m_tol_resolved", 5e-3)

    sns.set_theme(style="ticks", context="talk")
    pal = dict(zip([b for b, _ in BANDS], sns.color_palette("husl", len(BANDS))))

    fig, (ax_m, ax_R, ax_c) = plt.subplots(3, 1, figsize=(7.6, 12.5), sharex=True)

    for band, label in BANDS:
        r, m, R11, ac1 = series(cells, band)
        kw = dict(color=pal[band], marker="o", ms=7, lw=2, label=label)
        ax_m.plot(r, m, **kw)
        ax_R.plot(r, R11, **kw)
        ax_c.plot(r, ac1, **kw)

    # m panel: linear, y-range widened to hold the moderate/low-S/N climb while
    # the +/-tol band stays drawn as a sliver at 0.
    ax_m.axhspan(-m_tol, m_tol, color="0.75", alpha=0.4, zorder=0,
                 label=r"$|m| < 5\times10^{-3}$ tol.")
    ax_m.axhline(0, color="0.5", lw=0.8, zorder=0)
    ax_m.set_ylim(-0.05, 1.5)
    ax_m.set_ylabel("multiplicative bias  $m$")
    ax_m.set_title("Multiplicative bias")

    ax_R.axhline(1.0, color="0.5", lw=0.8, ls="--", zorder=0)
    ax_R.set_ylim(0, 1.15)
    ax_R.set_ylabel("metacal response  $R_{11}$")
    ax_R.set_title("Response")

    ax_c.set_yscale("log")
    ax_c.set_ylabel("additive bias  $|c_1|$")
    ax_c.set_title("Additive bias")
    ax_c.set_xlabel("resolution ratio  $r = \\mathrm{gal\\ hlr}\\,/\\,\\mathrm{psf\\ fwhm}$")

    for ax in (ax_m, ax_R, ax_c):
        ax.margins(x=0.04)

    ax_m.legend(frameon=False, fontsize=12, loc="upper right")
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(a.output, dpi=140, bbox_inches="tight")
    print("wrote", a.output)


if __name__ == "__main__":
    main()
