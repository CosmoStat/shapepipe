"""One defect per row: the raw epoch stamp, then the image and weights ngmix
receives under DEFECT_FILL = noise and DEFECT_FILL = interpolate.

Run from the shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src:scripts/validation/masked_pixels \
    python scripts/validation/masked_pixels/defect_gallery.py OUT.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from figstyle import (AQUA, INK, ORANGE, WEIGHT_LEGEND, clean, outline,
                      show_image, show_weight)
from shapepipe.modules.ngmix_package.defect_interpolation import (
    interpolable_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    OFF_TILE_FLAG,
    defect_mask,
    prepare_ngmix_weights,
)

N, C = 51, 25
yy, xx = np.indices((N, N))
RMS = 1.0


def blob(r0, c0, sigma, flux, q=1.0, theta=0.0):
    dy, dx = yy - r0, xx - c0
    ct, st = np.cos(theta), np.sin(theta)
    u, v = ct * dx + st * dy, -st * dx + ct * dy
    r = np.sqrt(u ** 2 + (v / q) ** 2)
    return flux * np.exp(-r / sigma) / (2 * np.pi * sigma ** 2 * q)


def case(kind):
    """Raw stamp, weight and flag for one defect type."""
    rng = np.random.RandomState(7)
    gal = blob(C, C, 2.6, 4000, q=0.7, theta=0.5) + rng.normal(0, RMS, (N, N))
    weight = np.ones((N, N))
    flag = np.zeros((N, N), np.int32)
    if kind == "column":
        flag[:, C + 5] = 1
        gal[:, C + 5] = 60.0
    elif kind == "bleed":
        flag[C - 2:, C - 7:C - 4] = 1
        gal[C - 2:, C - 7:C - 4] = 400.0
    elif kind == "cosmic":
        track = [(C + 4 + i, C + 3 + i // 2) for i in range(9)]
        for r, c in track:
            flag[r, c] = 1
            gal[r, c] = 300.0
        flag[C - 9, C - 4] = 1
        gal[C - 9, C - 4] = 500.0
    elif kind == "edge":
        flag[:6, :] = OFF_TILE_FLAG
    elif kind == "hole":
        weight[C - 3:C + 3, C + 4:C + 10] = 0.0
        gal[C - 3:C + 3, C + 4:C + 10] = 0.0
    return gal, weight, flag


ROWS = [
    ("column", "bad column"),
    ("bleed", "3 px bleed"),
    ("cosmic", "cosmic ray +\nsingle hot pixel"),
    ("edge", "edge band\n(off-tile, 6 rows)"),
    ("hole", "6×6 hole\n(too wide)"),
]
COLS = ["raw stamp", "image", "weight", "image", "weight"]

plt.rcParams.update({"font.size": 10})
fig, axs = plt.subplots(len(ROWS), 5, figsize=(10.2, 2.08 * len(ROWS)),
                        dpi=150, gridspec_kw=dict(hspace=0.06, wspace=0.05))
for i, (kind, label) in enumerate(ROWS):
    gal, weight, flag = case(kind)
    defect = defect_mask(weight, flag)
    interp = interpolable_defects(defect)
    noisefilled = defect & ~interp
    bkg_rms = np.full((N, N), RMS)
    out = {}
    for fill in ("noise", "interpolate"):
        img, w, _ = prepare_ngmix_weights(
            gal, weight, flag, np.random.RandomState(1), bkg_rms=bkg_rms,
            defect_fill=fill,
        )
        out[fill] = img, w
    show_image(axs[i, 0], gal)
    outline(axs[i, 0], defect, INK, lw=0.9)
    show_image(axs[i, 1], out["noise"][0])
    show_weight(axs[i, 2], out["noise"][1])
    show_image(axs[i, 3], out["interpolate"][0])
    outline(axs[i, 3], interp, AQUA, lw=0.9)
    outline(axs[i, 3], noisefilled, ORANGE, lw=0.9)
    show_weight(axs[i, 4], out["interpolate"][1], interp)
    axs[i, 0].set_ylabel(label, color=INK, fontsize=10)
for j, t in enumerate(COLS):
    axs[0, j].set_title(t, color=INK, fontsize=10)
for ax in axs.flat:
    clean(ax)

top = axs[0, 1].get_position().y1 + 0.035
for (a, b), t in (((1, 2), "DEFECT_FILL = noise"),
                  ((3, 4), "DEFECT_FILL = interpolate")):
    x0 = axs[0, a].get_position().x0
    x1 = axs[0, b].get_position().x1
    fig.text((x0 + x1) / 2, top, t, ha="center", va="bottom", color=INK,
             fontsize=11, family="monospace")
    fig.add_artist(Line2D([x0 + 0.005, x1 - 0.005], [top - 0.004] * 2,
                          color=INK, lw=0.8))

handles = [
    Line2D([], [], color=INK, lw=1.6, label="defect (raw stamp)"),
    Line2D([], [], color=ORANGE, lw=1.6, label="defect, noise-filled"),
    Line2D([], [], color=AQUA, lw=1.6, label="defect, interpolated"),
] + WEIGHT_LEGEND
fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
           fontsize=9, handlelength=1.4,
           bbox_to_anchor=(0.5, axs[-1, 0].get_position().y0 - 0.005))
fig.savefig(sys.argv[1], bbox_inches="tight")
