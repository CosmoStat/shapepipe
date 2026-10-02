"""One real DR6 epoch per defect type: the raw epoch stamp, then the image
and weights ngmix receives under DEFECT_FILL = noise and
DEFECT_FILL = interpolate (BLEND_HANDLING = noisefill throughout).

Each epoch is an npz written by ``dr6_epochs.py cut``. Run from the
shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src:scripts/validation/masked_pixels \
    python scripts/validation/masked_pixels/defect_gallery.py OUT.png \
      "bad column=ep_4202.npz" "3-column cluster=ep_28647.npz" ...
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from figstyle import (AQUA, BLUE, INK, ORANGE, WEIGHT_LEGEND, clean, outline,
                      show_image, show_weight)
from shapepipe.modules.ngmix_package.defect_interpolation import (
    interpolable_defects,
)
from shapepipe.modules.ngmix_package.ngmix import (
    defect_mask,
    prepare_ngmix_weights,
)

COLS = ["raw stamp", "image", "weight", "image", "weight"]


def main(out_path, rows):
    plt.rcParams.update({"font.size": 10})
    fig, axs = plt.subplots(len(rows), 5, figsize=(11.2, 2.2 * len(rows)),
                            dpi=150,
                            gridspec_kw=dict(hspace=0.06, wspace=0.05))
    for i, (kind, path) in enumerate(rows):
        ep = dict(np.load(path, allow_pickle=True))
        rms = float(np.median(ep["bkg_rms"]))
        stretch = dict(vmax=120.0 * rms, soft=3.0 * rms)
        defect = defect_mask(ep["weight"], ep["flag"], ep["bkg_rms"])
        interp = interpolable_defects(defect)
        out = {}
        for fill in ("noise", "interpolate"):
            out[fill] = prepare_ngmix_weights(
                ep["gal"], ep["weight"], ep["flag"],
                np.random.RandomState(1), bkg_rms=ep["bkg_rms"],
                neighbour=ep["neighbour"], defect_fill=fill,
            )[:2]
        show_image(axs[i, 0], ep["gal"], **stretch)
        outline(axs[i, 0], ep["neighbour"], BLUE, lw=0.9)
        outline(axs[i, 0], defect, INK, lw=0.9)
        show_image(axs[i, 1], out["noise"][0], **stretch)
        show_weight(axs[i, 2], out["noise"][1])
        show_image(axs[i, 3], out["interpolate"][0], **stretch)
        outline(axs[i, 3], interp, AQUA, lw=0.9)
        outline(axs[i, 3], defect & ~interp, ORANGE, lw=0.9)
        show_weight(axs[i, 4], out["interpolate"][1], interp)
        axs[i, 0].set_ylabel(
            f"{kind}\n{ep['label']}\n{ep['source']}", color=INK, fontsize=8,
            linespacing=1.15,
        )
    for j, t in enumerate(COLS):
        axs[0, j].set_title(t, color=INK, fontsize=10)
    for ax in axs.flat:
        clean(ax)

    top = axs[0, 1].get_position().y1 + 0.03
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
        Line2D([], [], color=BLUE, lw=1.6,
               label="neighbour (raw stamp; noise-filled, zero weight)"),
        Line2D([], [], color=ORANGE, lw=1.6, label="defect, noise-filled"),
        Line2D([], [], color=AQUA, lw=1.6, label="defect, interpolated"),
    ] + WEIGHT_LEGEND
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
               fontsize=9, handlelength=1.4,
               bbox_to_anchor=(0.5, axs[-1, 0].get_position().y0 - 0.005))
    fig.savefig(out_path, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]])
