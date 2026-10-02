"""One real DR6 epoch per defect type: the raw epoch stamp, then the image
and weights ngmix receives (BLEND_HANDLING = noisefill, the default).

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
from shapepipe.modules.ngmix_package.ngmix import (
    defect_mask,
    interpolated_defects,
    prepare_ngmix_weights,
)

COLS = ["raw stamp", "image ngmix sees", "weight"]


def main(out_path, rows):
    plt.rcParams.update({"font.size": 10})
    fig, axs = plt.subplots(len(rows), 3, figsize=(7.4, 2.35 * len(rows)),
                            dpi=150,
                            gridspec_kw=dict(hspace=0.06, wspace=0.05))
    for i, (kind, path) in enumerate(rows):
        ep = dict(np.load(path, allow_pickle=True))
        rms = float(np.median(ep["bkg_rms"]))
        stretch = dict(vmax=120.0 * rms, soft=3.0 * rms)
        defect = defect_mask(ep["weight"], ep["flag"], ep["bkg_rms"])
        interp = interpolated_defects(defect, ep["neighbour"], "noisefill")
        img, w, _ = prepare_ngmix_weights(
            ep["gal"], ep["weight"], ep["flag"], np.random.RandomState(1),
            bkg_rms=ep["bkg_rms"], neighbour=ep["neighbour"],
        )
        show_image(axs[i, 0], ep["gal"], **stretch)
        outline(axs[i, 0], ep["neighbour"], BLUE, lw=0.9)
        outline(axs[i, 0], defect, INK, lw=0.9)
        show_image(axs[i, 1], img, **stretch)
        outline(axs[i, 1], interp, AQUA, lw=0.9)
        outline(axs[i, 1], defect & ~interp, ORANGE, lw=0.9)
        show_weight(axs[i, 2], w, interp, defect | ep["neighbour"])
        axs[i, 0].set_ylabel(
            f"{kind}\n{ep['label']}\n{ep['source']}", color=INK, fontsize=8,
            linespacing=1.15,
        )
    for j, t in enumerate(COLS):
        axs[0, j].set_title(t, color=INK, fontsize=10)
    for ax in axs.flat:
        clean(ax)

    handles = [
        Line2D([], [], color=INK, lw=1.6, label="defect"),
        Line2D([], [], color=BLUE, lw=1.6,
               label="neighbour (noise-filled, zero weight)"),
        Line2D([], [], color=AQUA, lw=1.6, label="defect, interpolated"),
        Line2D([], [], color=ORANGE, lw=1.6, label="defect, noise-filled"),
    ] + WEIGHT_LEGEND
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               fontsize=8.5, handlelength=1.4,
               bbox_to_anchor=(0.5, axs[-1, 0].get_position().y0 - 0.005))
    fig.savefig(out_path, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]])
