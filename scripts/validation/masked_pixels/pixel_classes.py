"""Pixel classes of one epoch stamp and what ngmix receives: a synthetic
blend under BLEND_HANDLING = noisefill and uberseg, then a real DR6 blend
under uberseg.

Run from the shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src:scripts/validation/masked_pixels \
    python scripts/validation/masked_pixels/pixel_classes.py REAL.npz OUT.png

REAL.npz is one real epoch written by dr6_epochs.py.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch

from figstyle import (AQUA, BLUE, INK, INK2, ORANGE, PAPER, TARGET,
                      WEIGHT_LEGEND, clean, show_image, show_weight)
from shapepipe.modules.ngmix_package.ngmix import (
    EPOCH_CENTRAL_DEFECT_RADIUS,
    EPOCH_INTERPOLATED_DEFECT_RADIUS,
    defect_mask,
    interpolated_defects,
    prepare_ngmix_weights,
    uberseg_mask,
)


def synthetic():
    """The synthetic epoch: a target, a neighbour 15 px away, a bad column,
    a hot pixel and an off-tile band."""
    N, C = 51, 25
    rng = np.random.RandomState(4)
    yy, xx = np.indices((N, N))

    def blob(r0, c0, sigma, flux, q=1.0, theta=0.0):
        dy, dx = yy - r0, xx - c0
        ct, st = np.cos(theta), np.sin(theta)
        u, v = ct * dx + st * dy, -st * dx + ct * dy
        r = np.sqrt(u ** 2 + (v / q) ** 2)
        return flux * np.exp(-r / sigma) / (2 * np.pi * sigma ** 2 * q)

    rms = 1.0
    target = blob(C, C, 2.2, 3000, q=0.7, theta=0.5)
    nbr = blob(C - 10, C - 11, 2.0, 4000, q=0.8, theta=-0.3)
    gal = target + nbr + rng.normal(0, rms, (N, N))
    # Segmentation: each pixel above 2 sigma joins the brighter profile.
    above = (target + nbr) > 2 * rms
    seg = np.where(above, np.where(target >= nbr, 1, 2), 0).astype(np.int32)
    flag = np.zeros((N, N), np.int32)
    flag[:, C + 12] = 1               # bad column, 12 px right of centre
    flag[C + 9, C - 9] = 1            # hot pixel, 12.7 px from centre
    off_tile = np.zeros((N, N), bool)
    off_tile[N - 5:, :] = True        # off-tile band at the stamp's edge
    gal[:, C + 12] = 400.0
    gal[C + 9, C - 9] = 400.0
    return dict(gal=gal, weight=np.where(off_tile, 0.0, 1.0), flag=flag,
                bkg_rms=np.full((N, N), rms), seg=seg, object_number=1,
                neighbour=seg == 2, off_tile=off_tile)


def run(ep, blend_handling):
    kw = {}
    if blend_handling == "uberseg":
        kw.update(seg=ep["seg"], object_number=int(ep["object_number"]))
    img, w, _ = prepare_ngmix_weights(
        ep["gal"], ep["weight"], ep["flag"], np.random.RandomState(1),
        bkg_rms=ep["bkg_rms"], neighbour=ep["neighbour"],
        blend_handling=blend_handling, **kw,
    )
    defect = defect_mask(ep["weight"], ep["flag"], ep["bkg_rms"])
    interp = interpolated_defects(defect, ep["neighbour"], blend_handling)
    neighbour = (ep["neighbour"] if blend_handling == "noisefill" else
                 uberseg_mask(ep["seg"], int(ep["object_number"])))
    return img, w, interp, defect | neighbour


CLASS_CMAP = ListedColormap([PAPER, BLUE, ORANGE, AQUA, TARGET])


def classes(ep):
    cls = np.zeros(ep["gal"].shape, int)
    target = ep["seg"] == int(ep["object_number"])
    cls[target] = 4
    cls[ep["neighbour"]] = 1
    defect = defect_mask(ep["weight"], ep["flag"], ep["bkg_rms"])
    cls[defect] = 2
    cls[ep["off_tile"]] = 3
    return cls


def class_panel(ax, ep, radii=True):
    ax.imshow(classes(ep), origin="lower", cmap=CLASS_CMAP, vmin=-0.5,
              vmax=4.5, interpolation="nearest")
    if radii:
        n = ep["gal"].shape[0]
        c = (n - 1) / 2
        for r, ls in ((EPOCH_CENTRAL_DEFECT_RADIUS, "--"),
                      (EPOCH_INTERPOLATED_DEFECT_RADIUS, ":")):
            ax.add_patch(Circle((c, c), r, fill=False, ec=INK, ls=ls,
                                lw=1.0))


def main(real_path, out_path):
    syn = synthetic()
    real = dict(np.load(real_path, allow_pickle=True))
    # The synthetic stamp has unit noise; scale the real one's stretch by
    # its background RMS so both read alike.
    rms = float(np.median(real["bkg_rms"]))
    vmax_real, soft_real = 120.0 * rms, 3.0 * rms

    plt.rcParams.update({"font.size": 10})
    fig, axs = plt.subplots(
        2, 5, figsize=(13.6, 6.1), dpi=150,
        gridspec_kw=dict(hspace=0.08, wspace=0.05,
                         width_ratios=[1, 1, 1, 1, 1]),
    )
    # Open a gap between the synthetic and the real columns.
    for ax in axs[:, 3:].flat:
        p = ax.get_position()
        ax.set_position([p.x0 + 0.025, p.y0, p.width, p.height])

    # (a) synthetic epoch and its pixel classes
    show_image(axs[0, 0], syn["gal"])
    class_panel(axs[1, 0], syn)
    axs[0, 0].set_title("(a) one epoch's stamp", color=INK, fontsize=10)
    axs[0, 0].set_ylabel("image ngmix sees", color=INK, fontsize=11)
    axs[1, 0].set_ylabel("pixel class / weight", color=INK, fontsize=11)
    c = 25
    rn = EPOCH_CENTRAL_DEFECT_RADIUS
    axs[1, 0].text(c, c + rn + 0.8, f"{rn:g} px", ha="center", va="bottom",
                   fontsize=8, color=INK)
    axs[1, 0].text(c, c - 6.3, "7 px", ha="center", va="bottom",
                   fontsize=8, color=INK)

    for j, (blend, title) in enumerate(
        (("noisefill", "(b) noisefill (default)"),
         ("uberseg", "(c) uberseg")), start=1):
        img, w, interp, zeroed = run(syn, blend)
        show_image(axs[0, j], img)
        show_weight(axs[1, j], w, interp, zeroed)
        axs[0, j].set_title(title, color=INK, fontsize=10)

    # (d) the real DR6 epoch and its classes, (e) under uberseg
    show_image(axs[0, 3], real["gal"], vmax=vmax_real, soft=soft_real)
    class_panel(axs[1, 3], real)
    axs[0, 3].set_title(f"(d) {real['label']}\n{real['source']}",
                        color=INK, fontsize=10)
    img, w, interp, zeroed = run(real, "uberseg")
    show_image(axs[0, 4], img, vmax=vmax_real, soft=soft_real)
    show_weight(axs[1, 4], w, interp, zeroed)
    axs[0, 4].set_title("(e) uberseg", color=INK, fontsize=10)
    for ax in axs.flat:
        clean(ax)

    y = axs[0, 0].get_position().y1 + 0.075
    for (a, b), t in (((0, 2), "synthetic blend"),
                      ((3, 4), "real DR6 blend")):
        x0 = axs[0, a].get_position().x0
        x1 = axs[0, b].get_position().x1
        fig.text((x0 + x1) / 2, y, t, ha="center", va="bottom", color=INK,
                 fontsize=11)
        fig.add_artist(Line2D([x0 + 0.004, x1 - 0.004], [y - 0.006] * 2,
                              color=INK2, lw=0.8))

    handles = [
        Patch(color=TARGET, label="target footprint"),
        Patch(color=BLUE, label="neighbour (−1e30 in tile)"),
        Patch(color=ORANGE, label="defect (flag / weight / RMS)"),
        Patch(color=AQUA, label="off-tile (a defect)"),
        Line2D([], [], color=INK, ls="--", lw=1, label=f"{EPOCH_CENTRAL_DEFECT_RADIUS:g} px veto"),
        Line2D([], [], color=INK, ls=":", lw=1,
               label="7 px veto (interpolated)"),
    ] + WEIGHT_LEGEND
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False,
               fontsize=9, handlelength=1.4,
               bbox_to_anchor=(0.5, axs[1, 0].get_position().y0 - 0.01))
    fig.savefig(out_path, bbox_inches="tight")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
