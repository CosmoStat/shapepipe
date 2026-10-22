"""Shared look of the masked-pixel figures (image stretch, weight classes)."""
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from shapepipe.modules.ngmix_package.defect_interpolation import fourfold

INK, INK2 = "#0b0b0b", "#52514e"
PAPER, DARK = "#f4f3ef", "#3a3a38"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#8e5bd0"
TARGET = "#c3c2b7"


def stretch(a, soft=3.0):
    return np.arcsinh(np.asarray(a, dtype=float) / soft)


def show_image(ax, img, vmax=120.0, soft=3.0):
    ax.imshow(stretch(img, soft), origin="lower", cmap="gray_r",
              vmin=stretch(-3 * soft / 3.0, soft), vmax=stretch(vmax, soft),
              interpolation="nearest")


# Weight panel: 0 weighted, 1 zero weight, 2 zero weight on a rotated copy
# of an interpolated pixel (its light is kept).
WEIGHT_CMAP = ListedColormap([PAPER, DARK, VIOLET])


def weight_classes(w, interpolated, zeroed=None):
    """0 = weighted, 1 = zero weight, 2 = zero weight only as a quarter-turn
    copy of an interpolated pixel. ``zeroed`` marks the pixels that have
    zero weight for another reason (defects, neighbours); they stay 1."""
    cls = np.where(w > 0, 0, 1)
    copies = fourfold(interpolated) & ~interpolated
    if zeroed is not None:
        copies &= ~zeroed
    cls[(w == 0) & copies] = 2
    return cls


def show_weight(ax, w, interpolated=None, zeroed=None):
    if interpolated is None:
        interpolated = np.zeros(w.shape, bool)
    ax.imshow(weight_classes(w, interpolated, zeroed), origin="lower",
              cmap=WEIGHT_CMAP, vmin=-0.5, vmax=2.5, interpolation="nearest")


WEIGHT_LEGEND = [
    Patch(color=PAPER, ec=INK2, lw=0.5, label="weighted"),
    Patch(color=DARK, label="zero weight"),
    Patch(color=VIOLET, label="zero weight, light kept: rotated copy of an "
          "interpolated pixel"),
]


def outline(ax, mask, color, lw=1.0):
    """Draw the pixel-edge outline of a boolean mask."""
    mask = np.asarray(mask, bool)
    ny, nx = mask.shape
    pad = np.pad(mask, 1)
    segs = []
    for r in range(ny):
        for c in range(nx):
            if not mask[r, c]:
                continue
            if not pad[r, c + 1]:
                segs.append(((c - .5, c + .5), (r - .5, r - .5)))
            if not pad[r + 2, c + 1]:
                segs.append(((c - .5, c + .5), (r + .5, r + .5)))
            if not pad[r + 1, c]:
                segs.append(((c - .5, c - .5), (r - .5, r + .5)))
            if not pad[r + 1, c + 2]:
                segs.append(((c + .5, c + .5), (r - .5, r + .5)))
    for xs, ys in segs:
        ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="butt")


def clean(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(INK2)
        s.set_linewidth(0.6)
