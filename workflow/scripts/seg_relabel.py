"""Map a SExtractor SEGMENTATION map into the tile catalogue's NUMBER space.

WHY THIS EXISTS. UberSeg identifies the central object of a stamp BY LABEL: it
keeps the pixels whose nearest segmentation footprint carries the object's own
``NUMBER`` and zeros every other pixel's weight
(``ngmix_package/ngmix.py::uberseg_weight``, called with ``object_number =
obj_id``, the catalogue's ``NUMBER``). It never reads the label under the
object's position. So the seg stamp ngmix overlays must be labelled in the
SAME numbering as the catalogue it measures, or every mask is inverted.

Under ``tile_detection: unions_catalogue`` the two numberings are different by
construction: the sample is Steven Gwyn's per-tile catalogue, renumbered 1..N
by ``read_ext_sexcat``, while the segmentation map comes from a separate
SExtractor run whose labels are ITS OWN detection numbers. This script is the
bridge, and it is the whole of the guarantee:

  * each catalogue object claims the footprint its position falls in -- the one
    lookup that is meaningful, since the two runs share the tile's pixel grid;
  * that footprint is relabelled with the object's catalogue ``NUMBER``;
  * every footprint no catalogue object claims is relabelled ``NEIGHBOUR_LABEL``
    (negative, so it can never collide with a ``NUMBER``): UberSeg's partition
    only ever asks "self or not self", so neighbours need no identity;
  * an object whose position falls on sky, or on a footprint another object
    claimed first, is given a small disc of its own ``NUMBER`` at its position.

THE FALLBACK IS NOT COSMETIC. Without a footprint of its own an object has no
"self" in the Voronoi partition: ``uberseg_weight`` would hand it an all-zero
weight, and ``Ngmix._check_central_seg_label`` raises on a stamp that does not
contain the object's label at all. The disc gives such an object a defined,
conservative core -- small, centred, and (in the shared-footprint case) taking
its pixels from the neighbour that claimed the blend, which is the right way
round: the claimant keeps the bulk, the unclaimed object keeps a centre.

Claiming is in ``NUMBER`` order and first come first served, so the mapping is
deterministic and reproducible from the same two inputs.

Stdlib + numpy/astropy (both in the container). Writes ``int32``, preserving
the check image's header so the seg map stays on the tile's WCS -- vignetmaker
run 3 cuts the stamps by position, not by row.
"""

import argparse
import sys

import numpy as np
from astropy.io import fits

# The label every unclaimed footprint carries. Negative so that no catalogue
# NUMBER (1..N) can ever collide with it; UberSeg tests `seg != object_number`,
# so one shared label for all neighbours is enough.
NEIGHBOUR_LABEL = -1

# Radius, in pixels, of the disc given to an object with no footprint of its
# own. ~0.56" at the CFIS pixel scale: a plausible minimum galaxy core, small
# enough not to take a blend away from the object that owns its footprint.
FALLBACK_RADIUS = 3


def relabel(seg, number, x_image, y_image, fallback_radius=FALLBACK_RADIUS):
    """Relabel ``seg`` into the numbering of a catalogue.

    Parameters
    ----------
    seg : numpy.ndarray
        Segmentation map, 0 for sky and one positive label per detection.
    number : numpy.ndarray
        The catalogue's ``NUMBER`` column.
    x_image, y_image : numpy.ndarray
        The catalogue's ``XWIN_IMAGE``/``YWIN_IMAGE`` columns, FITS 1-indexed
        pixel coordinates on the same grid as ``seg``.
    fallback_radius : int, optional
        Radius of the disc painted for an object with no footprint of its own.

    Returns
    -------
    (numpy.ndarray, dict)
        The relabelled map (``int32``) and a count of each outcome:
        the four claim outcomes, which partition the catalogue --
        ``matched``, ``unclaimed`` (fell on sky), ``shared`` (fell on a
        footprint already claimed) and ``off_image`` -- plus
        ``shared_pixel``, a separate axis counting objects that round to the
        same pixel as a lower-numbered one.
    """
    number = np.asarray(number)
    # FITS pixel centres are 1-based; round to the pixel the centroid is in.
    col = np.rint(np.asarray(x_image, dtype=float)).astype(np.int64) - 1
    row = np.rint(np.asarray(y_image, dtype=float)).astype(np.int64) - 1
    n_row, n_col = seg.shape
    inside = (col >= 0) & (col < n_col) & (row >= 0) & (row < n_row)

    # Claim, in NUMBER order, the footprint each object's centre falls in.
    owner = {}
    counts = dict(matched=0, unclaimed=0, shared=0, off_image=0,
                  shared_pixel=0)
    fallback = []
    for i in np.argsort(number, kind="stable"):
        if not inside[i]:
            counts["off_image"] += 1
            continue
        label = int(seg[row[i], col[i]])
        if label == 0:
            counts["unclaimed"] += 1
            fallback.append(i)
        elif label in owner:
            counts["shared"] += 1
            fallback.append(i)
        else:
            owner[label] = int(number[i])
            counts["matched"] += 1

    # Every footprint becomes a neighbour unless an object claimed it. A lookup
    # table over the labels present is O(pixels) once, rather than one pass per
    # object.
    labels = np.unique(seg)
    table = np.full(int(labels.max()) + 1, NEIGHBOUR_LABEL, dtype=np.int32)
    table[0] = 0
    for label, num in owner.items():
        table[label] = num
    out = table[np.clip(seg, 0, None)]

    # The fallback discs go on before the centre pass below, so an object that
    # shares a claimed footprint takes its centre back from the claimant.
    if fallback:
        offsets = np.argwhere(
            np.add.outer(
                np.arange(-fallback_radius, fallback_radius + 1) ** 2,
                np.arange(-fallback_radius, fallback_radius + 1) ** 2,
            ) <= fallback_radius ** 2
        ) - fallback_radius
        for i in fallback:
            rows = np.clip(row[i] + offsets[:, 0], 0, n_row - 1)
            cols = np.clip(col[i] + offsets[:, 1], 0, n_col - 1)
            out[rows, cols] = int(number[i])

    # EVERY OBJECT ENDS WITH ITS CENTRE PIXEL, and this last pass is what makes
    # that true rather than usually true. A disc painted for one object can lie
    # over another's pixels -- over a small footprint, or over an earlier disc
    # -- and an object with no pixel of its own is not a soft failure: UberSeg
    # would hand it an all-zero weight and Ngmix._check_central_seg_label
    # raises on a stamp that does not contain the label. Re-stamping the centre
    # costs the overlapping object one pixel and nothing else.
    #
    # Two catalogue objects rounding to the SAME pixel is the one contest a
    # pixel cannot settle twice: the lower NUMBER keeps it and the other is
    # counted as ``shared_pixel``. The loser still holds the rest of its disc,
    # so it keeps a self of its own; the count is there because two objects
    # that close are worth seeing in the stage log.
    claimed_pixel = {}
    for i in np.argsort(number, kind="stable"):
        if not inside[i]:
            continue
        pixel = (int(row[i]), int(col[i]))
        if pixel in claimed_pixel:
            counts["shared_pixel"] += 1
            continue
        claimed_pixel[pixel] = True
        out[row[i], col[i]] = int(number[i])

    return out, counts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sexcat", required=True,
                        help="the tile catalogue (FITS-LDAC) being measured")
    parser.add_argument("--segmentation", required=True,
                        help="the SExtractor SEGMENTATION check image")
    parser.add_argument("--output", required=True,
                        help="where to write the relabelled map")
    parser.add_argument("--fallback-radius", type=int, default=FALLBACK_RADIUS)
    args = parser.parse_args(argv)

    with fits.open(args.sexcat) as hdus:
        cat = hdus["LDAC_OBJECTS"].data
        number = np.array(cat["NUMBER"])
        x_image = np.array(cat["XWIN_IMAGE"])
        y_image = np.array(cat["YWIN_IMAGE"])
    with fits.open(args.segmentation) as hdus:
        seg = hdus[0].data
        header = hdus[0].header

    out, counts = relabel(seg, number, x_image, y_image, args.fallback_radius)
    fits.PrimaryHDU(data=out, header=header).writeto(args.output,
                                                     overwrite=True)
    print(f"seg_relabel: {len(number)} objects, " + ", ".join(
        f"{k}={v}" for k, v in counts.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
