"""READ EXTERNAL SEXCAT.

Convert an ASCII SExtractor-format catalogue to FITS-LDAC for use in
ShapePipe tile processing, measuring the windowed positions the catalogue
lacks on the tile image.

:Author: Martin Kilbinger

"""

import numpy as np
from astropy.io import ascii as asc
from astropy.io import fits
from astropy.wcs import WCS

from shapepipe.modules.read_ext_sexcat_package import windowed_position as wp


def _build_ldac_imhead(img_header):
    """Build LDAC_IMHEAD extension from a FITS image header.

    The LDAC_IMHEAD binary table has one row and one column containing
    the image header cards as an array of 80-character strings.  The
    HISTORY cards are the ones used by ``make_post_process`` to identify
    the single exposures that contributed to the tile.

    Parameters
    ----------
    img_header : astropy.io.fits.Header
        Header of the tile image FITS file

    Returns
    -------
    astropy.io.fits.BinTableHDU
        LDAC_IMHEAD extension

    """
    cards = [str(card).ljust(80)[:80] for card in img_header.cards]
    n_cards = len(cards)
    header_str = "".join(cards)

    arr = np.array([header_str.encode()], dtype=f"S{len(header_str)}")
    col = fits.Column(
        name="Field Header Card",
        format=f"{len(header_str)}A",
        array=arr,
    )
    hdu = fits.BinTableHDU.from_columns([col])
    hdu.header["TDIM1"] = f"(80,{n_cards})"
    hdu.name = "LDAC_IMHEAD"
    return hdu


# SExtractor's VIGNET value for pixels that are not the object's: off the
# image and on a neighbour's segmentation footprint. ngmix flags these pixels
# (tile VIGNET == -1e30) and noise-fills them at zero weight.
BIG = -1e30

# The label of a footprint no catalogue object claims in the relabelled
# segmentation map. Negative, so it never collides with a NUMBER; VIGNET
# marking only asks "self or not self", so neighbours need no identity.
NEIGHBOUR_LABEL = -1

# Radius, in pixels, of the disc of its own NUMBER painted for an object with
# no footprint of its own (on sky, or on a footprint claimed by another).
FALLBACK_RADIUS = 3


def _centre_pixels(x_pos, y_pos):
    """0-based (column, row) of the pixel holding each 1-based position."""
    col = np.rint(np.asarray(x_pos, dtype=float)).astype(np.int64) - 1
    row = np.rint(np.asarray(y_pos, dtype=float)).astype(np.int64) - 1
    return col, row


def relabel_seg(seg, number, x_image, y_image, fallback_radius=FALLBACK_RADIUS):
    """Relabel a segmentation map into a catalogue's NUMBER.

    Each object claims, in NUMBER order and first come first served, the
    footprint its centre pixel falls in; that footprint becomes its NUMBER.
    Footprints nobody claims become ``NEIGHBOUR_LABEL``; sky stays 0. An
    object on sky or on an already-claimed footprint gets a disc of radius
    ``fallback_radius`` of its NUMBER, and every object on the image ends
    holding its own centre pixel (the lower NUMBER wins a shared pixel).

    Parameters
    ----------
    seg : numpy.ndarray
        Segmentation map, 0 for sky and one positive label per footprint
    number : array_like
        Catalogue ``NUMBER``
    x_image, y_image : array_like
        1-based pixel positions on the grid of ``seg``
    fallback_radius : int, optional
        Radius of the disc painted for an object without a footprint

    Returns
    -------
    numpy.ndarray
        Relabelled map, ``int32``
    dict
        Counts of the claim outcomes ``matched``, ``unclaimed`` (on sky),
        ``shared`` and ``off_image``, which partition the catalogue, plus
        ``shared_pixel`` (objects rounding to a lower NUMBER's pixel)

    @sc [decision:detection.catalogue_neighbour_marking]
    """
    number = np.asarray(number)
    col, row = _centre_pixels(x_image, y_image)
    n_row, n_col = seg.shape
    inside = (col >= 0) & (col < n_col) & (row >= 0) & (row < n_row)
    order = np.argsort(number, kind="stable")

    counts = dict(matched=0, unclaimed=0, shared=0, off_image=0,
                  shared_pixel=0)
    table = np.full(max(int(seg.max()), 0) + 1, NEIGHBOUR_LABEL, np.int32)
    table[0] = 0
    claimed, fallback = set(), []
    for i in order:
        if not inside[i]:
            counts["off_image"] += 1
            continue
        label = int(seg[row[i], col[i]])
        if label == 0 or label in claimed:
            counts["unclaimed" if label == 0 else "shared"] += 1
            fallback.append(i)
        else:
            claimed.add(label)
            table[label] = number[i]
            counts["matched"] += 1
    out = table[np.clip(seg, 0, None)]

    r = np.arange(-fallback_radius, fallback_radius + 1)
    disc = np.argwhere(np.add.outer(r**2, r**2) <= fallback_radius**2)
    disc -= fallback_radius
    for i in fallback:
        out[np.clip(row[i] + disc[:, 0], 0, n_row - 1),
            np.clip(col[i] + disc[:, 1], 0, n_col - 1)] = number[i]

    # Centres last, so no disc can take an object's centre away.
    taken = set()
    for i in order[inside[order]]:
        pixel = (row[i], col[i])
        if pixel in taken:
            counts["shared_pixel"] += 1
            continue
        taken.add(pixel)
        out[pixel] = number[i]
    return out, counts


def _extract_vignets(image_data, x_pos, y_pos, stamp_size, seg=None,
                     number=None):
    """Extract postage stamps from a tile image array.

    For each object position, a ``stamp_size x stamp_size`` cutout is
    extracted, centred on the pixel holding the position. Pixels off the
    image are set to ``BIG``, as SExtractor does. With a segmentation map in
    the catalogue's numbering (:func:`relabel_seg`), pixels on any footprint
    other than the object's own are set to ``BIG`` too, which is how
    SExtractor's VIGNET marks neighbours.

    Parameters
    ----------
    image_data : numpy.ndarray
        2-D tile image array, shape ``(ny, nx)``
    x_pos : array_like
        X pixel positions, 1-based (SExtractor convention)
    y_pos : array_like
        Y pixel positions, 1-based (SExtractor convention)
    stamp_size : int
        Side length of the square postage stamp (should be odd)
    seg : numpy.ndarray, optional
        Segmentation map on the image grid, labelled with ``number``
    number : array_like, optional
        Catalogue ``NUMBER``, required with ``seg``

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n_obj, stamp_size, stamp_size)``, dtype float32

    @sc [decision:detection.catalogue_neighbour_marking]
    """
    ny, nx = image_data.shape
    half = stamp_size // 2
    col, row = _centre_pixels(x_pos, y_pos)
    vignets = np.full((len(col), stamp_size, stamp_size), BIG, np.float32)
    seg_stamp = np.zeros((stamp_size, stamp_size), np.int32)

    for i, (xi, yi) in enumerate(zip(col, row)):
        x0, y0 = xi - half, yi - half
        xc0, xc1 = max(0, x0), min(nx, x0 + stamp_size)
        yc0, yc1 = max(0, y0), min(ny, y0 + stamp_size)
        if xc0 >= xc1 or yc0 >= yc1:
            continue
        inner = np.s_[yc0 - y0:yc1 - y0, xc0 - x0:xc1 - x0]
        vignets[i][inner] = image_data[yc0:yc1, xc0:xc1]
        if seg is not None:
            seg_stamp[:] = 0
            seg_stamp[inner] = seg[yc0:yc1, xc0:xc1]
            vignets[i][(seg_stamp != 0) & (seg_stamp != number[i])] = BIG

    return vignets


def measure_windowed_positions(cat_data, image_data, img_header, seg=None):
    """``XWIN_*``, ``YWIN_*`` and ``FLAGS_WIN`` for a catalogue without them.

    The catalogue is the sample; its isophotal barycentre (``X_IMAGE``,
    ``Y_IMAGE``) and half-light radius (``FLUX_RADIUS``) start and size
    SExtractor's windowed centroid, measured on the tile image less a
    SExtractor-like background map (``default_tile.sex``: ``BACK_SIZE``
    512, ``BACK_FILTERSIZE`` 9). Pixels exactly 0 or not finite carry no
    data in a MegaPipe tile: they are left out of the background and
    replaced by their mirror image in the window. Given the relabelled
    segmentation map, neighbours are masked as SExtractor's ``MASK_TYPE
    CORRECT`` masks them. Where the window fails, the position stays the
    barycentre and ``FLAGS_WIN`` says why
    (:mod:`windowed_position`). World positions come from the tile WCS,
    1-based as SExtractor's are.

    Parameters
    ----------
    cat_data : astropy.table.Table
        Catalogue with ``NUMBER``, ``X_IMAGE``, ``Y_IMAGE``, ``FLUX_RADIUS``,
        and optionally ``A_WORLD``, ``B_WORLD``, ``THETA_J2000``
    image_data : numpy.ndarray
        Tile image
    img_header : astropy.io.fits.Header
        Tile header, for the WCS
    seg : numpy.ndarray, optional
        Segmentation map relabelled to ``NUMBER`` (:func:`relabel_seg`)

    Returns
    -------
    list of (str, numpy.ndarray, str)
        ``XWIN_IMAGE``, ``YWIN_IMAGE``, ``XWIN_WORLD``, ``YWIN_WORLD``,
        ``FLAGS_WIN`` columns with their FITS formats

    @sc [decision:preparation.object_position_columns]
    """
    missing = {"X_IMAGE", "Y_IMAGE", "FLUX_RADIUS"} - set(cat_data.colnames)
    if missing:
        raise ValueError(
            f"The catalogue lacks {sorted(missing)}, which the windowed"
            + " centroid needs."
        )
    x = np.asarray(cat_data["X_IMAGE"], dtype=np.float64)
    y = np.asarray(cat_data["Y_IMAGE"], dtype=np.float64)
    wcs = WCS(img_header)
    cxx = cyy = cxy = None
    if {"A_WORLD", "B_WORLD", "THETA_J2000"} <= set(cat_data.colnames):
        cxx, cyy, cxy = wp.isophotal_ellipse(
            cat_data["A_WORLD"], cat_data["B_WORLD"], cat_data["THETA_J2000"],
            wcs, x, y,
        )
    has_data = (image_data != 0) & np.isfinite(image_data)
    measured = image_data - wp.sextractor_background(
        image_data, good=has_data
    )
    xwin, ywin, flags = wp.windowed_positions(
        measured, x, y, cat_data["FLUX_RADIUS"], seg=seg,
        number=np.asarray(cat_data["NUMBER"]), cxx=cxx, cyy=cyy, cxy=cxy,
        valid=has_data,
    )
    ra, dec = wcs.all_pix2world(xwin, ywin, 1)
    return [
        ("XWIN_IMAGE", xwin, "D"),
        ("YWIN_IMAGE", ywin, "D"),
        ("XWIN_WORLD", ra, "D"),
        ("YWIN_WORLD", dec, "D"),
        ("FLAGS_WIN", flags, "I"),
    ]


def _build_ldac_objects(cat_data, vignets, windowed=()):
    """Build LDAC_OBJECTS extension from an astropy table.

    Parameters
    ----------
    cat_data : astropy.table.Table
        Catalogue data read from the ASCII SExtractor file
    vignets : numpy.ndarray
        Array of shape ``(n_obj, stamp_size, stamp_size)``
    windowed : sequence of (str, numpy.ndarray, str), optional
        Extra ``(name, array, FITS format)`` columns, the windowed positions

    Returns
    -------
    astropy.io.fits.BinTableHDU
        LDAC_OBJECTS extension

    """
    stamp_size = vignets.shape[1]
    n_pix = stamp_size * stamp_size

    fits_cols = []
    for colname in cat_data.colnames:
        arr = np.array(cat_data[colname])
        kind = arr.dtype.kind
        if kind in ("i", "u"):
            fmt = "J"
        elif kind == "f":
            fmt = "E" if arr.dtype.itemsize <= 4 else "D"
        else:
            fmt = f"{arr.dtype.itemsize}A"
        fits_cols.append(fits.Column(name=colname, format=fmt, array=arr))

    for name, arr, fmt in windowed:
        fits_cols.append(fits.Column(name=name, format=fmt, array=arr))

    vignet_col = fits.Column(
        name="VIGNET",
        format=f"{n_pix}E",
        array=vignets.reshape(len(vignets), n_pix),
        dim=f"({stamp_size},{stamp_size})",
    )
    fits_cols.append(vignet_col)

    hdu = fits.BinTableHDU.from_columns(fits_cols)
    hdu.name = "LDAC_OBJECTS"
    return hdu


def make_ldac_from_ascii(
    input_cat_path,
    image_path,
    output_cat_path,
    stamp_size=51,
    seg_path=None,
    w_log=None,
):
    """Convert an external ASCII catalogue to FITS-LDAC format.

    Reads an ASCII catalogue in SExtractor format (column definitions as
    ``#  N  NAME  ...`` comment lines followed by data rows), attaches the
    tile image header as an ``LDAC_IMHEAD`` extension so that
    ``make_post_process`` can recover the contributing single exposures, and
    writes a standard FITS-LDAC file compatible with all downstream ShapePipe
    modules.

    The input columns, ``NUMBER`` included, are copied unchanged; a
    ``VIGNET`` column (postage stamps extracted from the tile image) is
    added to ``LDAC_OBJECTS``, and, when the catalogue has no windowed
    positions, ``XWIN_IMAGE``, ``YWIN_IMAGE``, ``XWIN_WORLD``,
    ``YWIN_WORLD`` and ``FLAGS_WIN`` measured on the tile
    (:func:`measure_windowed_positions`). Given the catalogue's segmentation map, which
    shares the tile's pixel grid, the map is relabelled to the catalogue's
    ``NUMBER`` (:func:`relabel_seg`) and used to set neighbours' pixels in each ``VIGNET`` to ``BIG``, as
    SExtractor does.

    Parameters
    ----------
    input_cat_path : str
        Path to input ASCII catalogue in SExtractor format
    image_path : str
        Path to tile image FITS file
    output_cat_path : str
        Path to the output FITS-LDAC catalogue
    stamp_size : int, optional
        Side length of the square postage stamp in pixels, default 51
    seg_path : str, optional
        Path to the catalogue's segmentation map (FITS, compressed or not)
    w_log : logging.Logger, optional
        Pipeline logger

    """
    cat_data = asc.read(input_cat_path, format="sextractor")
    n_obj = len(cat_data)
    if w_log:
        w_log.info(f"Read {n_obj} objects from {input_cat_path}")

    with fits.open(image_path) as hdul:
        img_header = hdul[0].header
        image_data = hdul[0].data.astype(np.float32)

    seg = None
    if seg_path is not None:
        with fits.open(seg_path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            seg_raw = hdu.data
        if seg_raw.shape != image_data.shape:
            raise ValueError(
                f"Segmentation map {seg_path} has shape {seg_raw.shape}, the"
                + f" image {image_data.shape}; they must share one grid."
            )
        seg, counts = relabel_seg(
            seg_raw, cat_data["NUMBER"], cat_data["X_IMAGE"],
            cat_data["Y_IMAGE"],
        )
        del seg_raw
        if w_log:
            w_log.info(
                f"Relabelled {seg_path} to NUMBER: "
                + ", ".join(f"{k}={v}" for k, v in counts.items())
            )

    if w_log:
        w_log.info(
            f"Extracting {stamp_size}x{stamp_size} vignets from {image_path}"
        )
    vignets = _extract_vignets(
        image_data,
        cat_data["X_IMAGE"],
        cat_data["Y_IMAGE"],
        stamp_size,
        seg=seg,
        number=np.asarray(cat_data["NUMBER"]),
    )

    windowed = []
    if "XWIN_IMAGE" not in cat_data.colnames:
        windowed = measure_windowed_positions(
            cat_data, image_data, img_header, seg=seg
        )
        if w_log and windowed:
            flags = windowed[-1][1]
            w_log.info(
                "Measured windowed positions; FLAGS_WIN counts: "
                + ", ".join(
                    f"{v}={c}" for v, c in zip(*np.unique(flags,
                                                          return_counts=True))
                )
            )

    ldac_imhead = _build_ldac_imhead(img_header)
    ldac_objects = _build_ldac_objects(cat_data, vignets, windowed)

    hdul_out = fits.HDUList([fits.PrimaryHDU(), ldac_imhead, ldac_objects])
    hdul_out.writeto(output_cat_path, overwrite=True)

    if w_log:
        w_log.info(f"Written FITS-LDAC catalogue to {output_cat_path}")
