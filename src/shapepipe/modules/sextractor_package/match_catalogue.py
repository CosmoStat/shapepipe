"""MATCH CATALOGUE.

Join a SExtractor tile catalogue to an external catalogue of the same image,
so the tile's rows carry the external catalogue's object numbers.

The UNIONS tile catalogue (Stephen Gwyn's MegaPipe SExtractor run on the DR6
tile) defines the object list and the ``NUMBER`` that the photometry and
photo-z catalogues share. ShapePipe runs its own SExtractor on the same image
with the same detection configuration, which reproduces that catalogue for
all but a few objects in a thousand, and takes membership and ``NUMBER`` from
it here: every
measurement (windowed positions, ``VIGNET`` and its neighbour marking, the
SExtractor columns) stays SExtractor's own.

:Author: Cail Daley

"""

import numpy as np
from astropy.io import ascii as asc
from astropy.io import fits
from scipy.spatial import cKDTree

# The SEG_VIGNET label of a footprint whose SExtractor row has no partner in
# the external catalogue, and so leaves the catalogue. Negative, so it never
# collides with a NUMBER (match_catalogue requires positive ones); UberSeg only
# asks "self or not self".
UNMATCHED_LABEL = -1

# The largest NUMBER the int32 NUMBER and SEG_VIGNET columns hold.
MAX_NUMBER = np.iinfo(np.int32).max


def mutual_nearest(x_a, y_a, x_b, y_b, radius):
    """One-to-one pairs of mutual nearest neighbours closer than ``radius``.

    Parameters
    ----------
    x_a, y_a : array_like
        Positions of the first set
    x_b, y_b : array_like
        Positions of the second set, in the same units
    radius : float
        Largest separation of a pair

    Returns
    -------
    numpy.ndarray
        Indices into the first set of the paired rows, ascending
    numpy.ndarray
        Index into the second set of each one's partner
    numpy.ndarray
        Separation of each pair

    """
    a = np.column_stack([x_a, y_a]).astype(np.float64)
    b = np.column_stack([x_b, y_b]).astype(np.float64)
    if len(a) == 0 or len(b) == 0:
        empty = np.zeros(0, np.int64)
        return empty, empty, np.zeros(0)
    dist, nearest_b = cKDTree(b).query(a)
    _, nearest_a = cKDTree(a).query(b)
    idx_a = np.arange(len(a))
    paired = (dist < radius) & (nearest_a[nearest_b] == idx_a)
    return idx_a[paired], nearest_b[paired], dist[paired]


def relabel(seg_vignets, old_number, new_number):
    """Map segmentation labels from one numbering to another.

    Parameters
    ----------
    seg_vignets : numpy.ndarray
        Segmentation stamps labelled with the old numbering, 0 on sky
    old_number : array_like of int
        Every old label in use (the catalogue's ``NUMBER`` before the join)
    new_number : array_like of int
        New label of each ``old_number``; ``UNMATCHED_LABEL`` for one that
        leaves the catalogue

    Returns
    -------
    numpy.ndarray
        The stamps in the new numbering, int32; labels outside
        ``old_number`` become ``UNMATCHED_LABEL``

    """
    old_number = np.asarray(old_number, np.int64)
    table = np.full(int(max(old_number.max(initial=0),
                            seg_vignets.max(initial=0))) + 1,
                    UNMATCHED_LABEL, np.int32)
    table[0] = 0
    table[old_number] = new_number
    if seg_vignets.min(initial=0) < 0:
        raise ValueError("Segmentation stamps hold negative labels.")
    return table[seg_vignets]


def match_catalogue(cat_path, ext_cat_path, radius=1.0, min_fraction=0.98,
                    tolerated_unpaired=20, w_log=None):
    """Take membership and ``NUMBER`` from an external catalogue.

    Each SExtractor row is paired with its mutual nearest neighbour in the
    external catalogue (``X_IMAGE``, ``Y_IMAGE``, same image grid) within
    ``radius`` pixels. Paired rows take the external ``NUMBER``; unpaired
    rows leave the catalogue. A ``SEG_VIGNET`` column, when present, is
    relabelled through the whole map from old to new numbers, with the
    footprints of rows that left marked ``UNMATCHED_LABEL``: as the external
    numbers are unique integers in ``[1, MAX_NUMBER]``, each row's own
    footprint carries its new ``NUMBER`` and no other footprint can,
    whatever the two numberings share.
    The catalogue is rewritten in place; every other HDU and column is kept.

    Both catalogues come from the same pixels, so pairs agree to ~1e-4
    pixel, and on eight DR6 tiles at least 99.2% of each side pairs; the
    shortfall is deblending that differs near very large objects (a cD
    galaxy, a bright star's halo). On different pixels (a DR5 image against
    the DR6 catalogue) only ~95% pair. The guard is there to catch that: the
    run stops when the unpaired rows of either side exceed both
    ``(1 - min_fraction)`` of that side and ``tolerated_unpaired``, so a
    near-empty edge tile with a couple of unpaired children passes. On the
    SExtractor side it catches a catalogue from other pixels; on the
    external side, a detection configuration that finds fewer objects than
    the catalogue's, which would leave its objects without shapes.

    Parameters
    ----------
    cat_path : str
        Path to the SExtractor FITS-LDAC catalogue, rewritten in place
    ext_cat_path : str
        Path to the external ASCII catalogue in SExtractor format, with
        ``NUMBER``, ``X_IMAGE`` and ``Y_IMAGE``
    radius : float, optional
        Largest separation of a pair, in pixels; default 1
    min_fraction : float, optional
        Smallest acceptable fraction of SExtractor rows, and of external
        objects, with a partner; default 0.98
    tolerated_unpaired : int, optional
        Number of unpaired rows on either side accepted whatever the
        fraction; default 20
    w_log : logging.Logger, optional
        Pipeline logger

    Returns
    -------
    dict
        ``n_sextractor``, ``n_external``, ``n_paired``, ``n_dropped``
        (SExtractor rows without a partner) and ``n_external_only``
        (external objects without a SExtractor row)

    Raises
    ------
    ValueError
        If the external ``NUMBER`` repeats or is not an integer in
        ``[1, MAX_NUMBER]``, or if more
        than ``tolerated_unpaired`` rows of either side, and more than
        ``1 - min_fraction`` of it, have no partner

    @sc [decision:detection.tile_detection]
    """
    ext = asc.read(ext_cat_path, format="sextractor",
                   include_names=["NUMBER", "X_IMAGE", "Y_IMAGE"])
    ext_number = np.asarray(ext["NUMBER"])
    if (not np.issubdtype(ext_number.dtype, np.integer)
            or (ext_number <= 0).any() or (ext_number > MAX_NUMBER).any()
            or len(np.unique(ext_number)) < len(ext)):
        raise ValueError(
            f"{ext_cat_path} has a NUMBER that is not an integer in"
            + f" [1, {MAX_NUMBER}] or that repeats; the join needs unique"
            + " numbers that fit the int32 NUMBER and SEG_VIGNET columns, so"
            + " that no relabelled footprint takes another object's number."
        )
    with fits.open(cat_path) as hdul:
        hdus = [hdu.copy() for hdu in hdul]
    objects = next(h for h in hdus if h.name == "LDAC_OBJECTS")
    data = objects.data

    i_sex, i_ext, dist = mutual_nearest(
        data["X_IMAGE"], data["Y_IMAGE"], ext["X_IMAGE"], ext["Y_IMAGE"],
        radius,
    )
    counts = dict(
        n_sextractor=len(data), n_external=len(ext), n_paired=len(i_sex),
        n_dropped=len(data) - len(i_sex),
        n_external_only=len(ext) - len(i_sex),
    )
    summary = ", ".join(f"{k}={v}" for k, v in counts.items())
    def too_many_unpaired(n_total):
        unpaired = n_total - len(i_sex)
        return (unpaired > tolerated_unpaired
                and unpaired > (1 - min_fraction) * n_total)

    fraction = len(i_sex) / max(len(data), 1)
    if too_many_unpaired(len(data)):
        raise ValueError(
            f"Only {fraction:.4f} of the SExtractor rows of {cat_path} pair"
            + f" with {ext_cat_path} within {radius} px (minimum"
            + f" {min_fraction} beyond {tolerated_unpaired} unpaired;"
            + f" {summary}). The catalogue must come from the"
            + " same pixels as the image: check that the tile image is the"
            + " release of the catalogue (DR6)."
        )
    fraction_ext = len(i_sex) / max(len(ext), 1)
    if too_many_unpaired(len(ext)):
        raise ValueError(
            f"Only {fraction_ext:.4f} of the objects of {ext_cat_path} pair"
            + f" with a row of {cat_path} within {radius} px (minimum"
            + f" {min_fraction} beyond {tolerated_unpaired} unpaired;"
            + f" {summary}). SExtractor finds fewer objects"
            + " than the catalogue: check that the detection configuration"
            + " is the catalogue's."
        )

    old_number = np.asarray(data["NUMBER"])
    new_number = np.full(len(data), UNMATCHED_LABEL, np.int64)
    new_number[i_sex] = ext_number[i_ext]

    columns = []
    for col in objects.columns:
        array = data[col.name][i_sex]
        if col.name == "NUMBER":
            array = new_number[i_sex].astype(array.dtype)
        elif col.name == "SEG_VIGNET":
            array = relabel(np.asarray(data[col.name]), old_number,
                            new_number)[i_sex]
        columns.append(fits.Column(
            name=col.name, format=col.format, unit=col.unit, dim=col.dim,
            array=array,
        ))
    new = fits.BinTableHDU.from_columns(
        columns, header=objects.header, name="LDAC_OBJECTS",
    )
    fits.HDUList(
        [new if h is objects else h for h in hdus]
    ).writeto(cat_path, overwrite=True)

    if w_log:
        sep = np.quantile(dist, [0.5, 0.99]) if len(dist) else [np.nan] * 2
        w_log.info(
            f"Matched {cat_path} to {ext_cat_path} within {radius} px:"
            + f" {summary}; pair separation p50 {sep[0]:.2e} px,"
            + f" p99 {sep[1]:.2e} px"
        )
    return counts
