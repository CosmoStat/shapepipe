"""CCD_FOOTPRINT

Project a single CCD onto the sky: image shape from a header, sky corners
from a WCS.

The two functions here are the geometry the coverage mask rests on, and they
are all that survives of the pre-Snakemake header scrape. The ``exp_footprint``
rule (``workflow/scripts/exp_footprint.py``) calls them once per CCD of an
exposure, reading the headers out of ``headers-<exp>.npy`` rather than out of
text files downloaded from VOSpace; the ``coverage_map`` rule then stamps the
corners into the campaign's HealSparse map.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""


def _image_shape(header):
    """Return the CCD image ``(nx, ny)`` pixel shape from a header.

    For fpack tile-compressed HDUs (``ZIMAGE = T``), ``NAXIS1``/``NAXIS2``
    describe the compressed *binary table* (row width in bytes, row count),
    not the image, and astropy's WCS never maps the true dimensions into
    ``pixel_shape``. The true image size is carried by ``ZNAXIS1``/``ZNAXIS2``,
    which are preferred here; a plain image HDU falls back to
    ``NAXIS1``/``NAXIS2``.

    Both branches are live. ``headers-<exp>.npy`` carries the *decompressed*
    header astropy hands back for a tile-compressed HDU, so it takes the plain
    ``NAXIS`` path; a header read straight off a fpacked file takes the ``Z*``
    one.

    Parameters
    ----------
    header : astropy.io.fits.Header
        single-HDU header

    Returns
    -------
    tuple
        ``(nx, ny)`` image pixel dimensions

    Raises
    ------
    ValueError
        if neither the ``Z*`` nor the plain ``NAXIS`` dimensions are present
    """
    if header.get("ZIMAGE", False):
        keys = ("ZNAXIS1", "ZNAXIS2")
    else:
        keys = ("NAXIS1", "NAXIS2")

    if keys[0] not in header or keys[1] not in header:
        raise ValueError(
            f"Header is missing image dimensions ({keys[0]}/{keys[1]})"
        )

    return int(header[keys[0]]), int(header[keys[1]])


def _ccd_corners(w, shape):
    """Return (ra, dec) of the 4 corners of a single CCD.

    The polygon is built from the CCD's pixel bounds ``shape = (nx, ny)`` using
    pixel edges ``(-0.5, n - 0.5)`` so the quadrilateral covers the full CCD
    area rather than pixel centres. Corners are returned in a consistent
    counterclockwise pixel order (bottom-left, bottom-right, top-right,
    top-left); since ``pixel_to_world`` maps this order to a simple,
    non-self-intersecting sky quadrilateral for a CCD-sized field, and
    HealSparse ``Polygon`` is orientation-agnostic, the resulting polygon is a
    valid convex footprint.

    Parameters
    ----------
    w : astropy.wcs.WCS
        WCS of a single CCD
    shape : tuple
        ``(nx, ny)`` CCD image pixel dimensions

    Returns
    -------
    tuple
        ``(ra_list, dec_list)`` each of length 4, in degrees
    """
    nx, ny = shape

    # Pixel-edge corners, counterclockwise: BL, BR, TR, TL
    px = [-0.5, nx - 0.5, nx - 0.5, -0.5]
    py = [-0.5, -0.5, ny - 0.5, ny - 0.5]

    sky = w.pixel_to_world(px, py)
    return list(sky.ra.deg), list(sky.dec.deg)
