"""READ EXTERNAL SEXCAT.

Convert an ASCII SExtractor-format catalogue to FITS-LDAC for use in
ShapePipe tile processing.

:Author: Martin Kilbinger

"""

import numpy as np
from astropy.io import ascii as asc
from astropy.io import fits


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


def _extract_vignets(image_data, x_pos, y_pos, stamp_size):
    """Extract postage stamps from a tile image array.

    For each object position, a ``stamp_size x stamp_size`` cutout is
    extracted.  Objects whose stamp falls partially outside the image are
    padded with zeros, matching SExtractor behaviour.

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

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n_obj, stamp_size, stamp_size)``, dtype float32

    """
    ny, nx = image_data.shape
    half = stamp_size // 2
    n_obj = len(x_pos)
    vignets = np.zeros((n_obj, stamp_size, stamp_size), dtype=np.float32)

    for i, (x, y) in enumerate(zip(x_pos, y_pos)):
        xi = int(round(float(x))) - 1
        yi = int(round(float(y))) - 1

        x0, x1 = xi - half, xi + half + 1
        y0, y1 = yi - half, yi + half + 1

        xc0, xc1 = max(0, x0), min(nx, x1)
        yc0, yc1 = max(0, y0), min(ny, y1)

        dx0, dx1 = xc0 - x0, xc0 - x0 + (xc1 - xc0)
        dy0, dy1 = yc0 - y0, yc0 - y0 + (yc1 - yc0)

        vignets[i, dy0:dy1, dx0:dx1] = image_data[yc0:yc1, xc0:xc1]

    return vignets


def _build_ldac_objects(cat_data, vignets):
    """Build LDAC_OBJECTS extension from an astropy table.

    Parameters
    ----------
    cat_data : astropy.table.Table
        Catalogue data read from the ASCII SExtractor file
    vignets : numpy.ndarray
        Array of shape ``(n_obj, stamp_size, stamp_size)``

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

    _aliases = {
        "XWIN_IMAGE": "X_IMAGE",
        "YWIN_IMAGE": "Y_IMAGE",
        "XWIN_WORLD": "ALPHA_J2000",
        "YWIN_WORLD": "DELTA_J2000",
    }
    col_map = {c.name: c for c in fits_cols}
    for alias, source in _aliases.items():
        if source in col_map and alias not in col_map:
            src_col = col_map[source]
            fits_cols.append(
                fits.Column(
                    name=alias,
                    format=src_col.format,
                    array=np.array(cat_data[source]),
                )
            )

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
    added to ``LDAC_OBJECTS``.

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

    if w_log:
        w_log.info(
            f"Extracting {stamp_size}x{stamp_size} vignets from {image_path}"
        )
    vignets = _extract_vignets(
        image_data,
        cat_data["X_IMAGE"],
        cat_data["Y_IMAGE"],
        stamp_size,
    )

    ldac_imhead = _build_ldac_imhead(img_header)
    ldac_objects = _build_ldac_objects(cat_data, vignets)

    hdul_out = fits.HDUList([fits.PrimaryHDU(), ldac_imhead, ldac_objects])
    hdul_out.writeto(output_cat_path, overwrite=True)

    if w_log:
        w_log.info(f"Written FITS-LDAC catalogue to {output_cat_path}")
