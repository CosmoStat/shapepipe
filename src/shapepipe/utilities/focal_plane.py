"""FOCAL PLANE GEOMETRY.

The sky footprint of a MegaCam exposure, read from its image headers. Both
star-catalogue producers need exactly this and used to carry their own copy of
it -- ``workflow/scripts/star_cats.py`` (the HEALPix chunk store's ``cut``) and
``scripts/python/create_star_cat.py`` (the one-cone-per-exposure path the store
replaced). They must agree on which sky an exposure covers, so there is one
definition of it.

:Author: consolidated from the two star-catalogue scripts

"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS


def get_wcs(header):
    """Build the WCS by hand, from the linear terms only.

    Deliberately NOT ``WCS(header)``: it sidesteps distortion-convention
    incompatibilities between these headers and astropy, and a footprint needs
    nothing finer than the linear terms.

    Parameters
    ----------
    header : astropy.io.fits.Header
        Image header

    Returns
    -------
    astropy.wcs.WCS
        WCS object

    """
    final_wcs = WCS(naxis=2)
    final_wcs.wcs.ctype = [header["CTYPE1"], header["CTYPE2"]]
    try:
        final_wcs.wcs.cunit = [header["CUNIT1"], header["CUNIT2"]]
    except KeyError:
        final_wcs.wcs.cunit = ["deg", "deg"]
    final_wcs.wcs.crpix = [header["CRPIX1"], header["CRPIX2"]]
    final_wcs.wcs.crval = [header["CRVAL1"], header["CRVAL2"]]
    final_wcs.wcs.cd = [
        [header["CD1_1"], header["CD1_2"]],
        [header["CD2_1"], header["CD2_2"]],
    ]

    return final_wcs


def ccd_center_and_radius(header):
    """Return ``(ra_deg, dec_deg, radius_deg)`` for a single CCD.

    The radius is the half-diagonal: centre to the ``(0, 0)`` corner.

    """
    w = get_wcs(header)
    (ra_c, dec_c), (ra_0, dec_0) = w.all_pix2world(
        [[header["NAXIS1"] / 2.0, header["NAXIS2"] / 2.0], [0, 0]], 1)
    center = SkyCoord(ra_c * u.deg, dec_c * u.deg)
    radius = center.separation(SkyCoord(ra_0 * u.deg, dec_0 * u.deg)).deg
    return float(ra_c), float(dec_c), float(radius)


def focal_plane_disc(image, n_ccd=40):
    """Return ``(ra_deg, dec_deg, radius_deg)`` covering all CCDs of one exposure.

    The centre is the mean of the CCD centres; the radius is the largest
    centre-to-CCD-centre distance plus that CCD's own half-diagonal.

    ONE ``fits.open`` for all the extensions: ``fits.getheader(image, ext)``
    opens the file, walks the HDU list to ``ext`` and closes again, so a loop
    over it costs ``n_ccd`` opens and O(n^2) header seeks.

    """
    centers, radii = [], []
    with fits.open(image) as hdul:
        for ext in range(1, n_ccd + 1):
            ra_c, dec_c, radius = ccd_center_and_radius(hdul[ext].header)
            centers.append((ra_c, dec_c))
            radii.append(radius)

    ras = np.array([c[0] for c in centers])
    decs = np.array([c[1] for c in centers])
    center = SkyCoord(ras.mean() * u.deg, decs.mean() * u.deg)
    seps = center.separation(SkyCoord(ras * u.deg, decs * u.deg)).deg
    return (float(center.ra.deg), float(center.dec.deg),
            float(np.max(seps + np.array(radii))))
