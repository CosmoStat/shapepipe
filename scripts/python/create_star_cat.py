#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script create_star_cat.py

:Description: Create reference star catalogue for masking of
bright star halos and diffraction spikes

:Authors: Axel Guinot, Martin Kilbinger

"""


import os
import random
import re
import sys
import time

from cs_util import args as cs_args
from cs_util import logging as cs_logging

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from astroquery.vizier import Vizier


# GSC 2.3 catalog ID
CDS_CAT_ID = "I/305/out"

VIZIER_SERVERS = [
    "vizier.cds.unistra.fr",
    "vizier.cfa.harvard.edu",
    "vizier.iucaa.in",
]

VIZIER_TIMEOUTS = [10, 20, 40]


def _get_wcs(header):
    """Get WCS.

    Compute the astropy WCS from header manually.
    (The purpose of this is to avoid possible incompatibility on distortion
    convention)

    Parameters
    ----------
    header : astropy.header
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
    except:
        final_wcs.wcs.cunit = ["deg", "deg"]
    final_wcs.wcs.crpix = [header["CRPIX1"], header["CRPIX2"]]
    final_wcs.wcs.crval = [header["CRVAL1"], header["CRVAL2"]]
    final_wcs.wcs.cd = [
        [header["CD1_1"], header["CD1_2"]],
        [header["CD2_1"], header["CD2_2"]],
    ]

    return final_wcs


def _sphere_dist_arcmin(ra1, dec1, ra2, dec2):
    """Compute angular distance between two sky positions in arcmin."""
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)
    return c1.separation(c2).arcmin


def _ccd_center_and_radius(header):
    """Return (ra, dec, radius_arcmin) for a single CCD."""
    w = _get_wcs(header)
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    cx, cy = nx / 2.0, ny / 2.0
    ra_c, dec_c = w.all_pix2world([[cx, cy]], 1)[0]
    # radius = half-diagonal to CCD corner
    ra_corner, dec_corner = w.all_pix2world([[0, 0]], 1)[0]
    radius = _sphere_dist_arcmin(ra_c, dec_c, ra_corner, dec_corner)
    return ra_c, dec_c, radius


def _focal_plane_center_and_radius(f, n_ccd=40):
    """Return (ra, dec, radius_arcmin) covering all CCDs of an exposure."""
    ras, decs, radii = [], [], []
    for ind in range(1, n_ccd + 1):
        h = fits.getheader(f, ind)
        ra, dec, r = _ccd_center_and_radius(h)
        ras.append(ra)
        decs.append(dec)
        radii.append(r)

    ras = np.array(ras)
    decs = np.array(decs)
    radii = np.array(radii)

    ra_center = np.mean(ras)
    dec_center = np.mean(decs)

    # Radius = max distance from focal plane center to any CCD center + that CCD's half-diagonal
    dists = np.array([
        _sphere_dist_arcmin(ra_center, dec_center, ras[i], decs[i])
        for i in range(len(ras))
    ])
    radius = np.max(dists + radii)

    return ra_center, dec_center, radius


def query_vizier(ra, dec, radius_arcmin):
    """Query Vizier (GSC 2.3) with retries over timeouts and servers.

    Parameters
    ----------
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    radius_arcmin : float
        Search radius in arcminutes

    Returns
    -------
    astropy.table.Table
        Star catalog table

    """
    p = np.array([ra, dec], dtype="single")
    coord = SkyCoord(ra=p[0] * u.deg, dec=p[1] * u.deg, frame="icrs")

    time.sleep(random.uniform(0, 5))

    result = []
    for attempt, timeout in enumerate(VIZIER_TIMEOUTS):
        for server in VIZIER_SERVERS:
            v = Vizier(row_limit=-1, timeout=timeout, vizier_server=server)
            result = v.query_region(
                coord, radius=radius_arcmin * u.arcmin, catalog=CDS_CAT_ID
            )
            if len(result) > 0:
                break
            print(
                f"Vizier returned empty list at {coord}, {radius_arcmin:.2f} arcmin, "
                f"server={server}, timeout={timeout}s"
            )
        if len(result) > 0:
            print("Vizier query successful")
            break
        if attempt < len(VIZIER_TIMEOUTS) - 1:
            wait = 10 * 2**attempt
            print(
                f"All servers failed, retrying in {wait}s "
                f"(attempt {attempt+1}/{len(VIZIER_TIMEOUTS)}, "
                f"next timeout={VIZIER_TIMEOUTS[attempt+1]}s)"
            )
            time.sleep(wait)

    if len(result) == 0:
        raise RuntimeError(
            f"Vizier astroquery returned empty list at {coord}, {radius_arcmin:.2f} arcmin"
        )

    return result[0]


def main(input_dir, output_dir, kind):

    file_list = os.listdir(input_dir)

    for f in file_list:
        if "image" not in f:
            continue

        img_number = re.split("image", os.path.splitext(f)[0])[1]
        fpath = os.path.join(input_dir, f)

        if kind == "exp":
            # One query covering the full MegaCam focal plane
            output_name = f"{output_dir}/star_cat{img_number}.fits"
            if os.path.isfile(output_name):
                continue

            ra, dec, radius = _focal_plane_center_and_radius(fpath)
            print(
                f"Focal plane center: ra={ra:.4f}, dec={dec:.4f}, radius={radius:.2f} arcmin"
            )
            table = query_vizier(ra, dec, radius)
            table.write(output_name, overwrite=True)

        else:
            h = fits.getheader(fpath, 0)
            w = _get_wcs(h)
            nx, ny = h["NAXIS1"], h["NAXIS2"]
            cx, cy = nx / 2.0, ny / 2.0
            ra, dec = w.all_pix2world([[cx, cy]], 1)[0]
            ra_corner, dec_corner = w.all_pix2world([[0, 0]], 1)[0]
            radius = _sphere_dist_arcmin(ra, dec, ra_corner, dec_corner)

            output_name = f"{output_dir}/star_cat{img_number}.fits"
            if os.path.isfile(output_name):
                continue

            table = query_vizier(ra, dec, radius)
            table.write(output_name, overwrite=True)

    return 0


def params_default():
    """Return default parameters, short options, types, and help strings."""
    _params = {
        "input_dir": ".",
        "output_dir": ".",
        "kind": "exp",
    }
    _short_options = {
        "input_dir": "-i",
        "output_dir": "-o",
        "kind": "-k",
    }
    _types = {}
    _help_strings = {
        "input_dir": "input directory containing image files; default is {}",
        "output_dir": "output directory for star catalogues; default is {}",
        "kind": "processing kind, 'exp' for full MegaCam focal plane, 'tile' for single image; default is {}",
    }
    return _params, _short_options, _types, _help_strings


if __name__ == "__main__":

    _params, _short_options, _types, _help_strings = params_default()

    options = cs_args.parse_options(
        _params,
        _short_options,
        _types,
        _help_strings,
    )

    cs_logging.log_command(sys.argv)

    main(options["input_dir"], options["output_dir"], options["kind"])
