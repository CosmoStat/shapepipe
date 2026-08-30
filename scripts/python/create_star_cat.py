#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script create_star_cat.py

:Description: Create reference star catalogue for masking of
bright star halos and diffraction spikes

:Authors: Axel Guinot, Martin Kilbinger

"""


import os
import re
import sys

from cs_util import args as cs_args
from cs_util import logging as cs_logging

from astropy.io import fits

from shapepipe.utilities.file_io import write_atomic
from shapepipe.utilities.focal_plane import ccd_center_and_radius, focal_plane_disc
from shapepipe.utilities.vizier import query_vizier as _query_vizier


# GSC 2.3 catalog ID
CDS_CAT_ID = "I/305/out"


def query_vizier(ra, dec, radius_arcmin):
    return _query_vizier(ra, dec, radius_arcmin, CDS_CAT_ID)


def main(input_dir, output_dir, kind):

    file_list = os.listdir(input_dir)

    for f in file_list:
        if "image" not in f:
            continue

        img_number = re.split("image", os.path.splitext(f)[0])[1]
        fpath = os.path.join(input_dir, f)

        output_name = f"{output_dir}/star_cat{img_number}.fits"
        if os.path.isfile(output_name):
            continue

        if kind == "exp":
            # One query covering the full MegaCam focal plane
            ra, dec, radius_deg = focal_plane_disc(fpath)
            radius = radius_deg * 60.0
            print(
                f"Focal plane center: ra={ra:.4f}, dec={dec:.4f}, radius={radius:.2f} arcmin"
            )
        else:
            # A single image: its own centre and half-diagonal.
            ra, dec, radius_deg = ccd_center_and_radius(fits.getheader(fpath, 0))
            radius = radius_deg * 60.0

        table = query_vizier(ra, dec, radius)
        write_atomic(table, output_name)

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
