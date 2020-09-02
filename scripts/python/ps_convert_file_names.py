#!/usr/bin/env python


"""Script ps_convert_file_names.py

Convert file names of (downloaded from canfar) Pan-STARRS
file names, according to filter read from FITS header.

:Authors: Martin Kilbinger

:Date: 02/09/2020
"""


# Compability with python2.x for x>6
from __future__ import print_function

import sys
import os
import re
import copy
import glob

import numpy as np
import pylab as plt

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup

from astropy.io import fits
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle, SkyCoord

from shapepipe.pipeline import file_io as io


def main(argv=None):
    """Main program.
    """

    # Parameters
    dirs_Git = glob.glob('output/run_sp_Git_*')
    pattern = 'CFIS.V0.skycell.'
    ext = 'fits'
    hdu_no = 1

    # Get input file names
    last_Git = '{}/get_images_runner/output'.format(dirs_Git[-1])

    ps_files = glob.glob('{}/{}*.{}'.format(last_Git, pattern, ext))

    # Read headers
    for psf in ps_files:
        header = fits.getheader(psf, hdu_no) 
        filter_long = header['HIERARCH FPA.FILTERID']
        filter_letter = filter_long[0]
        print(filter_letter)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
