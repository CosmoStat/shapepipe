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
from shapepipe.modules.get_images_runner import in2out_pattern


def main(argv=None):
    """Main program.
    """

    # Parameters
    dirs_Git = glob.glob('output/run_sp_Git_*')
    pattern = 'CFIS.V0.skycell.'
    ext = 'fits'
    hdu_no = 1

    output_base = 'UNIONS_'
    pattern_map = {
        'unconv' : '_image-',
        'wt' : '_weight-',
        'mask' : '_flag-'
    }


    last_Git = False

    if last_Git:
        print('Converting PS image names in last Git run dir')
    else:
        print('Converting PS image names in all Git run dirs')

    for dir_Git in dirs_Git:

        if last_Git:
            # Only process last Git run
            if dir_Git != dirs_Git[-1]:
                continue

        dir_Git_out = '{}/get_images_runner/output'.format(dir_Git)
        dir_Git_out = os.path.abspath(dir_Git_out)
        print(dir_Git_out)

        ps_fnames = glob.glob('{}/{}*.{}'
                              ''.format(dir_Git_out, pattern, ext))
        for psfn in ps_fnames:

            # Get filter name
            header = fits.getheader(psfn, hdu_no) 
            filter_long = header['HIERARCH FPA.FILTERID']
            filter_letter = filter_long[0]

            # Get tile ID
            input_name = os.path.basename(psfn)
            m = re.match('.*(\d{3}\.\d{3}).*', input_name)
            number = m[1]
            number_final = in2out_pattern(number)

            # Get image type
            mm = re.match('.*\.(.*)\.fits', input_name)
            output_type = pattern_map[mm[1]]

            # Assemble output file name
            output_name = '{}{}{}{}.{}'.format(output_base, filter_letter,
                                               output_type, number_final, ext)

            output_path = '{}/{}'.format(dir_Git_out, output_name)
            #os.rename(psfn, output_path)
            try:
                print(' {} -> {}'.format(input_name, output_name))
                os.symlink(psfn, output_path)
            except FileExistsError:
                print(' {} already exists, skipping'.format(output_name))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
