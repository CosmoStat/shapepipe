#!/usr/bin/env python
  

"""Script cfis_check_weights.py

Checking weight files.

:Authors: Martin Kilbinger

:Date: 13/02/2018
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

from astropy.io import ascii
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup
import datetime

import cfis
import stuff



def get_file_list(pattern, verbose=False):
    """Read present working directory and return list of files accounting for input pattern.
    """

    files = glob.glob('*')

    pattern = cfis.get_file_pattern(pattern, 'r', 'weight')

    dst_list = []
    for f in files:
        m = re.findall(pattern, f)
        if len(m) != 0:
            dst_list.append(f)

    return dst_list


def diagnostics(files, verbose=False):

    limits = (0.0, 0.5)
    nbin   = 10

    diagn = []

    for f in files:

        if verbose:
            print('Checking file {}'.format(f))

        hdu  = fits.open(f)

        dat  = hdu[0].data
        y, x = np.histogram(dat, bins=nbin, range=limits)  

        # Ratio of numer of close-to-zero pixels to total
        ratio = y[0] / sum(y)

        # Diagnostic
        ok    = ratio > 0.01

        date  = hdu[0].header['DATE']

        hdu.close()

        diagn.append([f, ratio, ok, date])

    # Sort according to date in header
    diagn_s = sorted(diagn, key=lambda x: datetime.datetime.strptime(x[3], '%Y-%m-%dT%H:%M:%S'))

    return diagn_s



def output(diagn):

    print('{:30s} {:7s} {:15s} {:20s}'.format('# name', 'ratio[#0/all]', 'pass/fail', 'date'))
    for d in diagn:
        print('{:30s} {:13.5f} {:5s} {:20s}'.format(d[0], d[1], str(d[2]), d[3]))



def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class tuff.param
        parameter values
    """

    p_def = stuff.param(
        pattern  = '.',
    )

    return p_def



def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class tuff.param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage, formatter=stuff.IndentedHelpFormatterWithNL())

    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern,
         help='file pattern to match, e.g.~\'^21\d{5}p\', default=none (=all match)')

    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

    options, args = parser.parse_args()

    return options, args



def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    erg: bool
        Result of option check. False if invalid option value.
    """

    see_help = 'See option \'-h\' for help.'

    return True



def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class mkstuff.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class mkstuff.param
        updated paramter values
    """

    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    return param



def main(argv=None):
    """Main program.
    """


    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)


    # Save calling command
    stuff.log_command(argv)
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')


    if param.verbose is True:
        print('Start of program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    files = get_file_list(param.pattern, verbose=param.verbose)

    diagn = diagnostics(files, verbose=param.verbose)

    output(diagn)

    ### End main program

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

