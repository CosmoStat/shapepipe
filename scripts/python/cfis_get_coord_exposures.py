#!/usr/bin/env python

"""Script cfis_get_coord_exposures.py

Extracts coordinates from CFIS individual exposures (FITS image)

:Author: Martin Kilbinger

*Date*: 2018

"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import copy
import glob

import numpy as np
import pylab as plt

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column

from optparse import OptionParser
from optparse import OptionGroup

from astropy import units
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord

import stuff 



def get_pointings_list(input_names, verbose=False):
    """Return list of images.

    Parameters
    ----------
    input_names: string
        file names, wild cards possible
    verbose: bool, optional
        verbose output if True, default=False

    Return
    ------
    img_list: list of strings
        image list
    """

    file_list = glob.glob(input_names)

    if verbose == True:
        print('{} files found'.format(len(file_list)))

    return file_list



def get_coordinatess(pointings, verbose=False):
    """Return coordinates of pointings.

    Parameters
    ----------
    pointings: list of strings
        images names of pointings
    verbose: bool, optional, default=False
        verbose output

    Returns
    -------
    sc: list of SkyCoord
        list of coordinates
    """

    ihdu = 1

    keys  = ['RA_DEG', 'DEC_DEG']
    sc    = []

    for p in pointings:
        hdu = fits.open(p)
        header = hdu[ihdu].header
        coord = []
        for key in keys:
            if key in header:
                coord.append('{}deg'.format(header[key]))
            else:
                stuff.error('Key \'{}\' not found in header of file {}'.format(key, p))

        sc.append(SkyCoord(coord[0], coord[1]))

    return sc



def output(pointings, coords, output):
    """Print list of pointing file names and associated coordinates.

    Parameters
    ----------
    pointings: list of strings
        images names of pointings
    coords: list of SkyCoord
        list of coordinates

    Returns
    -------
    None
    """

    if output is None:
        fout = sys.stdout
    else:
        fout = open(output, 'w')

    print('# Pointing R.A.[degree] Declination[degree]', file=fout)
    for i, p in enumerate(pointings):
        pb = os.path.basename(p)
        print(pb, coords[i].ra.degree, coords[i].dec.degree, file=fout)



def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class stuff.param
        parameter values
    """

    p_def = stuff.param(
    )

    return p_def



def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class stuff.param
        parameter values

    Returns
    -------
    options: tuple
        Command line options
    args: string
        Command line string
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # I/O
    parser.add_option('-i', '--input_names', dest='input_names', type='string', default=None,
         help='input file name(s), wild cards possible')
    parser.add_option('-o', '--output', dest='output', type='string', default=None,
         help='output file name, default is stdout')

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

    if options.input_names is None:
        stuff.error("No input file names given")

    see_help = 'See option \'-h\' for help.'

    return True



def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class stuff.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class stuff.param
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

    # Do extra stuff if necessary
    # ...

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
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    pointings = get_pointings_list(param.input_names, verbose=param.verbose)

    coords    = get_coordinatess(pointings, verbose=param.verbose)

    output(pointings, coords, param.output)



    ### End main program

    if param.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

