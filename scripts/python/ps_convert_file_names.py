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

from optparse import OptionParser

from astropy.io import fits
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle, SkyCoord

from shapepipe.pipeline import file_io as io
from shapepipe.modules.get_images_runner import in2out_pattern

import cfis


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class cfis.param
        parameter values
    """

    p_def = cfis.param(
        indir = 'output/run_sp_Git_*', 
        last_Git = False,
    )

    return p_def


def parse_options(p_def):
    """Parse command line options.

    Parameters
    ----------
    p_def: class cfis.param
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
    parser.add_option('-i', '--indir', dest='indir', type='string', default=p_def.indir,
         help='input directory (pattern), default=\'{}\''.format(p_def.indir))
    parser.add_option('-o', '--outdir', dest='outdir', type='string', default=None,
         help='output directory, if not given: create links in input dir(s)')
    parser.add_option('-O', '--original', dest='outdir_orig', type='string', default=None,
         help='output dir for original file names; no output if not given')
    parser.add_option('-s', '--sp_format', dest='sp_format', action='store_true',
         help='output number in ShapePipe format (000-000)')
    parser.add_option('-l', '--last_Git', dest='last_Git', action='store_true',
        help='use only last run of \'get_images_runner\'')

    # Monitoring
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', help='verbose output')

    options, args = parser.parse_args()

    return options, args


def check_options(options):
    """Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Raises
    ------
    ValueError
    """

    if options.outdir_orig and not options.outdir:
        raise ValueError('Option \'-O\' only valid with \'-o\'')


def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class param
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


def ln_s(src, dst, src_base=None, verbose=False):
    """Create symbolic mimicking to shell cmd `ln -s`.

    Parameters
    ----------
    src : string
        source file path
    dst : string
        destination link name
    src_base : string, optional, default=None
        source base name, for verbose output only
    verbose : bool, optional, default=False
        verbose output if True
    """

    try:
       if verbose:
            if src_base:
                src_verbose = src_base
            else:
                src_verbose = src
            print(' {} -> {}'.format(src_verbose, dst))
       os.symlink(src, dst)
    except FileExistsError:
        if verbose:
            print(' {} already exists, skipping'.format(dst))


def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    check_options(options)

    param = update_param(p_def, options)

    # Save calling command
    cfis.log_command(argv)
    if param.verbose:
        cfis.log_command(argv, name='sys.stdout')

    if param.verbose:
        print('Start of program {}'.format(os.path.basename(argv[0])))


    # Parameters
    dirs_Git = glob.glob(param.indir)
    pattern = 'CFIS.V0.skycell.'
    ext = 'fits'
    hdu_no = 1

    output_base = 'UNIONS_'
    pattern_map = {
        'unconv' : '_image-',
        'wt' : '_weight-',
        'mask' : '_flag-'
    }


    ### Start main program ###

    for d in [param.outdir, param.outdir_orig]:
        if d and not os.path.isdir(d):
            raise IOError('Output path \'{}\' not a valid directory'
                        ''.format(d))

    if param.verbose:
        if param.last_Git:
            print('Converting PS image names in last Git run dir')
        else:
            print('Converting PS image names in all Git run dirs')

    for dir_Git in dirs_Git:

        if param.last_Git:
            # Only process last Git run
            if dir_Git != dirs_Git[-1]:
                continue

        dir_Git_out = '{}/get_images_runner/output'.format(dir_Git)
        dir_Git_out = os.path.abspath(dir_Git_out)

        if param.outdir:
            outdir = param.outdir
        else:
            outdir = dir_Git_out

        if param.verbose:
            print(dir_Git_out, outdir)
    
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
            if param.sp_format:
                # 000.000 -> 000-000
                number_final = in2out_pattern(number)
            else:
                number_final = number

            # Get image type
            mm = re.match('.*\.(.*)\.fits', input_name)
            if not mm[1] in pattern_map:
                print('Pattern \'{}\' in file name \'{}\' not matched, continuing...'
                      ''.format(mm[1], input_name))
                continue
            output_type = pattern_map[mm[1]]

            # Assemble output file name
            output_name = '{}{}{}{}.{}'.format(output_base, filter_letter,
                                               output_type, number_final, ext)

            output_path = '{}/{}'.format(outdir, output_name)

            ln_s(psfn, output_path, src_base=input_name, verbose=param.verbose)

            if param.outdir_orig:
                output_name = input_name
                output_path = '{}/{}'.format(param.outdir_orig, output_name)
                ln_s(psfn, output_path, psfn, verbose=param.verbose)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
