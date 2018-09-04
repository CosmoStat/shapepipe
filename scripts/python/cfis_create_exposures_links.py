#!/usr/bin/env python
  
"""Script cfis_create_exposures_links.py

Create links to exposures that are used in stacks.

:Authors: Martin Kilbinger

:Date: 4/09/2018
"""


# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import re
import copy
import glob

import numpy as np

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup
from astropy.io import fits

import stuff
import cfis


def params_default():
    """Set default parameter values.

    Parameters
    ----------
    None

    Returns
    -------
    p_def: class param
        parameter values
    """

    input_dir = '.'

    p_def = stuff.param(
        input_dir_tiles = input_dir,
        input_dir_exp   = '{}/astro/data/CFIS/pitcairn'.format(os.environ['HOME']),
        input_dir_exp_weights = '{}/astro/data/CFIS/weights'.format(os.environ['HOME']),
        output_dir      = '{}/exposures'.format(input_dir),
        band            = 'r',
        pattern_base    = 'CFIS-',
        exp_base_new    = 'CFISexp',
        exp_weight_base_new = 'CFISexp.weight',
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

    parser.add_option('-i', '--input_dir_tiles', dest='input_dir_tiles', type='string', default=p_def.input_dir_tiles,
         help='input directory for tiles, default=\'{}\''.format(p_def.input_dir_tiles))
    parser.add_option('-I', '--input_dir_exp', dest='input_dir_exp', type='string', default=p_def.input_dir_exp,
         help='input directory for exposures, default=\'{}\''.format(p_def.input_dir_exp))
    parser.add_option('', '--input_dir_exp_weights', dest='input_dir_exp_weights', type='string', default=p_def.input_dir_exp_weights,
         help='input directory for exposure weight maps, default=\'{}\''.format(p_def.input_dir_exp_weights))
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, where links will be created, default=\'{}\''.format(p_def.output_dir))

    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern_base,
        help='file pattern to match, default=\'{}\''.format(p_def.pattern_base))
    parser.add_option('', '--exp_base_new', dest='exp_base_new', type='string', default=p_def.exp_base_new,
         help='exposure base name of links to be created, default=\'{}\''.format(p_def.exp_base_new))
    parser.add_option('', '--exp_weight_base_new', dest='exp_weight_base_new', type='string', default=p_def.exp_weight_base_new,
         help='exposure weight map base name of links to be created, default=\'{}\''.format(p_def.exp_weight_base_new))


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

    return param



def get_tile_list(input_dir, pattern_base, band, verbose=False):
    """Return list of CFIS tiles filenames

    Parameters
    ----------
    input_dir: string
        input directory
    pattern_base: string
        base of file name pattern
    band: string
        band, one of 'r', 'u, only used if pattern_base=''
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    dst_list: list of string
        file name list
    """

    files = glob.glob('{}/*'.format(input_dir))

    #pattern = cfis.get_file_pattern(pattern_base, band, 'tile')
    pattern = pattern_base

    dst_list = []
    for f in files:

        # Test if file matches pattern
        m = re.findall(pattern, f)
        if len(m) != 0:
            dst_list.append(f)

    if len(dst_list) == 0:
        stuff.error('No files found in \'{}\' that matches pattern \'{}\''.format(input_dir, pattern))
    if verbose == True:
        print('Found {} tiles'.format(len(dst_list)))

    return dst_list



def get_exposure_list(tiles, verbose=False):
    """Return list of exposures that are used in the tiles stacks.

    Parameters
    ----------
    tiles: list of strings
        file names of tiles
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    exposures: list of strings
        file names of exposures
    """

    exp_list = []
    for f in tiles:

        try:
            hdu   = fits.open(f)
            hist  = hdu[0].header['HISTORY']

        except:

            if verbose:
                print('Error while reading FITS file {}, continuing...'.format(f))

        for h in hist:
            temp = h.split(' ')
            name = '{}.fz'.format(temp[3])
            exp_list.append(name)
    
    if verbose:
        print('Found {} exposures'.format(len(exp_list)))

    return exp_list



def create_links(exp_list, input_dir, input_dir_weights, output_dir, exp_base, exp_weight_base, verbose=False):
    """Create links to exposures in pipeline format.

    Parameters
    ----------
    exp_list: list of strings
        list of exposure file names
    input_dir: string
        input directory for exposures
    input_dir_weights: string
        input directory for exposure weight maps
    output_dir: string
        output directory
    exp_base: string
        base name of exposure link names in pipeline format
    exp_weight_base: string
        base name of exposure weight link names in pipeline format
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    if not os.path.isdir(output_dir):
        stuff.error('Path {} does not exist'.format(output_dir))

    num = 0
    ext = 'fits.fz'
    band = 'r'
    for exp in exp_list:

        print(exp)
        m = re.findall('(.*)\.{}'.format(ext), exp)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(exp))

        # Look for correponding weight image
        weight_name = cfis.get_file_pattern(m[0], band, 'exposure_weight.fz', want_re=False)
        weight_path = '{}/{}'.format(input_dir_weights, weight_name)
        if not os.path.isfile(weight_path):
            stuff.error('Weight file \'{}\' not found'.format(weight_path))

        # Link to image
        source = '{}/{}'.format(input_dir, exp)
        link   = '{}/{}-{:04d}-0.{}'.format(output_dir, exp_base, num, ext)
        os.symlink(source, link)
        print('symlink {} <- {}'.format(source, link))

        # Link to weight
        source = weight_path
        link   = '{}/{}-{:04d}-0.{}'.format(output_dir, exp_weight_base, num, ext)
        os.symlink(source, link)
        print('symlink {} <- {}'.format(source, link))

        num = num + 1

    if verbose:
        print('Create {} links'.format(num))



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    tiles     = get_tile_list(param.input_dir_tiles, param.pattern, param.band, verbose=param.verbose)

    exposures = get_exposure_list(tiles, verbose=param.verbose)
    create_links(exposures, param.input_dir_exp, param.input_dir_exp_weights, param.output_dir, param.exp_base_new, param.exp_weight_base_new, verbose=param.verbose)

    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

