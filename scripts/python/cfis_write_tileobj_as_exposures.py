#!/usr/bin/env python
  
"""Script cfis_write_tileobj_as_exposures.py

For objects in catalogue detected on tiles, write
files according to the exposures where the object appears.

:Authors: Martin Kilbinger

:Date: 12/10/2018
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


from generic import stuff
from cfis import cfis


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
        input_dir_cat_tiles    = input_dir,
        output_dir_cat_exp     = '{}/exposures'.format(input_dir),
        log_path               = '{}/log_exposures.txt'.format(input_dir),
        cat_tiles_pattern_base = 'CFIS-',
        cat_exp_pattern_base   = 'cfisexp-obj',
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

    # Input
    parser.add_option('-i', '--input_dir_cat_tiles', dest='input_dir_cat_tiles', type='string', default=p_def.input_dir_cat_tiles,
         help='input directory for tiles catalogues, default=\'{}\''.format(p_def.input_dir_cat_tiles))

    # Output
    parser.add_option('-o', '--output_dir_cat_exp', dest='output_dir_cat_exp', type='string', default=p_def.output_dir_cat_exp,
         help='output directory for exposure catalogues, default=\'{}\''.format(p_def.output_dir_cat_exp))
    parser.add_option('-l', '--log_path', dest='log_path', type='string', default=p_def.log_path,
         help='log file name, default=\'{}\''.format(p_def.log_path))

    parser.add_option('-p', '--cat_tiles_pattern', dest='cat_tiles_pattern', type='string', default=p_def.cat_tiles_pattern,
        help='file pattern to match input tiles catalogues, default=\'{}\''.format(p_def.cat_tiles_pattern))
    parser.add_option('-P', '--cat_exp_pattern', dest='cat_exp_tile_pattern', type='string', default=p_def.cat_exp_pattern,
         help='file pattern for output exposure catalogues, default=\'{}\''.format(p_def.cat_exp_pattern))

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

    if not options.output_format in ('links', 'hdu'):
        stuff.error('Option -O (--output_format) needs to be \'links\' or \'hdu\'')

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

    if param.exp_weight_base_new is None:
        param.exp_weight_base_new = '{}_weight'.format(param.exp_base_new)
    if param.exp_flag_base_new is None:
        param.exp_flag_base_new = '{}_flag'.format(param.exp_base_new)
        

    return param



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

    # Save calling command
    stuff.log_command(argv)
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')

    if param.verbose is True:
        print('Start of program {}'.format(os.path.basename(argv[0])))

    ### Start main program ###

    # Get list of catalogues of objects selected on tiles
    cat_tiles = get_cat_list(param.input_dir_cat_tiles, param.cat_tiles_pattern_base, verbose=param.verbose)

    # The log file lists all exposures for each tile
    log = get_log_file(param.log_path, verbose=param.verbose)

    write_exposure_files(cat_tiles, log, param.output_dir_cat_exp, param.cat_exp_base, verbose=param.verbose)

    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

