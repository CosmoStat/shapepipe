#!/usr/bin/env python

"""Script cfis_select_mexp.py

Select objects according to their multi-exposure data.

:Authors: Martin Kilbinger

:Date: 07/01/2019

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

import scatalog as sc

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

    base_dir = '{}/data'.format(os.environ['HOME'])

    p_def = stuff.param(
        input_dir_cat_mexp  = '{}/tiles'.format(base_dir),
        mexp_pattern        = 'CFIS_MOBJ-',
        method              = 'PSF_size',
        param_method        = '',
        output_dir          = '{}/tiles'.format(base_dir),
        outcat_pattern      = 'CFIS_GAL-',
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
    parser = OptionParser(usage=usage, formatter=stuff.IndentedHelpFormatterWithNL())

    # Input
    parser.add_option('', '--input_dir_cat_mexp', dest='input_dir_cat_mexp', type='string', default=p_def.input_dir_cat_mexp,
         help='input directory for multi-exposure catalogues, default=\'{}\''.format(p_def.input_dir_cat_mexp))
    parser.add_option('-p', '--mexp_pattern', dest='mexp_pattern', type='string', default=p_def.mexp_pattern,
         help='input multi-exposure file pattern, default=\'{}\''.format(p_def.mexp_pattern))
    parser.add_option('-m', '--method', dest='method', type='string', default=p_def.method,
         help='selection method, default=\'{}\''.format(p_def.method))
    parser.add_option('-P', '--param', dest='param_method', type='string', default=p_def.param_method,
         help='selection method parameters'.format(p_def.param_method))

    # Output
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-O', '--outcat_pattern', dest='outcat_pattern', type='string', default=p_def.outcat_pattern,
         help='output catalogue file pattern, default=\'{}\''.format(p_def.outcat_pattern))

    parser.add_option('-c', '--check_consistency', dest='check_consistency', action='store_true', help='check_consistency')
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



def check_consistency(gal, nexp):
    """Checks consistency of multi-exposure galaxy data

    Parameters
    ----------
    gal: list of FITS_rec
        galaxy and PSF data from multiple exposures
    nexp: int
        number of exposures

    Returns
    -------
    None
    """

    if len(gal) != nexp:
        raise stuff.MyError('Length of galaxy data {} != number of exposures {}'.format(len(gal), nexp))



def select_one(gal, nexp, method, params):
    """Select galaxy according to multi-exposure information.

    Parameters
    ----------
    gal: list of FITS_rec
        galaxy and PSF data from multiple exposures
    nexp: int
        number of exposures
    method: string
        selection method name
    params: dict
        method parameters

    Returns
    -------
    is_gal, params_out: bool
        True if object is selected as galaxy
    """

    is_gal    = False
    params_out = {}

    if method == 'PSF_size':

        fwhm = gal['FWHM_IMAGE'][0]
        n_larger = 0
        for g in gal:
            if g['HSM_FLAG'] == 0 and fwhm > 2.355 * g['SIGMA_PSF_HSM']:
                n_larger += 1
        params_out['n_larger'] = n_larger
        if n_larger == nexp:
            is_gal = True

    else:

        raise stuff.MyError('Invalid selection method \'{}\''.format(method))

    return is_gal, params_out



def select_all(cat_tiles, method, params_in, do_check_consistency=False, verbose=False):
    """Select multi-exposure objects.

    Parameters
    ----------
    cat_tiles: list of string
        file names of mexp tile files
    method: string
        selection method name
    params_in: dict
        method parameters
    do_check_consistency: bool, optional, default=False
        check consistency of input galaxy multi-exposure catalogue if True
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    for tile in cat_tiles:

        if verbose:
            print('tile {}'.format(tile))

        # Open catalogue and get data
        f_tile = sc.FITSCatalog(tile, SEx_catalog=False)
        f_tile.open()
        dat = f_tile.get_data()

        col_names = f_tile.get_col_names()

        # Create empty data array for selected galaxies
        dat_gal = {}
        for c in col_names: 
            dat_gal[c] = []
        dat_gal['nexp'] = []

        IDs = dat['ID']
        IDs_unique, nexp = np.unique(IDs, return_counts=True)

        # MKDEBUG
        n_is_gal = 0
        for u, n in zip(IDs_unique, nexp):
            gal = dat[IDs == u]

            if do_check_consistency:
                check_consistency(gal, n)
            is_gal, params_out = select_one(gal, n, method, params_in)
            if is_gal:
                n_is_gal += 1
                for c in col_names: 
                    dat_gal[c].append(gal[0][c])
                dat_gal['nexp'].append(n)

        print(n_is_gal)
        if verbose:
            print('{}/{} galaxies  selected'.format(len(dat_gal['nexp']), len(dat)))



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
    cat_tiles = stuff.get_file_list(param.input_dir_cat_mexp, param.mexp_pattern, ext='.fits', verbose=param.verbose)

    select_all(cat_tiles, param.method, param.param_method, do_check_consistency=param.check_consistency, verbose=param.verbose)


    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))
