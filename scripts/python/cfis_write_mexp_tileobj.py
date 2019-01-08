#!/usr/bin/env python

"""Script cfis_write_mexp_tileobj.py

For objects in catalogue detected on tiles, write
new multi-exposure files according to the exposures where the object appears.

:Authors: Martin Kilbinger

:Date: 11/12/2018

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
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
from astropy.table import Table

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
        input_dir_cat_exp   = '{}/hdu'.format(base_dir),
        input_dir_psf       = '{}/hdu'.format(base_dir),
        log_path            = '{}/log_exposure.txt'.format(base_dir),
        cat_exp_pattern     = 'cfisexp-obj-',
        psf_pattern         = 'galaxy_psf-',
        output_dir          = 'outdir_gal_psf_mexp',
        outcat_pattern      = 'CFIS_MOBJ-',
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
    parser.add_option('', '--input_dir_cat_exp', dest='input_dir_cat_exp', type='string', default=p_def.input_dir_cat_exp,
         help='input directory for object exposure catalogues, default=\'{}\''.format(p_def.input_dir_cat_exp))
    parser.add_option('', '--input_dir_psf', dest='input_dir_psf', type='string', default=p_def.input_dir_psf,
         help='input directory for PSF files, default=\'{}\''.format(p_def.input_dir_psf))
    parser.add_option('-l', '--log_path', dest='log_path', type='string', default=p_def.log_path,
         help='log file name, default=\'{}\''.format(p_def.log_path))
    parser.add_option('-P', '--cat_exp_pattern', dest='cat_exp_tile_pattern', type='string', default=p_def.cat_exp_pattern,
         help='input file pattern for object exposure catalogues, default=\'{}\''.format(p_def.cat_exp_pattern))
    parser.add_option('-p', '--psf_pattern', dest='psf_pattern', type='string', default=p_def.psf_pattern,
         help='input psf file pattern, default=\'{}\''.format(p_def.psf_pattern))

    # Output
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-O', '--outcat_pattern', dest='outcat_pattern', type='string', default=p_def.outcat_pattern,
         help='output catalogue file pattern, default=\'{}\''.format(p_def.outcat_pattern))
    parser.add_option('', '--vignet', dest='vignet', action='store_true', help='Add PSF vignets to output files')

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


 
def collect_info_gal_psf(log, input_dir_cat_exp, cat_exp_pattern, input_dir_psf, psf_pattern, output_dir, outcat_pattern,
                         verbose=False, vignet=False):
    """Collect information about galaxies and PSF from exposures used in tiles
       given in log file.
    Parameters
    ----------
    log: list of string
        log file lines
    input_dir_cat_exp: string
        input directory for exposure catalogues
    cat_exp_pattern: string
        file name pattern for exposure catalogues
    input_dir_psf: string
        input directory for PSF files
    psf_pattern: string
        file name pattern for PSF files
    outut_dir: string
        output directory
    outcat_pattern: string
        name pattern for output files
    verbose: bool, optional, default=False
        verbose output if True
     vignet: bool, optional, default=False
        if True, adds PSF vignet to output files

    Returns
    -------
    None
    """

    tile_num_list = cfis.log_get_tile_nums(log)

    # loop over tiles
    for tile_num in tile_num_list:

        tile_data = None

        if verbose:
            print('tile {}'.format(tile_num))

        # Get all exposure numbers for this tile from log file
        exp_num_list = cfis.log_get_exp_nums_for_tiles_num(log, tile_num)

        for exp_num in exp_num_list:

            # Object exposure catalogue
            exp_cat_name = '{}/{}{:03d}-0.fits'.format(input_dir_cat_exp, cat_exp_pattern, int(exp_num))
            hdu_no       = 2
            exp_cat      = sc.FITSCatalog(exp_cat_name, hdu_no=hdu_no)
            try:
                exp_cat.open()
            except:
                print('Object exposure catalogue \'{}\' not found, continuing...'.format(exp_cat_name))
                continue

            print('Found \'{}\''.format(exp_cat_name))
            exp_cat_data = exp_cat.get_data()
            n_data = len(exp_cat_data)

            # PSF file
            psf_name = '{}/{}{:03d}-0.fits'.format(input_dir_psf, psf_pattern, int(exp_num))
            psf      = sc.FITSCatalog(psf_name, hdu_no=2)
            try:
                psf.open()
            except:
                print('PSF file \'{}\' not found, please check later...'.format(psf_name))
                continue

            psf_data = psf.get_data()

            if n_data != len(psf_data):
                raise stuff.MyError('Lengh of objects ({}) and PSF ({}) have to be the same'.format(n_data, len(psf_data)))

            # Filter objects that come from the current tile
            ind_tile = exp_cat_data['tile_num'] == int(tile_num)
            if ind_tile.any():

                exp_cat_data_plus = {}

                # Copy object data
                for c in exp_cat.get_col_names():
                    exp_cat_data_plus[c] = exp_cat_data[c][ind_tile]

                # Add PSF information to object catalogue
                for c in psf.get_col_names():
                    if vignet or c != 'VIGNET':
                        # Add vignet only if argument 'vignet' is True
                        exp_cat_data_plus[c] = psf_data[c][ind_tile]

                if tile_data is None:
                    tile_data = exp_cat_data_plus
                else:
                    for c in exp_cat_data_plus:
                        tile_data[c] = np.concatenate((exp_cat_data_plus[c], tile_data[c]))

            exp_cat.close()
            psf.close()

        #import ipdb; ipdb.set_trace()

        # TODO (here?) selecting objects according to size

        # Saving to file maybe only for testing. Note that tile_num here is string
        out_path = '{}/{}{}-0.fits'.format(output_dir, outcat_pattern, tile_num)
        print('Writing tile data to file \'{}\''.format(out_path))
        if os.path.exists(out_path):
            os.unlink(out_path)
        output = sc.FITSCatalog(out_path, open_mode=sc.BaseCatalog.OpenMode.ReadWrite)
        output.save_as_fits(tile_data)



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

    # The log file lists all exposures for each tile
    log = cfis.get_log_file(param.log_path, verbose=param.verbose)
 
    collect_info_gal_psf(log, param.input_dir_cat_exp, param.cat_exp_pattern, param.input_dir_psf, param.psf_pattern, param.output_dir, param.outcat_pattern,
                         verbose=param.verbose, vignet=param.vignet)


    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))
