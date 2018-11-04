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
        input_dir_cat_tiles = '{}/tiles'.format(base_dir),
        input_dir_img_exp   = '{}/hdu'.format(base_dir),
        img_exp_pattern     = 'cfisexp-',
        output_dir_cat_exp  = 'out_cat_exp',
        log_path            = '{}/log_exposure.txt'.format(base_dir),
        cat_tiles_pattern   = 'CFIS-',
        cat_exp_pattern     = 'cfisexp-obj-',
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



def get_log_file(path, verbose=False):
    """Return log file content

    Parameters
    ----------
    path: string
        log file path
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    log: list of strings
        log file lines
    """

    if not os.path.isfile(path):
        stuff.error('Log file \'{}\' not found'.format(path))

    f_log = open(path, 'r')
    log   = f_log.readlines()
    if verbose:
        print('Reading log file, {} lines found'.format(len(log)))
    f_log.close()

    return log



def write_exposure_files(cat_tiles, log, cat_tiles_pattern, input_dir_img_exp, img_exp_pattern, output_dir_cat_exp, cat_exp_pattern, verbose=False):
    """Write catalogues corresponding to exposure coordinates with object info from corresponding tiles.

    Parameters
    ----------
    cat_tiles: list of string
        tiles catalogue file names
    log: list of string
        log file lines
    cat_tiles_pattern: string
        base tiles catalogue file name
    input_dir_img_exp: string
        input directory for expoure images
    img_exp_pattern:
        input exposure image file name base
    output_dir_cat_exp: string
        output directory for exposure catalogues
    cat_exp_base:
        output exposure catalogue file name base
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    cols  = ('X_IMAGE', 'Y_IMAGE', 'X_WORLD', 'Y_WORLD')
    dtype = [(c, float) for c in cols]
    dt    = [float for c in cols]

    exp_wcs = {}
    exp_cat = {}

    # First loop over tiles: Initialise exposure catalogues
    if verbose:
        print('First loop over tiles')
    for tile in cat_tiles:

        if verbose:
            print('tile {}'.format(tile))

        tile_num = stuff.get_pipe_file_number(cat_tiles_pattern, tile)

        # Get all exposure numbers for this tile from log file
        exp_num_list = cfis.log_get_exp_nums_for_tiles_num(log, tile_num)

        # If not already done (for previous tile): Get WCS header
        for exp_num in exp_num_list:
            if not exp_num in exp_wcs:
                exp_img_name = '{}/{}{}-0.fits'.format(input_dir_img_exp, img_exp_pattern, exp_num)
                exp_img      = sc.FITSCatalog(exp_img_name, hdu_no=0)
                exp_img.open()
                header        = exp_img.get_header(hdu_no=0)
                ny, nx = exp_img.get_data().shape
                exp_img.close()

                exp_wcs[exp_num] = wcs.WCS(header)

    if verbose:
        print('{} exposures found'.format(len(exp_wcs)))

    # Second loop over tiles: Distribute objects on tiles to exposure catalogues
    if verbose:
        print('Second loop over tiles')
    for tile in cat_tiles:

        if verbose:
            print('tile {}'.format(tile))

        # Open catalogue and get data
        f_tile = sc.FITSCatalog(tile, SEx_catalog=True)
        f_tile.open()
        tmp   = f_tile.get_data()

        # Use only columns given above
        dat_tile = Table([tmp[:][c] for c in cols], names=cols, dtype=dt)
        size  = len(dat_tile)

        tile_num     = stuff.get_pipe_file_number(cat_tiles_pattern, tile)
        exp_num_list = cfis.log_get_exp_nums_for_tiles_num(log, tile_num)

        # For all objects find exposure, set corresponding mask to True
        all_coord_tile_wcs = SkyCoord(ra=dat_tile['X_WORLD']*u.degree, dec=dat_tile['Y_WORLD']*u.degree)
        for exp_num in exp_num_list:
            all_coord_tile_xy  = exp_wcs[exp_num].all_world2pix(all_coord_tile_wcs.ra, all_coord_tile_wcs.dec, 0)
            ind_in_range       = ((all_coord_tile_xy[0] >= 0) & (all_coord_tile_xy[0] < nx) & \
                                  (all_coord_tile_xy[1] >= 0) & (all_coord_tile_xy[1] < ny))
            if ind_in_range.any():
                # Append objects within range to exposure catalogue
                if exp_num in exp_cat:
                    for d in dat_tile[ind_in_range]:
                        exp_cat[exp_num].add_row(d)
                else:
                    exp_cat[exp_num] = dat_tile[ind_in_range]

    # Write exposure catalogues to disk
    for exp_num in exp_cat:
        output_path = '{}/{}{}-0.fits'.format(output_dir_cat_exp, cat_exp_pattern, exp_num)
        print(output_path)
        exp_cat_file = sc.FITSCatalog(output_path, open_mode=sc.BaseCatalog.OpenMode.ReadWrite)
        #exp_cat_file.save_as_fits(data=np.array(exp_cat[exp_num]), names=cols, ext_name='LDAC_OBJECTS')
        exp_cat_file.save_as_fits(data=exp_cat[exp_num], names=cols, ext_name='LDAC_OBJECTS')

    if verbose:
        print('{} object files written'.format(len(exp_cat)))


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
    cat_tiles = stuff.get_file_list(param.input_dir_cat_tiles, param.cat_tiles_pattern, ext='.cat', verbose=param.verbose)

    # The log file lists all exposures for each tile
    log = get_log_file(param.log_path, verbose=param.verbose)

    write_exposure_files(cat_tiles, log, param.cat_tiles_pattern, param.input_dir_img_exp, param.img_exp_pattern, \
                         param.output_dir_cat_exp, param.cat_exp_pattern, verbose=param.verbose)

    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

