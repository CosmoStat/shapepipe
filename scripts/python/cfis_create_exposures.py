#!/usr/bin/env python
  
"""Script cfis_create_exposures.py

For exposures that are used in stacks, create either:
  (1) links to the exposure files
  (2) FITS files for each HDU

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
        input_dir_tiles = input_dir,
        input_dir_exp   = '{}/astro/data/CFIS/pitcairn'.format(os.environ['HOME']),
        input_dir_exp_weights = '{}/astro/data/CFIS/weights'.format(os.environ['HOME']),
        input_dir_exp_flags = '{}/astro/data/CFIS/flags'.format(os.environ['HOME']),
        output_dir      = '{}/exposures'.format(input_dir),
        output_format   = 'hdu',
        log_path        = '{}/log_exposures.txt'.format(input_dir),
        band            = 'r',
        pattern_base    = 'CFIS-',
        exp_base_new    = 'CFISexp',
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
    parser.add_option('-i', '--input_dir_tiles', dest='input_dir_tiles', type='string', default=p_def.input_dir_tiles,
         help='input directory for tiles, default=\'{}\''.format(p_def.input_dir_tiles))
    parser.add_option('-I', '--input_dir_exp', dest='input_dir_exp', type='string', default=p_def.input_dir_exp,
         help='input directory for exposures, default=\'{}\''.format(p_def.input_dir_exp))
    parser.add_option('', '--input_dir_exp_weights', dest='input_dir_exp_weights', type='string', default=p_def.input_dir_exp_weights,
         help='input directory for exposure weight maps, default=\'{}\''.format(p_def.input_dir_exp_weights))
    parser.add_option('', '--input_dir_exp_flags', dest='input_dir_exp_flags', type='string', default=p_def.input_dir_exp_flags,
         help='input directory for exposure flag maps, default=\'{}\''.format(p_def.input_dir_exp_flags))

    # Output
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, where links will be created, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-O', '--output_format', dest='output_format', type='string', default=p_def.output_format,
         help='output format, one in \'links\' (create links), \'hdu\' (write FITS files for each HDU; default)')
    parser.add_option('-l', '--log_path', dest='log_path', type='string', default=p_def.log_path,
         help='log file name, default=\'{}\''.format(p_def.log_path))

    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-p', '--pattern', dest='pattern', type='string', default=p_def.pattern_base,
        help='file pattern to match, default=\'{}\''.format(p_def.pattern_base))
    parser.add_option('', '--exp_base_new', dest='exp_base_new', type='string', default=p_def.exp_base_new,
         help='exposure base name of links to be created, default=\'{}\''.format(p_def.exp_base_new))
    parser.add_option('', '--exp_weight_base_new', dest='exp_weight_base_new', type='string',
         help='exposure weight map base name of links to be created, default=\'<EXP_BASE_NEW>_weight\'')
    parser.add_option('', '--exp_flag_base_new', dest='exp_flag_base_new', type='string',
         help='exposure flag map base name of links to be created, default=\'<EXP_BASE_NEW>_flag\'')

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



def get_exposure_list(tiles, pattern_base, verbose=False):
    """Return list of exposures that are used in the tiles stacks.

    Parameters
    ----------
    tiles: list of strings
        file names of tiles
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    exposures: list of tuples of strings
        tupel of file names of exposures, tile number
    """

    exp_list = []
    for f in tiles:

        try:
            hdu   = fits.open(f)
            hist  = hdu[0].header['HISTORY']

        except:
            if verbose:
                print('Error while reading FITS file {}, continuing...'.format(f))

        tile_num = stuff.get_pipe_file_number(pattern_base, f)

        # Get exposure file names
        for h in hist:
            temp     = h.split(' ')
            exp_name = '{}.fz'.format(temp[3])
            exp_list.append((exp_name, tile_num))

    exp_list_uniq = stuff.list_unique(exp_list)
    
    if verbose:
        print('Found {} exposures used in {} tiles'.format(len(exp_list_uniq), len(tiles)))
        print('{} duplicates were removed'.format(len(exp_list) - len(exp_list_uniq)))

    return exp_list_uniq



def create_links(num, output_dir, exp_path, weight_path, flag_path, exp_base, \
                 exp_weight_base, exp_flag_base, ext, verbose=False):
    """Create links image, weight, and flag file
    """

    source = exp_path
    link   = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_base, num, ext)
    os.symlink(source, link)
    if verbose:
        print('symlink {} <- {}'.format(source, link))

    # Link to weight
    source = weight_path
    link   = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_weight_base, num, ext)
    os.symlink(source, link)
    if verbose:
        print('symlink {} <- {}'.format(source, link))

    # Link to flag
    source = flag_path
    link   = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_flag_base, num, ext)
    os.symlink(source, link)
    if verbose:
        print('symlink {} <- {}'.format(source, link))

    return num + 1

    

def get_n_hdu_from_log(path, log):

    n_hdu = 0
    for line in log:
        if path in line:
            n_hdu = n_hdu + 1

    return n_hdu



def write_hdu(k_img, k_weight, k_flag, img_file, weight_file, flag_file, output_dir, exp_base, exp_weight_base, exp_flag_base, \
              ext, num, exp_path, verbose=False):
    """Write HDUs.
    """

    # Change coordinates to astropy-readable format
    import sip_tpv as stp

    h = img_file._cat_data[k_img].header
    stp.pv_to_sip(h)

    d = img_file._cat_data[k_img].data
    new_fits = fits.PrimaryHDU(data=d, header=h)

    from astropy import wcs
    
    out_name = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_base, num, ext)
    if not os.path.isfile(out_name):
        new_fits.writeto(out_name)
        img_str = 'written'
    else:
        # TODO: error, or earlier check whether all images in log file are written to disk
        img_str = 'skipped'

    h_weight = weight_file._cat_data[k_weight].header
    d_weight = weight_file._cat_data[k_weight].data
    new_fits = fits.PrimaryHDU(data=d_weight, header=h_weight)
    out_name = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_weight_base, num, ext)
    if not os.path.isfile(out_name):
        new_fits.writeto(out_name)
        weight_str = 'written'
    else:
        # TODO: error, or earlier check whether all images in log file are written to disk
        weight_str = 'skipped'

    h_flag = flag_file._cat_data[k_flag].header
    d_flag = flag_file._cat_data[k_flag].data
    d_flag = d_flag.astype(np.int16)
    new_fits = fits.PrimaryHDU(data=d_flag, header=h_flag)
    out_name = '{}/{}-{:03d}-0.{}'.format(output_dir, exp_flag_base, num, ext)
    if not os.path.isfile(out_name):
        new_fits.writeto(out_name)
        flag_str = 'written'
    else:
        # TODO: error, or earlier check whether all images in log file are written to disk
        flag_str = 'skipped'

    if verbose:
         print('Image/weight/flag file name/number {}/{} hdu {}/{}/{} {}/{}/{}'.\
                format(exp_path, num, k_img, k_weight, k_flag, img_str, weight_str, flag_str))




def create_hdus(num, output_dir, exp_path, weight_path, flag_path, \
                exp_base, exp_weight_base, exp_flag_base, ext, tile_num, log, verbose=False):
    """Create FITS files for each hdu in the image (CCDs)
    """

    import scatalog as sc

    log_app = []

    n_hdu = get_n_hdu_from_log(exp_path, log)
    if n_hdu != 0:
        if verbose:
            print('Skipping image {}, found {} HDUs in log file'.format(exp_path, n_hdu))
        num = num + n_hdu
        return num, log_app

    img_file = sc.FITSCatalog(exp_path, hdu_no=1)
    img_file.open()
    weight_file = sc.FITSCatalog(weight_path, hdu_no=1)
    weight_file.open()
    flag_file = sc.FITSCatalog(flag_path, hdu_no=1)
    flag_file.open()
    
    hdu_max = [len(cat._cat_data) for cat in (img_file, weight_file, flag_file)]
    all_hdu = True
    if not all(np.diff(hdu_max) == 0):
        all_hdu = False
        warning.warn('Inconsistent #hdus={}/{}/{} of image/weight/flag, writing only some hdus'.format(hdu_max[0], hdu_max[1], hdu_max[2]))


    # Write FITS files
    for k_img in range(1, 41):

        if all_hdu:
            k_weight = k_img
            k_flag   = k_img
        else:
            h_img = img_file._cat_data[k_img].header
            coord_img = re.findall(r"[\w]+", h_img['DETSEC'])

            k_weight_match = -1
            for k_weight in range(1, hdu_max[1]):
                h_weight = img_file._cat_data[k_weight].header
                coord_weight = re.findall(r"[\w]+", h_weight['DETSEC'])
                if coord_img == coord_weight:
                    k_weight_match = k_weight

            k_flag_match = -1
            for k_flag in range(1, hdu_max[2]):
                h_flag = img_file._cat_data[k_flag].header
                coord_flag = re.findall(r"[\w]+", h_flag['DETSEC'])
                if coord_img == coord_flag:
                    k_flag_match = k_flag
                    break

            if k_weight_match == -1 and k_flag_match == -1:
                print('No matching weight and flag HDUs found for image hdu {}'.format(k_img))
                continue

            k_weight = k_weight_match
            k_flag   = k_flag_match
             
        write_hdu(k_img, k_weight, k_flag, img_file, weight_file, flag_file, output_dir, exp_base, exp_weight_base, exp_flag_base, \
                  ext, num, exp_path, verbose=verbose)

        log_app = cfis.log_append_to_tiles_exp(log_app, exp_path, tile_num, k_img, k_weight, k_flag, num)

        num = num + 1

    img_file.close()
    weight_file.close()
    flag_file.close()

    return num, log_app



def create_output(exp_list, input_dir, input_dir_weights, input_dir_flags, output_dir, \
                 exp_base, exp_weight_base, exp_flag_base, log_path, output_format='hdu', verbose=False):
    """Create links/write FITS files for exposures in pipeline format.

    Parameters
    ----------
    exp_list: list of tupels
        list of tupels with exposure file names, tile numbers
    input_dir: string
        input directory for exposures
    input_dir_weights: string
        input directory for exposure weight maps
    input_dir_flags: string
        input directory for exposure flag maps
    output_dir: string
        output directory
    exp_base: string
        base name of exposure link names in pipeline format
    exp_weight_base: string
        base name of exposure weight link names in pipeline format
    exp_flag_base: string
        base name of exposure flag link names in pipeline format
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    None
    """

    if not os.path.isdir(output_dir):
        stuff.error('Path {} does not exist'.format(output_dir))

    if not os.path.isfile(log_path):
        # Create empty log file
        open(log_path, 'a').close()

    # Open log file for append-write
    f_log = open(log_path, 'a+')
    log   = f_log.readlines()
    if verbose:
        print('Reading log file, {} lines found'.format(len(log)))


    num = 0
    ext = 'fits'
    band = 'r'
    for (exp, tile_num) in exp_list:

        m = re.findall('(.*)\.{}'.format(ext), exp)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(exp))

        exp_path = '{}/{}'.format(input_dir, exp)

        # Look for correponding weight image
        weight_name = cfis.get_file_pattern(m[0], band, 'exposure_weight.fz', want_re=False)
        weight_path = '{}/{}'.format(input_dir_weights, weight_name)
        if not os.path.isfile(weight_path):
            stuff.error('Weight file \'{}\' not found'.format(weight_path))

        # Look for correponding flag image
        flag_name = cfis.get_file_pattern(m[0], band, 'exposure_flag.fz', want_re=False)
        flag_path = '{}/{}'.format(input_dir_flags, flag_name)
        if not os.path.isfile(flag_path):
            stuff.error('Flag file \'{}\' not found'.format(flag_path))

        if output_format == 'links':
            num = create_links(num, output_dir, exp_path, weight_path, flag_path, \
                               exp_base, exp_weight_base, exp_flag_base, ext, verbose=verbose)
        else:
            num, log_app = create_hdus(num, output_dir, exp_path, weight_path, flag_path, \
                                       exp_base, exp_weight_base, exp_flag_base, ext, tile_num, log, verbose=verbose)
            # Append newly written file info to log file
            if len(log_app) > 0:
                f_log.writelines(log_app)
                f_log.flush()

    if verbose:
        print('Created {} links'.format(num))

    f_log.close()


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

    tiles     = stuff.get_file_list(param.input_dir_tiles, param.pattern, verbose=param.verbose)

    exposures = get_exposure_list(tiles, param.pattern, verbose=param.verbose)

    print('MKDEBUG TODO: more than 3 digits for file number')
    create_output(exposures, param.input_dir_exp, param.input_dir_exp_weights, param.input_dir_exp_flags, param.output_dir, \
                  param.exp_base_new, param.exp_weight_base_new, param.exp_flag_base_new, param.log_path, output_format=param.output_format, verbose=param.verbose)

    ### End main program ###

    if param.verbose is True:
        print('End of program {}'.format(os.path.basename(argv[0])))


    return 0



if __name__ == "__main__":
    sys.exit(main(sys.argv))

