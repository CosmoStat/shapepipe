#!/usr/bin/env python


"""Script create_image_links.py

Create links to images with link names according to pipeline format

:Authors: Martin Kilbinger

:Date: 15/02/2018
"""

# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import glob
import re
import copy

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup


from cfis import cfis
from generic import stuff



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
        input            = '.',
        output_dir       = '.',
        band             = 'r',
        tile_base_new    = 'CFIS',
        weight_base_new  = 'CFIS_weight',
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
    parser = OptionParser(usage=usage)

    parser.add_option('-i', '--input', dest='input', type='string', default=p_def.input,
         help='input image list, can be ascii file or directory path, default=\'{}\''.format(p_def.input))
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, where links will be created, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')

    parser.add_option('', '--tile_base_new', dest='tile_base_new', type='string', default=p_def.tile_base_new,
         help='tile base name of link to be created, default=\'{}\''.format(p_def.tile_base_new))
    parser.add_option('', '--weight_base_new', dest='weight_base_new', type='string', default=p_def.weight_base_new,
         help='weight base name of link to be created, default=\'{}\''.format(p_def.weight_base_new))

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



def create_links(inp, output_dir, tile_base_new, weight_base_new, band, verbose=False):

    image_list = cfis.get_image_list(inp, band, 'tile', col='#Name', verbose=verbose)
    file_list  = [i.name for i in image_list]

    if verbose:
        print('Found {} files'.format(len(file_list)))

    # Filter file list to match CFIS image pattern for tiles
    img_list = []
    pattern = cfis.get_file_pattern('', band, 'tile')
    for img in file_list:

        m = re.findall(pattern, img)
        if len(m) != 0:
            img_list.append(img)

    if verbose:
        print('Found {} CFIS images (tiles)'.format(len(img_list)))

    # Create links
    if verbose:
        print('Creating links:')

    num = 0
    ext = 'fits'
    for img in img_list:

        base_name = img
        m = re.findall('(.*)\..*', base_name)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(img))

        # Look for correponding weight image
        weight_name = cfis.get_file_pattern(m[0], band, 'weight', want_re=False)
        m = re.findall('(.*)\..*', weight_name)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(img))
        weight_base = m[0]

        # Check whether weight image file exists
        #weight_path = '{}/{}'.format(input_dir, weight_name)
        weight_path = '{}'.format(weight_name)
        if not os.path.isfile(weight_path):
            stuff.error('Weight file \'{}\' not found'.format(weight_path))

        link_name_tile = '{}/{}-{:03d}-0.{}'.format(output_dir, tile_base_new, num, ext)
        if verbose:
            print(' {} <- {}'.format(img, link_name_tile))
        os.symlink(img, link_name_tile)

        link_name_weight = '{}/{}-{:03d}-0.{}'.format(output_dir, weight_base_new, num, ext)
        if verbose:
            print(' {} <- {}'.format(weight_path, link_name_weight))
        os.symlink(weight_path, link_name_weight)

        num = num + 1



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

    create_links(param.input, param.output_dir, param.tile_base_new, param.weight_base_new, param.band, verbose=param.verbose)


if __name__ == "__main__":
    sys.exit(main(sys.argv))


