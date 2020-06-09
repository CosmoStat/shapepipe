#!/usr/bin/env python


"""Script create_image_links.py

Create links to images with link names according to pipeline format

:Authors: Martin Kilbinger

:Date: 15/02/2018 (v 1.0)
       03/03/2020 (v 1.1: added links for exposures)

:Package: ShapePipe
"""

# Compability with python2.x for x>6
from __future__ import print_function


import sys
import os
import glob
import re
import copy

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup

import cfis


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

    p_def = cfis.param(
        input = '.',
        output_dir = '.',
        band = 'r',
        image_type = 'exposure',
        image_base_new = 'image',
        weight_base_new = 'weight',
        flag_base_new = 'flag',
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

    parser.add_option('-i', '--input', dest='input', type='string', default=p_def.input,
         help='input image list, can be ascii file or directory path, default=\'{}\''.format(p_def.input))
    parser.add_option('-o', '--output_dir', dest='output_dir', type='string', default=p_def.output_dir,
         help='output directory, where links will be created, default=\'{}\''.format(p_def.output_dir))
    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='image_type', type='string', default=p_def.image_type,
        help='input image type, here: one of \'exposure\' (default)|\'tile\'')

    parser.add_option('', '--image_base_new', dest='image_base_new', type='string', default=p_def.image_base_new,
         help='image file base name of link to be created, default=\'{}\''.format(p_def.image_base_new))
    parser.add_option('', '--weight_base_new', dest='weight_base_new', type='string', default=p_def.weight_base_new,
         help='weight file base name of link to be created, default=\'{}\''.format(p_def.weight_base_new))
    parser.add_option('', '--flag_base_new', dest='flag_base_new', type='string', default=p_def.flag_base_new,
         help='flag file base name of link to be created, default=\'{}\''.format(p_def.flag_base_new))

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

    if not options.image_type in ('exposure', 'tile'):
        print('image type (option \'-t\') must be one of \'exposure\'|\'tile\'')
        return False

    if not os.path.isdir(options.output_dir):
        print('Output directory \'{}\' does not exist'.format(options.output_dir))
        return False

    return True


def update_param(p_def, options):
    """Return default parameter, updated and complemented according to options.
    
    Parameters
    ----------
    p_def:  class cfis.param
        parameter values
    optiosn: tuple
        command line options
    
    Returns
    -------
    param: class cfis.param
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


def create_links(inp, output_dir, image_type, image_base_new, weight_base_new, flag_base_new, band, verbose=False):

    image_list = cfis.get_image_list(inp, band, image_type, verbose=verbose)
    file_list  = [i.name for i in image_list]

    if verbose:
        print('Found {} files'.format(len(file_list)))

    # Filter file list to match CFIS image pattern
    img_list = []
    pattern = cfis.get_file_pattern('', band, image_type)
    for img in file_list:

        m = re.findall(pattern, img)
        if len(m) != 0:
            img_list.append(img)

    if verbose:
        print('Found {} CFIS images (type {})'.format(len(img_list), image_type))

    # Create links
    if verbose:
        print('Creating links:')

    num = 0

    if image_type == 'exposure':
        ext = 'fitsfz'
        weight_type = 'exposure_weight.fz'
        ext_weight = ext
        flag_type = 'exposure_flag.fz'
        num_form = ''
    elif image_type == 'tile':
        ext = 'fits'
        weight_type = 'weight.fz'
        ext_weight = 'fits.fz'
        flag_type = None
        num_form = ':04d'

    for img in img_list:

        base_name = img
        m = re.findall('(.*)\.fits.?', base_name)
        if len(m) == 0:
            raise cfis.CfisError('Invalid file name \'{}\' found'.format(img))

        # File number: Use running number fo tile, and image digits for exposures
        if image_type == 'tile':
            num_str = num
        else:
            mm = re.findall('.*/(\d*).*', m[0])
            if len(mm) == 0:
                raise cfis.CfisError('Invalid file name \'{}\' found'.format(mm))
            num_str = mm[0]

        # Look for correponding weight image
        weight_name = cfis.get_file_pattern(m[0], band, weight_type, want_re=False)
        mw = re.findall('(.*)\..*', weight_name)
        if len(mw) == 0:
            raise cfis.CfisError('Invalid file name \'{}\' found'.format(img))

        # Check whether weight image file exists
        weight_path = '{}'.format(weight_name)
        if not os.path.isfile(weight_path):
            raise cfis.CfisError('Weight file \'{}\' not found'.format(weight_path))

        # Create links
        form = '{{}}/{{}}-{{{}}}.{{}}'.format(num_form)
        link_name_image = form.format(output_dir, image_base_new, num_str, ext)
        cfis.symlink(img, link_name_image, verbose=verbose)

        link_name_weight = form.format(output_dir, weight_base_new, num_str, ext_weight)
        cfis.symlink(weight_path, link_name_weight, verbose=verbose)

        # Create link for flag if it exists
        if flag_type: 
            flag_name = cfis.get_file_pattern(m[0], band, flag_type, want_re=False)
            m = re.findall('(.*)\..*', flag_name)
            if len(m) > 0:
                flag_path = '{}'.format(flag_name)
                if os.path.isfile(flag_path):
                    link_name_flag = form.format(output_dir, flag_base_new, num_str, ext)
                    cfis.symlink(flag_path, link_name_flag, verbose=verbose)

        num = num + 1

    if verbose:
        print('Links for {} images (+weights, +flags) created'.format(num))


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
    cfis.log_command(argv)
    if param.verbose:
        cfis.log_command(argv, name='sys.stdout')

    create_links(param.input, param.output_dir, param.image_type, param.image_base_new,
                 param.weight_base_new, param.flag_base_new,
                 param.band, verbose=param.verbose)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))


