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

import cfis
import stuff



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
        input_dir  = '../tiles',
        band       = 'r',
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

    return p_def


def create_links(input_dir, band):

    # Read all file names
    file_list = glob.glob('{}/*'.format(input_dir))
    file_list.sort()

    # Filter file list to match CFIS image pattern for tiles
    img_list = []
    pattern = cfis.get_file_pattern('', band, 'tile')
    for img in file_list:

        m = re.findall(pattern, img)
        if len(m) != 0:
            img_list.append(img)

    # Create links
    print('Creating links:')
    num = 0
    ext = 'fits'
    for img in img_list:

        # Tile name stripped of path and extension
        base_name = os.path.basename(img)
        m = re.findall('(.*)\..*', base_name)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(img))
        tile_base = m[0]

        # Look for correponding weight image
        weight_name = cfis.get_file_pattern(m[0], band, 'weight', want_re=False)
        m = re.findall('(.*)\..*', weight_name)
        if len(m) == 0:
            stuff.error('Invalid file name \'{}\' found'.format(img))
        weight_base = m[0]

        # Check whether weight image file exists
        weight_path = '{}/{}'.format(input_dir, weight_name)
        if not os.path.isfile(weight_path):
            stuff.error('Weight file \'{}\' not found'.format(weight_path))

        link_name_tile = '{}-{:03d}-0.{}'.format(tile_base, num, ext)
        print(' {} <- {}'.format(img, link_name_tile))
        os.symlink(img, link_name_tile)

        link_name_weight = '{}-{:03d}-0.{}'.format(weight_base, num, ext)
        print(' {} <- {}'.format(weight_path, link_name_weight))
        os.symlink(weight_path, link_name_weight)

        num = num + 1



def main(argv=None):
    """Main program.
    """

    # Set default parameters
    p_def = params_default()

    param = p_def

    create_links(param.input_dir, param.band)


if __name__ == "__main__":
    sys.exit(main(sys.argv))


