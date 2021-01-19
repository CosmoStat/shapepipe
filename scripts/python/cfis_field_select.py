#!/usr/bin/env python


"""Script cfis_field_select.py

Handling and selecting CFIS fields and pointings.

:Authors: Martin Kilbinger

:Date: 19/01/2018
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

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup

from astropy.io import fits
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle, SkyCoord

import cfis


def get_images_used_in_tiles(images, band, image_type):
    """Return exposures used in tiles.

    Parameters
    ----------
    images: list of class cfis.image
        list of images
    band: string
        optical band
    image_type: string
        image type ('exposure', 'exposure_weight.fz', \
        'exposure_flag', 'exposure_flag.fz', 'cat')
    """

    exp_list = []
    for img in images:

        try:
            hdu = fits.open(img.name)
            hist = hdu[0].header['HISTORY']
        except:
            raise cfis.CfisError('Error while reading tile FITS file {}'.format(img.name))

        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)p\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise cfis.CfisError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)
            exp_list.append(exp_name)

    exp_list_uniq = list(set(exp_list))

    return exp_list_uniq


def get_coord_at_image(number, band, image_type, images, no_cuts=False, verbose=False):
    """Return coordinate of image with given number.

    Parameters
    ----------
    number: string
        image number
    band: string
        optical band
    image_type: string
        image type ('tile', 'exposure', 'cat', weight')
    image: list of cfis.image
        list of images, used for type='exposure'
    no_cuts: bool, optional, default=False
        no cuts (of short exposure, validation flag) if True
    verbose: bool, optional
        verbose output if True, default=False

    Returns
    -------
    im_found: class cfis.image
        found image
    """

    img_found = None

    if image_type == 'tile':
        nix, niy = cfis.my_string_split(number, num=2, stop=True)
        tile_name = cfis.get_tile_name(nix, niy, band)

        if verbose == True:
            print('Looking for coordinates for tile with numbers ({},{})'.format(nix, niy))

        ra, dec   = cfis.get_tile_coord_from_nixy(nix, niy)
        img_found = cfis.image(tile_name, ra, dec)

    elif image_type == 'exposure':
        ra  = []
        dec = []
        for img in images:
            m = re.findall(number, img.name)
            if len(m) != 0:
                if img.cut(no_cuts=no_cuts) == False:
                    img_found = img

    else:
        raise cfis.CfisError('Image type \'{}\' not implemented yet'.format(image_type))

    return img_found

    

def test_tile_number():

    #ra  = Angle(180, unit=units.deg)
    #dec = Angle(35, unit=units.deg)
    ra  = Angle('10:44:00.0 hours')
    dec = Angle('30:00:00 degrees')
    print(ra, dec)
    print(ra.deg, dec.deg)

    nix, niy = cfis.get_tile_number_from_coord(ra, dec)
    print(nix, niy)

    ra, dec = cfis.get_tile_coord_from_nixy(nix, niy)
    print(ra, dec)


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
        input  = '.',
        input_format = 'full',
        mode  = 'c',
        band  = 'r',
        image_type  = 'tile',
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
    parser.add_option('-i', '--input', dest='input', type='string', default=p_def.input,
         help='input image list, can be ascii file or directory path')
    parser.add_option('-c', '--column', dest='col', type='string', default=None,
         help='column name if input is file, default=file has only one column)')
    parser.add_option('', '--input_format', dest='input_format', type='string',
         default=p_def.input_format,
         help='input format, one of \'full\', \'ID_only\', default=\'{}\''.format(p_def.input_format))
    parser.add_option('-o', '--outbase', dest='outbase', type='string', default=None,
         help='output file name base (\'.txt\' is added), default=stdout')
    parser.add_option('', '--plot', dest='plot', action='store_true',
         help='create plots')
    parser.add_option('', '--out_base_name', dest='out_base_name', action='store_true',
         help='output base names, not entire path if input is directory')
    parser.add_option('', '--out_name_only', dest='out_name_only', action='store_true',
         help='output only file names, not coordinates and metainfo')
    parser.add_option('', '--out_ID_only', dest='out_ID_only', action='store_true',
         help='output only file IDs, not full file names')
    parser.add_option('', '--interactive', dest='interactive', action='store_true',
         help='interactive mode (showing plots, recommended for call from jupyer notebook)')
    parser.add_option('-s', '--short', dest='short', action='store_true', help='short output')

    # Job control
    parser.add_option('', '--no_cuts', dest='no_cuts', action='store_true',
        help='output all exposures, no cuts (default: cut short and invalid exposures)')

    # Field and image options
    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='image_type', type='string', default=p_def.image_type,
        help='image type, one of \'tile\' (default)|\'weight\'|\'weight.fz\'|\'exposure\'|\'exposure_weight\''
             '|\'exposure_weight.fz\'|\'exposure_flag\'|\'exposure_flag.fz\'|\'cat\'')

    parser.add_option('', '--coord', dest='coord', type='string', default=None,
        help='(white-space or \'_\' separated) string of input coordinates, as astropy.coordinates.Angle')
    parser.add_option('', '--number', dest='number', type='string', default=None,
        help='input image number')
    parser.add_option('', '--area', dest='area', type='string', default=None,
        help='area corner coordinates ra0_dec0_ra1_dec1')
    parser.add_option('', '--tile', dest='tile', action='store_true',
        help='return exposures used in input tile(s)')

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

    Returns
    -------
    erg: bool
        Result of option check. False if invalid option value.
    """

    if int(options.number != None) + int(options.coord != None) \
        + int(options.area != None) + int(options.tile != None) > 1:
        raise cfis.CfisError('Only one option out of \'--number\', \'--coord\', \'--area\', \'--tile\' can be given')

    if options.image_type != 'exposure' and options.no_cuts == True:
        raise cfis.CfisError('option \'--no_cuts\' only possible for image_type=exposure')

    if options.input in ['{}.txt'.format(options.outbase), '{}.pdf',format(options.outbase)]:
        raise cfis.CfisError('Output base same as input, latter will be overwritten!')

    if options.input_format not in ['full', 'ID_only']:
        raise cfis.CfisError('input_format needs to be one of \'fulll\', \'ID_only\'')

    see_help = 'See option \'-h\' for help.'

    return True


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

    # Do extra stuff if necessary
    if param.outbase is not None:
        param.fout = open('{}.txt'.format(param.outbase), 'w')
    else:
        param.fout = sys.stdout

    return param


def run_mode(images, param):
    """Performs action according to run mode.

    Parameters
    ----------
    images: list of class cfis.images
        list of images
    param: class param
        parameter values

    Returns
    -------
    ex: int
        exit code, 0 if successful
    """

    # Default value
    ex = 1

    if param.coord:

        # Image number search: Return name of image(s) covering input coordinate
        images_found = cfis.find_image_at_coord(images, param.coord, param.band, param.image_type,
                                                no_cuts=param.no_cuts, input_format=param.input_format,
                                                verbose=param.verbose)
        if len(images_found) > 0:
            if not param.short:
                images_found[0].print_header(file=param.fout)
            for img in images_found:
                img.print(file=param.fout, base_name=param.out_base_name, name_only=param.out_name_only, ID_only=param.out_ID_only)
            ex = 0

    elif param.number:

        # Coordinate search: Return coordinate covered by image with input number/file name
        img_found = get_coord_at_image(param.number, param.band, param.image_type, images, no_cuts=param.no_cuts, verbose=param.verbose)
        if img_found != None:
            if not param.short:
                img_found.print_header(file=param.fout)
            img_found.print(file=param.fout, base_name=param.out_base_name, name_only=param.out_name_only, ID_only=param.out_ID_only)
            ex = 0
        else:
            if param.verbose:
                print('No image found, try with --no_cuts', file=sys.stderr)

    elif param.area:

        # Area search: Return images within input area
        angles = cfis.get_Angle_arr(param.area, num=4, verbose=param.verbose)
        images_found = cfis.find_images_in_area(images, angles, param.band, param.image_type, no_cuts=param.no_cuts, verbose=param.verbose)
        if len(images_found) > 0:
            if not param.short:
                images_found[0].print_header(file=param.fout)
            for img in images_found:
                img.print(file=param.fout, base_name=param.out_base_name, name_only=param.out_name_only, ID_only=param.out_ID_only)
            if param.plot == True:
                if param.verbose == True:
                    print('Creating plots')
                ra_c, dec_c, radius = cfis.plot_area(images_found, angles, param.image_type, param.outbase, param.interactive, show_numbers=True)

                if param.verbose:
                    print('RA_c[deg] DEC_c[deg] radius[argmin] = {:.2f} {:.2f} {:.2f}'.format(ra_c.deg, dec_c.deg, radius*60))
            ex = 0

    elif param.tile:

        # Search exposures used in input tile(s)
        images_found = get_images_used_in_tiles(images, param.band, param.image_type)
        if len(images_found) > 0:
            for img in images_found:
                print(img, file=param.fout)
            ex = 0


    else:
        raise cfis.CfisError('One of \'--coord\', \'--number\', '
                             '\'--area\', \'--tile\' needs to be specified')

    return ex


def main(argv=None):
    """Main program.
    """

    
    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)
    # Without option parsing, this would be: args = argv[1:]

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)


    # Save calling command
    cfis.log_command(argv)
    if param.verbose:
        cfis.log_command(argv, name='sys.stdout')


    if param.verbose and not param.short:
        print('Start of program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    if param.tile:
        # For search of exposures in tiles, image list from directory needs to be tiles
        image_type = 'tile'
    else:
        image_type = param.image_type
    images = cfis.get_image_list(param.input, param.band, image_type, col=param.col,
                                 input_format=param.input_format, verbose=param.verbose)


    # Check wether images have been found, if necessary
    if param.number:
        if images is None:
            raise cfis.CfisError('Input list file \'{}\' not found, neither '
                                 'existing file nor directory'.format(param.input))
        if len(images) == 0:
            raise cfis.CfisError('No corresponding image files found in input \'{}\''.format(param.input))


    # Run
    ex = run_mode(images, param)


    if param.outbase is not None:
        param.fout.close()
    

    ### End main program

    if param.verbose and not param.short:
        print('End of program {}'.format(os.path.basename(argv[0])))

    return ex



if __name__ == "__main__":
    sys.exit(main(sys.argv))

