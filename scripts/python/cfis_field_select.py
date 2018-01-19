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

from astropy.io import ascii
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle

from optparse import OptionParser, IndentedHelpFormatter
from optparse import OptionGroup

import cfis
from cfis import unitdef
import stuff


#unitdef = 'degree'



def get_image_list(inp, band, image_type, verbose=False):
    """Return list of images.

    Parameters
    ----------
    input: string
        file name or direcory path
    verbose: bool, optional
        verbose output if True, default=False

    Return
    ------
    img_list: list of strings
        image list
    """

    # Get file list
    if os.path.isfile(inp):
        inp_type  = 'file'
        file_list = cfis.read_list(inp)

    elif os.path.isdir(inp):
        inp_type  = 'dir'
        file_list = glob.glob(inp)

    else:
        file_list = None


    # Filter file list to match CFIS image pattern
    img_list = []
    pattern = cfis.get_file_pattern('', band, image_type)

    for f in file_list:

        m = re.findall(pattern, f)
        if len(m) != 0:
            img_list.append(m[0])

    if verbose == True and len(img_list) > 0:
        print('{} image files found in input {} \'{}\''.format(len(img_list), inp_type, inp))

    return img_list



def find_image_at_coord(images, coord, band, image_type, verbose=False):

    ra, dec   = cfis.get_Angle(coord)

    if verbose == True:
        print('Looking for image at coordinates {}, {}'.format(ra, dec))

    if image_type == 'tiles':
        nix, niy  = cfis.get_tile_number_from_coord(ra, dec, return_type=int)
        tile_name = cfis.get_tile_name(nix, niy, band)
        im_found = None

        for im in images:
            if im == tile_name:
                im_found = im
                break

        if im_found != None:
                pass
        else:
            if verbose == True:
                print('Image with numbers ({}, {}) not found'.format(nix, niy))

    else:
        stuff.error('Only implemented for image_type=tiles')

    return im_found



def find_images_in_area(images, angles, band, image_type, verbose=False):

    if verbose == True:
        print('Looking for all images within coordinates ', angles)

    if image_type == 'tiles':
        found = []
        for img in images:
            nix, niy = cfis.get_tile_number(img)
            ra, dec  = cfis.get_tile_coord_from_nixy(nix, niy)
            if ra.is_within_bounds(angles[0].ra, angles[1].ra) \
                and dec.is_within_bounds(angles[0].dec, angles[1].dec):
                img = cfis.image(img, ra, dec)
                found.append(img)

    else:
        stuff.error('Image type {} not implemented yet'.format(image_type))

    if verbose == True:
        print('{} images found in area'.format(len(found)))

    return found


def get_coord_at_image(number, image_type, verbose=False):


    if image_type == 'tiles':
	# TODO: Read entire image name, not just the two tile numbers
        nix, niy = stuff.my_string_split(number, num=2, stop=True)

        if verbose == True:
            print('Looking for coordinates for tile with numbers ({},{})'.format(nix, niy))

        ra, dec  = cfis.get_tile_coord_from_nixy(nix, niy)
    else:
        stuff.error('Image type {} not implemented yet'.format(image_type))

    return ra, dec

    

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


def square_from_centre(x, y, dx, dy):
    """Return coordinate vectors of corners cx, cy that define a closed square for plotting.
    """

    cx = [x-dx, x+dx, x+dx, x-dx, x-dx]
    cy = [y-dy, y-dy, y+dy, y+dy, y-dy]

    return cx, cy



def square_from_corners(ang0, ang1):
    """Return coordinate vectors of corners cx, cy that define a closed square for plotting.
    """

    cx = [ang0.ra, ang1.ra, ang1.ra, ang0.ra, ang0.ra]
    cy = [ang0.dec, ang0.dec, ang1.dec, ang1.dec, ang0.dec]

    cxd = [getattr(i, unitdef) for i in cx]
    cyd = [getattr(i, unitdef) for i in cy]

    return cxd, cyd



def plot_area(images, angles, outbase):
    """Plot images within area.

    Parameters
    ----------
    images: array of cfis.image
        images
    angles: array(SkyCoord, 2)
        Corner coordinates of area rectangle
    outbase: string
        output file name base
    """

    if outbase is None:
        outname = 'plot.pdf'
    else:
        outname = '{}.pdf'.format(outbase)

    fig, ax = plot_init()

    # Field center
    ra_c  = sum([img.ra for img in images])/len(images)
    dec_c = sum([img.dec for img in images])/len(images)
    plt.plot(ra_c, dec_c, 'or', mfc='none', ms=3)

    for img in images:
        # Image center
        x  = img.ra.degree
        y  = img.dec.degree
        plt.plot(x, y, 'b.', markersize=1)

        # Image boundary
        dx = 0.25
        dy = 0.25
        cx, cy = square_from_centre(x, y, dx, dy)
        plt.plot(cx, cy, 'g-', linewidth=0.5) 

    # Area border
    cx, cy = square_from_corners(angles[0], angles[1])
    plt.plot(cx, cy, 'r-.', linewidth=0.5)

    plt.xlabel('R.A. [degree]')
    plt.ylabel('Declination [degree]')
    if outbase is not None:
        plt.title(outbase)

    plt.axis('equal')

    print('Saving plot to {}'.format(outname))
    plt.savefig(outname)
    plt.show()



def plot_init():

    fs = 12
    fig = plt.figure()

    ax = plt.gca()
    ax.yaxis.label.set_size(fs)
    ax.xaxis.label.set_size(fs)

    plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.rcParams.update({'figure.autolayout': True})

    return fig, ax



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
        input  = '.',
        mode  = 'c',
        band  = 'r',
        image_type  = 'tiles',
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

    # I/O
    parser.add_option('-i', '--input', dest='input', type='string', default=p_def.input,
         help='input image list, can be ascii file or directory path')
    parser.add_option('-o', '--outbase', dest='outbase', type='string', default=None,
         help='output file name base (\'.txt\' is added), default=stdout')
    parser.add_option('', '--plot', dest='plot', action='store_true',
         help='create plots')

    # Job control
    parser.add_option('-m', '--mode', dest='mode', type='string', default=None,
         help='run mode, one of [n|c|a]:\n'
	      ' n: search image number given a coordinate (--coord)\n'
	      ' c: search for image coord. given a number (--number)\n'
	      ' a: area search, search images within area (--area)\n')

    # Field and image options
    parser.add_option('-b', '--band', dest='band', type='string', default=p_def.band,
        help='band, one of \'r\' (default)|\'u\'')
    parser.add_option('-t', '--type', dest='image_type', type='string', default=p_def.image_type,
        help='image type, one of \'tiles\' (default)| \'cat\'|\'weight\'|\'raw\'')

    parser.add_option('', '--coord', dest='coord', type='string', default=None,
        help='(white-space or \'_\' separated) string of input coordinates, as astropy.coordinates.Angle')
    parser.add_option('', '--number', dest='number', type='string', default=None,
        help='input image number')
    parser.add_option('', '--area', dest='area', type='string', default=None,
        help='area corner coordinates ra0_dec0_ra1_dec1')

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

    if options.mode is None:
        stuff.error('Run mode not given (option \'-m\')')

    if options.mode == 'n' and options.coord == None:
        stuff.error('No input coordinates given (option \'--coord\')')
    if options.mode == 'c' and options.number == None:
        stuff.error('No input image number given (option \'--number\')')
    if options.mode == 'a' and options.area == None:
        stuff.error('No input area given (option \'--area\')')

    see_help = 'See option \'-h\' for help.'

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
    images: list of strings
        list of images
    param: class stuff.param
        parameter values

    Returns
    -------
    ex: int
        exit code, 0 if successful
    """

    # Default value
    ex = 1

    if param.mode == 'n':

        # Image number search: Return number of tile that covers input coordinate
        image_name = find_image_at_coord(images, param.coord, param.band, param.image_type, verbose=param.verbose)
        if image_name != None:
            print('# name')
            print(image_name, file=param.fout) 
            ex =  0

    elif param.mode == 'c':

        # Coordinate search: Return coordinate covered by tile with input number
        ra, dec = get_coord_at_image(param.number, param.image_type, verbose=param.verbose)
        print('# ra[{0}] dec[{0}]'.format(unitdef), file=param.fout)
        print(getattr(ra, unitdef), getattr(dec, unitdef), file=param.fout)
        ex = 0

    elif param.mode == 'a':

        # Area search: Return tiles within input area
        ex = 0
        angles = cfis.get_Angle_arr(param.area, num=4, verbose=param.verbose)
        images = find_images_in_area(images, angles, param.band, param.image_type, verbose=param.verbose)
        print('# Name ra[{0}] dec[{0}]'.format(unitdef), file=param.fout)
        for img in images:
            print('{} {} {}'.format(img.name, getattr(img.ra, unitdef), getattr(img.dec, unitdef)), file=param.fout)
        if param.plot == True:
            if param.verbose == True:
                print('Creating plots')
            plot_area(images, angles, param.outbase)
        ex = 0

    else:

        stuff.error('Unknown run mode \'{}\''.format(param.mode))


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
    stuff.log_command(argv)
    if param.verbose:
        stuff.log_command(argv, name='sys.stderr')


    if param.verbose is True:
        print('Start program {}'.format(os.path.basename(argv[0])))


    ### Start main program ###

    images = get_image_list(param.input, param.band, param.image_type, verbose=param.verbose)


    # Check wether images have been found, if necessary for run mode
    if param.mode != 'c':
        if images is None:
            stuff.error('Input {} not found, neither existing file nor directory'.format(param.input))
        if len(images) == 0:
            stuff.error('No image files found in input \'{}\''.format(param.input))


    # Run
    ex = run_mode(images, param)


    if param.outbase is not None:
        param.fout.close()
    

    ### End main program

    if param.verbose is True:
        print('Finish program {}'.format(os.path.basename(argv[0])))

    return ex



if __name__ == "__main__":
    sys.exit(main(sys.argv))

