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

from astropy.io import fits
from astropy.table import Table, Column
from astropy import units
from astropy.coordinates import Angle, SkyCoord

from optparse import OptionParser, IndentedHelpFormatter, OptionGroup
import textwrap


from cfis import CfisError, param, unitdef, size
import cfis



class IndentedHelpFormatterWithNL(IndentedHelpFormatter):
  """Allows newline to have effect in option help.
     From https://groups.google.com/forum/#!msg/comp.lang.python/bfbmtUGhW8I/sZkGryaO8gkJ
     Usage: parser = OptionParser(usage=usage, formatter=IndentedHelpFormatterWithNL())
  """
  def format_description(self, description):
    if not description: return ""
    desc_width = self.width - self.current_indent
    indent = " "*self.current_indent
# the above is still the same
    bits = description.split('\n')
    formatted_bits = [
      textwrap.fill(bit,
        desc_width,
        initial_indent=indent,
        subsequent_indent=indent)
      for bit in bits]
    result = "\n".join(formatted_bits) + "\n"
    return result

  def format_option(self, option):
    # The help for each option consists of two parts:
    #   * the opt strings and metavars
    #   eg. ("-x", or "-fFILENAME, --file=FILENAME")
    #   * the user-supplied help string
    #   eg. ("turn on expert mode", "read data from FILENAME")
    #
    # If possible, we write both of these on the same line:
    #   -x    turn on expert mode
    #
    # But if the opt string list is too long, we put the help
    # string on a second line, indented to the same column it would
    # start in if it fit on the first line.
    #   -fFILENAME, --file=FILENAME
    #       read data from FILENAME
    result = []
    opts = self.option_strings[option]
    opt_width = self.help_position - self.current_indent - 2
    if len(opts) > opt_width:
      opts = "%*s%s\n" % (self.current_indent, "", opts)
      indent_first = self.help_position
    else: # start help on same line as opts
      opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
      indent_first = 0
    result.append(opts)
    if option.help:
      help_text = self.expand_default(option)
# Everything is the same up through here
      help_lines = []
      for para in help_text.split("\n"):
        help_lines.extend(textwrap.wrap(para, self.help_width))
# Everything is the same after here
      result.append("%*s%s\n" % (
        indent_first, "", help_lines[0]))
      result.extend(["%*s%s\n" % (self.help_position, "", line)
        for line in help_lines[1:]])
    elif opts[-1] != "\n":
      result.append("\n")
    return "".join(result)


def find_image_at_coord(images, coord, band, image_type, no_cuts=False, verbose=False):
    """Return image covering given coordinate.

    Parameters
    ----------
    images: list of class cfis.image
        list of images
    coord: string
        coordinate ra and dec with units
    band: string
        optical band
    image_type: string
        image type ('tile', 'weight', 'weight.fz', 'exposure', 'exposure_weight', \
        'exposure_weight.fz', 'exposure_flag', 'exposure_flag.fz', 'cat')
    no_cuts: bool, optional, default=False
        no cuts (of short exposure, validation flag) if True
    verbose: bool, optional
        verbose output if True, default=False

    Returns
    -------
    img_found: list of cfis.image
        Found image(s), None if none found.
    """

    ra, dec = cfis.get_Angle(coord)

    if verbose == True:
        print('Looking for image at coordinates {}, {}'.format(ra, dec))

    if image_type in ('tile', 'weight', 'weight.fz'):
        nix, niy  = cfis.get_tile_number_from_coord(ra, dec, return_type=int)
        tile_name = cfis.get_tile_name(nix, niy, band, image_type)

        img_found = []
        for img in images:
            if os.path.basename(img.name) == tile_name:
                # Update coordinate in image for tiles with central coordinates
                ra_c, dec_c = cfis.get_tile_coord_from_nixy(nix, niy)
                if img.ra is not None or img.dec is not None:
                    raise CfisError('Coordinates in image are already set to {}, {}, cannot update to {}, {}'.\
                                format(img.ra, img.dec, ra_c, dec_c))
                img.ra = ra_c
                img.dec = dec_c
                img_found.append(img)

        if len(img_found) != 0:
                pass
        else:
            if verbose == True:
                print('Tile with numbers ({}, {}) not found'.format(nix, niy))

        if len(img_found) > 1:
            raise CfisError('More than one tile ({}) found'.format(len(img_found)))

    elif image_type == 'exposure':
        sc_input = SkyCoord(ra, dec)

        img_found = []
        for img in images:
            # Check distance along ra and dec from image center
            sc_img_same_ra  = SkyCoord(ra, img.dec)
            sc_img_same_dec = SkyCoord(img.ra, dec)
            distance_ra  = sc_input.separation(sc_img_same_dec)
            distance_dec = sc_input.separation(sc_img_same_ra)
            if distance_ra.degree < size[image_type]/2 and distance_dec.degree < size[image_type]/2:
                if img.cut(no_cuts=no_cuts) == False:
                    img_found.append(img)

        if len(img_found) != 0:
                pass
        else:
            if verbose == True:
                print('No exposure image found')

    else:
        raise CfisError('Only implemented for image_type=tile')

    return img_found



def find_images_in_area(images, angles, band, image_type, no_cuts=False, verbose=False):
    """Return image list within coordinate area (rectangle)

    Parameters
    ----------
    images: list of class cfis.image
        list of images
    angles: string
        coordinates ra0_dec0_ra1_dec1 with units
    band: string
        optical band
    image_type: string
        image type ('tile', 'exposure', 'cat', 'weight', 'weight.fz')
    no_cuts: bool, optional, default=False
        no cuts (of short exposure, validation flag) if True
    verbose: bool, optional, default=False
        verbose output if True
`
    Returns
    -------
    found: list of cfis.image
        found images
    """

    if verbose == True:
        print('Looking for all images within rectangle, lower left=({},{}), upper right=({},{}) deg'.format(
              angles[0].ra.deg, angles[0].dec.deg, angles[1].ra.deg, angles[1].dec.deg))

    found = []

    if image_type in ('tile', 'weight', 'weight.fz'):
        for img in images:
            nix, niy = cfis.get_tile_number(img.name)
            ra, dec  = cfis.get_tile_coord_from_nixy(nix, niy)
            # MKDEBUG TODO: bounds around 0 deg
            if ra.is_within_bounds(angles[0].ra, angles[1].ra) \
                and dec.is_within_bounds(angles[0].dec, angles[1].dec):
                # Update coordinate in image class. This could be done for all images,
                # not just the found ones.
                if img.ra is not None or img.dec is not None:
                    raise CfisError('Coordinates in image are already set to {}, {}, cannot update to {}, {}'.\
                                format(img.ra, img.dec, ra, dec))
                img.ra  = ra
                img.dec = dec
                found.append(img)

    elif image_type == 'exposure':
        for img in images:
            if img.ra.is_within_bounds(angles[0].ra, angles[1].ra) \
                and img.dec.is_within_bounds(angles[0].dec, angles[1].dec):

                if img.cut(no_cuts=no_cuts) == False:
                    found.append(img)

    else:
        raise CfisError('Image type \'{}\' not implemented yet'.format(image_type))

    if verbose == True:
        print('{} images found in area'.format(len(found)))

    return found


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
            raise CfisError('Error while reading tile FITS file {}'.format(img.name))

        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)p\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise CfisError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

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
        raise CfisError('Image type \'{}\' not implemented yet'.format(image_type))

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



def plot_area(images, angles, image_type, outbase, interactive):
    """Plot images within area.

    Parameters
    ----------
    images: array of cfis.image
        images
    angles: array(SkyCoord, 2)
        Corner coordinates of area rectangle
    image_type: string
        image type ('tile', 'exposure', 'cat', weight')
    outbase: string
        output file name base
    interactive: bool
        show plot if True
    """

    if outbase is None:
        outname = 'plot.pdf'
    else:
        outname = '{}.pdf'.format(outbase)

    lw = 0.25
    color = {'tile': 'b', 'exposure': 'g', 'weight': 'r'}

    ax = plot_init()

    # Field center
    n_ima = len(images)
    if n_ima > 0:
        ra_c  = sum([img.ra for img in images])/float(n_ima)
        dec_c = sum([img.dec for img in images])/float(n_ima)
        plt.plot(ra_c, dec_c, 'or', mfc='none', ms=3)
    else:
        ra_c = 0
        dec_c = 0

    # Circle around field
    dx = abs(angles[0].ra - angles[1].ra)
    dy = abs(angles[0].dec - angles[1].dec)
    dx = getattr(dx, unitdef)
    dy = getattr(dy, unitdef)
    radius = max(dx, dy)/2 + (cfis.size['exposure'] + cfis.size['tile']) * np.sqrt(2)
    circle = plt.Circle((ra_c.deg, dec_c.deg), radius, color='r', fill=False)
    ax.add_artist(circle)

    for img in images:
        # Image center
        x  = img.ra.degree
        y  = img.dec.degree
        #plt.plot(x, y, 'b.', markersize=1)

        # Image boundary
        dx = size[image_type] / 2
        dy = size[image_type] / 2
        cx, cy = square_from_centre(x, y, dx, dy)
        plt.plot(cx, cy, '{}-'.format(color[image_type]), linewidth=lw)

    # Area border
    cx, cy = square_from_corners(angles[0], angles[1])
    plt.plot(cx, cy, 'r-.', linewidth=lw)

    plt.xlabel('R.A. [degree]')
    plt.ylabel('Declination [degree]')
    if outbase is not None:
        plt.title(outbase)

    # Limits
    border = 2
    xm = (angles[1].ra.degree + angles[0].ra.degree) / 2
    ym = (angles[1].dec.degree + angles[0].dec.degree) / 2
    dx = angles[1].ra.degree - angles[0].ra.degree
    dy = angles[1].dec.degree - angles[0].dec.degree
    lim = max(dx, dy)
    plt.xlim(xm - lim/2 - border, xm + lim/2 + border)
    plt.ylim(ym - lim/2 - border, ym + lim/2 + border)

    # Somehow this does not work (any more?)
    #limits = plt.axis('equal')
    #print(limits)

    print('Saving plot to {}'.format(outname))
    plt.savefig(outname)

    if interactive == True:
        plt.show()

    return ra_c, dec_c, radius


def plot_init():

    fs = 12
    fig, ax = plt.subplots()
    #fig = plt.figure()

    ax = plt.gca()
    ax.yaxis.label.set_size(fs)
    ax.xaxis.label.set_size(fs)

    plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.rcParams.update({'figure.autolayout': True})

    #return fig, ax
    return ax



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

    p_def = param(
        input  = '.',
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
    parser = OptionParser(usage=usage, formatter=IndentedHelpFormatterWithNL())

    # I/O
    parser.add_option('-i', '--input', dest='input', type='string', default=p_def.input,
         help='input image list, can be ascii file or directory path')
    parser.add_option('-c', '--column', dest='col', type='string', default=None,
         help='column name if input is file, default=file has only one column)')
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
        raise CfisError('Only one option out of \'--number\', \'--coord\', \'--area\', \'--tile\' can be given')

    if options.image_type != 'exposure' and options.no_cuts == True:
        raise CfisError('option \'--no_cuts\' only possible for image_type=exposure')

    if options.input in ['{}.txt'.format(options.outbase), '{}.pdf',format(options.outbase)]:
        raise CfisError('Output base same as input, latter will be overwritten!')

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
        images_found = find_image_at_coord(images, param.coord, param.band, param.image_type, no_cuts=param.no_cuts, verbose=param.verbose)
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
        images_found = find_images_in_area(images, angles, param.band, param.image_type, no_cuts=param.no_cuts, verbose=param.verbose)
        if len(images_found) > 0:
            if not param.short:
                images_found[0].print_header(file=param.fout)
            for img in images_found:
                img.print(file=param.fout, base_name=param.out_base_name, name_only=param.out_name_only, ID_only=param.out_ID_only)
            if param.plot == True:
                if param.verbose == True:
                    print('Creating plots')
                ra_c, dec_c, radius = plot_area(images_found, angles, param.image_type, param.outbase, param.interactive)

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
        raise CfisError('One of \'--coord\', \'--number\', \'--area\', \'--tile\' needs to be specified')

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
    images = cfis.get_image_list(param.input, param.band, image_type, col=param.col, verbose=param.verbose)


    # Check wether images have been found, if necessary
    if param.number:
        if images is None:
            raise CfisError('Input list file \'{}\' not found, neither existing file nor directory'.format(param.input))
        if len(images) == 0:
            raise CfisError('No corresponding image files found in input \'{}\''.format(param.input))


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

