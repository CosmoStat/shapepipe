#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module cfis.py

CFIS module

:Authors: Martin Kilbinger

:Date: 19/01/2018
"""

# Compability with python2.x for x>6
from __future__ import print_function


import re
import sys

import numpy as np

from astropy import units
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord

import stuff


unitdef = 'degree'

# Maybe define class for these constants?
size = {}
size['tile']     = 0.5
size['weight']   = 0.5
size['exposure'] = 1.0

# Cut criteria for exposures
exp_time_min     = 95
flag_valid       = 'V'


class image():

    def __init__(self, name, ra, dec, exp_time=-1, valid='Unknown'):
        """Create image information.

        Parameters
        ----------
        name: string
            file name
        ra: Angle
            right ascension
        dec: Angle
            declination
        exp_time: integer, optiona, default=-1
            exposure time
        valid: string, optional, default='Unknown'
            validation flag

        Returns
        -------
        self: class image
            image information
        """
            
        self.name     = name
        self.ra       = ra
        self.dec      = dec
        if exp_time == None:
            self.exp_time = -1
        else:
            self.exp_time = exp_time
        if valid == None:
            self.valid = 'Unknown'
        else:
            self.valid = valid


    def cut(self, no_cuts=False):
        """Return True (False) if image does (not) need to be cut from selection.

        Parameters
        ----------
        no_cuts: bool, optiona, default=False
            do not cut if True

        Returns
        -------
        cut: bool
            True (False) if image is (not) cut
        """

        # Do not cut if no_cuts flag is set
        if no_cuts == True:
            return False

        # Cut if exposure time smaller than minimum (and not flagged as unknown or n/a)
        if self.exp_time < exp_time_min and self.exp_time != -1:
            return True

        # Cut if validation flag is not valid (and not unknown)
        if self.valid != flag_valid and self.valid != 'Unknown':
            return True

        return False


    def print(self, file=sys.stdout):
    #def print(self, **kwds):
        """Print image information as ascii Table column

        Parameters
        ----------
        file: file handle, optional, default=sys.stdout
            output file handle

        Returns
        -------
        None
        """

        print('{} {:10.2f} {:10.2f} {:5d} {:8s}'.format(self.name, getattr(self.ra, unitdef), getattr(self.dec, unitdef), \
              self.exp_time, self.valid), file=file)


    def print_header(self, file=sys.stdout):
        """Print header for ascii Table output

        Parameters
        ----------
        file: file handle, optional, default=sys.stdout
            output file handle

        Returns
        -------
        None
        """

        print('# Name ra[{0}] dec[{0}] exp_time[s] validation'.format(unitdef), file=file)



def get_file_pattern(pattern, band, image_type, want_re=True):
    """Return file pattern of CFIS image file.
    """

    if pattern == '':
        if image_type == 'exposure':
            pattern_base = '\d{7}p'
        else:
            pattern_base  = 'CFIS.*\.{}'.format(band)
    else:
        pattern_base = pattern

    if image_type == 'exposure':
        pattern  = '{}\.fits.fz'.format(pattern_base)
    elif image_type == 'tile':
        pattern = '{}\.fits'.format(pattern_base)
    elif image_type == 'cat':
        pattern = '{}\.cat'.format(pattern_base)
    elif image_type == 'weight':
        pattern = '{}\.weight\.fits'.format(pattern_base)
    elif image_type == 'weight.fz':
        pattern = '{}\.weight\.fits.fz'.format(pattern_base)
    else:
        stuff.error('Invalid type \'{}\''.format(image_type))

    if want_re == False:
        pattern = pattern.replace('\\', '')

    return pattern


def get_tile_number_from_coord(ra, dec, return_type=str):
    """Return CFIS stacked image tile number covering input coordinates.
        This is the inverse to get_tile_coord_from_nixy.

    Parameters
    ----------
    ra: Angle
        right ascension
    dec: Angle
        declination
    return type: <type 'type'>
	return type, int or str

    Returns
    -------
    nix: string
        tile number for x
    niy: string
        tile number for y
    """

    y = (dec.degree + 90) * 2.0
    yi = int(np.rint(y))

    x = ra.degree * np.cos(dec.radian) * 2.0
    #x = ra.degree * 2 * np.cos(y/2 / 180 * np.pi - np.pi/2)
    xi = int(np.rint(x))
    if xi == 720:
        xi = 0

    if return_type == str:
        nix = '{:03d}'.format(xi)
        niy = '{:03d}'.format(yi)
    elif return_type == int:
        nix = xi
        niy = yi
    else:
        stuff.error('Invalid return type {}'.format(return_type))

    return nix, niy


def get_tile_coord_from_nixy(nix, niy):
    """ Return coordinates corresponding to tile with number (nix,niy).
        This is the inverse to get_tile_number_from_coord.
    Parameters
    ----------
    nix: string
        tile number for x
    niy: string
        tile number for y

    Returns
    -------
    ra: Angle
        right ascension
    dec: Angle
        declination
    """

    xi = int(nix)
    yi = int(niy)

    d = (yi/2.0 - 90)
    dec = Angle('{} degrees'.format(d))
    r = xi / 2.0 / np.cos(dec.radian)
    ra = Angle('{} degrees'.format(r))

    return ra, dec



def get_tile_name(nix, niy, band):
    """Return tile name for given tile numbers.

   Parameters
   ----------
    nix: string
        tile number for x
    niy: string
        tile number for y
    band: string
        band, one in 'r' or 'u'

    Returns
    -------
    tile_name: string
        tile name
    """

    if type(nix) is int and type(niy) is int:
    	tile_name = 'CFIS.{:03d}.{:03d}.{}.fits'.format(nix, niy, band)

    elif type(nix) is str and type(niy) is str:
    	tile_name = 'CFIS.{}.{}.{}.fits'.format(nix, niy, band)

    else:
        stuff.error('Invalid type for input tile numbers {}, {}'.format(nix, niy))


    return tile_name



def get_tile_number(tile_name):
    """Return tile number of given image tile name

    Parameters
    ----------
    tile_name: string
        tile name

    Returns
    -------
    nix: string
        tile number for x
    niy: string
        tile number for y
    """

    m = re.search('CFIS\.(\d{3})\.(\d{3})', tile_name)
    if m == None or len(m.groups()) != 2:
        stuff.error('Image name \'{}\' does not match tile name syntax'.format(tile_name))

    nix = m.groups()[0]
    niy = m.groups()[1]

    return nix, niy
    

def check_ra(ra):
    """Range check of right ascension.

    Parameters
    ----------
    ra: Angle
        right ascension

    Returns
    -------
    res: bool
        result of check (True if pass, False if fail)
    """

    print(ra.deg)
    if ra.deg < 0 or ra.deg > 360:
        stuff.error('Invalid ra, valid range is 0 < ra < 360 deg')
        return 1

    return 0


def check_dec(dec):
    """Range check of declination.

    Parameters
    ----------
    dec: Angle
        declination

    Returns
    -------
    res: bool
        result of check (True if pass, False if fail)
    """

    if dec.deg < -90 or dec.deg > 90:
        stuff.error('Invalid dec, valid range is -90 < dec < 90 deg')
        return 1

    return 0



def get_Angle(str_coord):
    """Return Angles ra, dec from coordinate string

    Parameters
    ----------
    str_coord: string
        string of input coordinates

    Returns
    -------
    ra: Angle
        right ascension
    dec: Angle
        declination
    """

    ra, dec = stuff.my_string_split(str_coord, num=2, stop=True)

    a_ra  = Angle(ra)
    a_dec = Angle(dec)

    return r_ra, a_dec



def get_Angle_arr(str_coord, num=-1, verbose=False):
    """Return array of Angles from coordinate string

    Parameters
    ----------
    str_coord: string
        string of input coordinates
    num: int, optional, default=-1
        expected number of coordinates (even number)
    verbose: bool, optional, default=False
        verbose output

    Returns
    -------
    angles: array of SkyCoord
        array of sky coordinates (pairs ra, dec)
    """

    angles_mixed = stuff.my_string_split(str_coord, num=num, verbose=verbose, stop=True)
    n = len(angles_mixed)
    n = int(n / 2)

    angles = []
    for i in range(n):
        c = SkyCoord(angles_mixed[2*i], angles_mixed[2*i+1])
        angles.append(c)

    return angles



def read_list(fname):
    """Read list of from ascii file.

    Parameters
    ----------
    fname: string
        ascii file name

    Returns
    -------
    file_list: list of strings
        list of file name
    """

    f = open(fname, 'rU')
    file_list = [x.strip() for x in f.readlines()]
    f.close()

    # from sapastro1:toolbox/python/CFIS.py
    #except IOError as exc:
        #if exc.errno == errno.ENOENT:
            #if verbose == True:
                #print('Not using exclude file')
            #out_list = []
        #else:
            #raise


    file_list.sort()
    return file_list


def create_image_list(fname, ra, dec, exp_time=[], valid=[]):
    """Return list of image information.

    Parameters
    ----------
    fname: list of strings
        file names
    ra: list of strings
        right ascension
    dec: list of strings
        declination
    exp_time: list of integers, optional, default=[]
        exposure time
    valid: list of strings, optional, default=[]
        QSO exposure validation flag

    Returns
    -------
    images: list of cfis.image
        list of image information
    """

    nf = len(fname)
    nr = len(ra)
    nd = len(dec)
    if nf == 0:
        stuff.error('No entries in file name list')
    if (nf != nr or nf != nd) and nr != 0 and nd != 0:
        stuff.error('Lists fname, ra, dec have not same length ({}, {}, {})'.format(nf, nr, nd))

    images = []
    for i in range(nf):
        if nr > 0 and nd > 0:
            r = Angle('{} {}'.format(ra[i], unitdef))
            d = Angle('{} {}'.format(dec[i], unitdef))
        else:
            r = None
            d = None
        if len(exp_time) > 0:
            e = exp_time[i]
        else:
            e = -1
        if len(valid) > 0:
            v = valid[i]
        else:
            v = None
        im = image(fname[i], r, d, exp_time=e, valid=v)
        images.append(im)

    return images


def get_exposure_info(logfile_name, verbose=False):  
    """Return information on run (single exposure) from log file.

    Parameters
    ----------
    logfile_name: string
        file name
    verbose: bool, optional, default=False
        verbose output

    Returns:
    images: list of class image
        list of exposures
    """

    images = []
    f = open(logfile_name)
    for line in f:
        dat = re.split(' |', line)
        name = dat[0]
        ra   = Angle(' hours'.format(dat[8]))
        dec  = Angle(' degree'.format(dat[9]))
        valid = dat[21]
    
        img = image(name, ra, dec, valid=valid)
        image.append(img)

    return image



def exclude(f, exclude_list):
    """Return True if f is on exclude_list

    Parameters
    ----------
    f: string
        file name
    exclude_list: list of strings
        list of files

    Returns
    -------
    is_in_exclude: bool
        True (False) if f is in list
    """

    return f in exclude_list

