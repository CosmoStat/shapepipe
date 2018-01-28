#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module cfis.py

CFIS module

:Authors: Martin Kilbinger

:Date: 19/01/2018
"""

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
size['exposure'] = 1


class image():

    def __init__(self, name, ra, dec):
        self.name = name
        self.ra   = ra
        self.dec  = dec

    def print(self, **kwds):
        print(self.__dict__)



def get_file_pattern(pattern, band, image_type):
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
        pattern = '{}\.weight\.fits\.fz'.format(pattern_base)
    else:
        stuff.error('Invalid type \'{}\''.format(image_type))

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

    return Angle(ra), Angle(dec)



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

    file_list.sort()
    return file_list


def create_image_list(fname, ra, dec):
    """Return list of image information.

    Parameters
    ----------
    fname: list of strings
        file names
    ra: list of strings
        right ascension
    dec: list of strings
        declination

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
        im = image(fname[i], r, d)
        images.append(im)

    return images



