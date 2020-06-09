#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script create_star_cat.py

Create reference star catalogue for masking of
bright stars and diffraction spikes

:Authors: Axel Guinot, Martin Kilbinger

:Date: 2019, 2020

:Package: ShapePipe
"""


from shapepipe.pipeline.execute import execute

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
from astropy import units as u
import re
import os
import sys


def _get_wcs(header):
    """Get WCS

    Compute the astropy WCS from header manually.
    (The purpose of this is to avoid possible incompatibility on distortion
    convention)

    Parameters
    ----------
    header : astropy.header
        Image header

    Returns
    -------
    final_wcs : astropy.wcs.WCS
        WCS object
    """

    final_wcs = WCS(naxis=2)
    final_wcs.wcs.ctype = [header['CTYPE1'], header['CTYPE2']]
    try:
        final_wcs.wcs.cunit = [header['CUNIT1'], header['CUNIT2']]
    except:
        final_wcs.wcs.cunit = ['deg', 'deg']
    final_wcs.wcs.crpix = [header['CRPIX1'], header['CRPIX2']]
    final_wcs.wcs.crval = [header['CRVAL1'], header['CRVAL2']]
    final_wcs.wcs.cd = [[header['CD1_1'], header['CD1_2']],
                        [header['CD2_1'], header['CD2_2']]]

    return final_wcs


def _get_image_radius(center, wcs):
    """Get image radius

    Compute the diagonal distance of the image in arcmin.

    Parameters
    ----------
    center : numpy.ndarray
        Coordinates of the center of the image (in pixel)

    Returns
    -------
    float
        The diagonal distance of the image in arcmin.
    """

    if center is None:
        return SphereDist(self._fieldcenter['pix'], np.zeros(2))/60.
    else:
        if type(center) is np.ndarray:
            return SphereDist(center, np.zeros(2), wcs)/60.
        else:
            raise TypeError('center has to be a numpy.ndarray')


def SphereDist(position1, position2, wcs):
    """Compute spherical distance

    Compute spherical distance between 2 points.

    Parameters
    ----------
    position1 : numpy.ndarray
        [x,y] first point (in pixel)
    position2 : numpy.ndarray
        [x,y] second point (in pixel)

    Returns
    -------
    float
        The distance in degree.
    """

    if (type(position1) is not np.ndarray) & (type(position2) is not np.ndarray):
        raise ValueError('Positions need to be a numpy.ndarray')

    p1 = (np.pi/180.)*np.hstack(wcs.all_pix2world(position1[0], position1[1], 1))
    p2 = (np.pi/180.)*np.hstack(wcs.all_pix2world(position2[0], position2[1], 1))

    dTheta = p1 - p2
    dLong = dTheta[0]
    dLat = dTheta[1]

    dist = 2*np.arcsin(np.sqrt(np.sin(dLat/2.)**2. + np.cos(p1[1])*np.cos(p2[1])*np.sin(dLong/2.)**2.))

    return dist*(180./np.pi)*3600.


def find_stars(position, output_name, radius=None):
    """Find stars

    Return GSC (Guide Star Catalog) objects for a field with center (ra,dec) and radius r.

    Parameters
    ----------
    position : numpy.ndarray
        Position of the center of the field
    radius : float
        Radius in which the query is done (in arcmin)

    Returns
    -------
    dict
        Stars dicotionnary for GSC objects in the field.
    """

    ra = position[0]
    dec = position[1]

    # check ra dec types

    if dec > 0.:
        sign = '+'
    else:
        sign = ''

    cmd_line = '{0} {1} {2}{3} -r {4} -n 1000000'.format('findgsc2.2', ra, sign, dec, radius)

    # output=subprocess.check_output(cmd_line, shell=True)
    CDS_stdout, CDS_stderr = execute(cmd_line)

    output_file = open(output_name, 'w')
    output_file.write(CDS_stdout)
    output_file.close()

    # if CDS_stderr != '':
    #     err = True
    #     return None
    #
    # # return self._make_star_cat(output.decode("utf-8"))
    # return CDS_stdout


def main(input_dir, output_dir, kind):

    file_list = os.listdir(input_dir)

    for f in file_list:
        if 'image' not in f:
            continue

        # img = fits.open(input_dir + '/' + f)

        if kind == 'exp':
            list_ind = range(1, 41)
        else:
            list_ind = [0]

        for ind in list_ind:

            if kind == 'exp':
                exp_suff = '-{}'.format(ind-1)
            else:
                exp_suff = ''

            output_name = output_dir + '/star_cat' + re.split('image', os.path.splitext(f)[0])[1] + exp_suff + '.cat'

            if os.path.isfile(output_name):
                continue

            # h = img[ind].header
            h = fits.getheader(input_dir + '/' + f, ind)

            w = _get_wcs(h)

            img_shape = (h['NAXIS2'], h['NAXIS1'])
            img_center = np.array([img_shape[1]/2., img_shape[0]/2.])
            wcs_center = w.all_pix2world([img_center], 1)[0]
            astropy_center = SkyCoord(ra=wcs_center[0], dec=wcs_center[1], unit='deg')

            rad = _get_image_radius(img_center, w)

            find_stars(np.array([astropy_center.ra.value, astropy_center.dec.value]), output_name, rad)


if __name__ == '__main__':

    argv = sys.argv

    try:
        input_path = argv[1]
    except:
        input_path = '.'

    try:
        output_path = argv[2]
    except:
        output_path = '.'

    # Temporary fix
    try:
        kind = argv[3]
    except:
        kind = ''

    main(input_path, output_path, kind)

