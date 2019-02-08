# -*- coding: utf-8 -*-

"""TILEOBJ_AS_EXP_SCRIPT

This script contains a class to handle processing for the combine_mexp_script
module: Combine information on objects and PSF from multiple
exposure-single-CCD catalogues.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: February 2019

:Version: 1.0

"""


import os
import re

import numpy as np

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy import units as u
from astropy.table import Table

import shapepipe.pipeline.file_io as io


class CombineMExpError(Exception):
   """ CombineMExpError

   Generic error that is raised by this script.

   """

   pass


class combine_mexp():

    """Class for combining multiple exposure-single-CCD catalogue and PSF information.

    Parameters
    ----------
    input_file_name: list of 2 string
        path to input files, [image, catalogue]
    config: config class
        config file content
    output_dir: string
        output directory
    image_number: string
        image number for input tile according to numbering scheme from config file
    w_log: log file class
        log file

    Returns
    -------
    None
    """

    def __init__(self, input_file_names, output_dir, image_number, config, w_log):

        self._img_tile_path = input_file_names[0]
        self._cat_tile_path = input_file_names[1]
        self._output_dir = output_dir
        self._image_number = image_number
        self._config = config
        self._w_log = w_log

    def process(self):

        # List of exposures contributing to this tile
        exp_list_uniq = self.get_exposure_list()

        # List of exposure-single-CCD catalogue file names that were created by the
        # tileobj_as_exp modules
        exp_CCD_cat_list, PSF_list = self.get_exp_CCD_cat_and _PSF_list(len(exp_list_uniq))

        print(exp_CCD_cat_list)
        print(PSF_list)


    def get_exposure_list(self):
        """Return list of exposure file used for the tile in process, from tiles FITS header
           MKDEBUG: Note: this function is identical to the find_exposures class routine

        Parameters
        ----------

        Returns
        -------
        exp_list_uniq: list of string
            list of exposure basenames
        """

        try:
            hdu = fits.open(self._img_tile_path)
            hist = hdu[0].header['HISTORY']

        except:
            self._w_log.info('Error while reading tile catalogue FITS file {}, continuing...'.format(self._img_tile_path))

        tile_num = self._image_number

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise FindExposureError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile #{}'.format(len(exp_list_uniq), tile_num))
        self._w_log.info('{} duplicates were removed'.format(len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq

    def get_exp_CCD_cat_and_PSF_list(self, n_exp):
        """Return lists of exposure-single-CCD catalogue file base names and
           PSF vignet file base names.
           MKDEBUG note: This function is identical to the tileob_as_exp
           class routine.

        Parameters
        ----------
        n_exp: int
            number of exposures

        Returns
        -------
        exp_CCD_list: list of strings
            exposure-single-CCD catalogue base names
        PSF_list: list of strings
            PSF vignet file base names
        """

        n_hdu = int(self._config.get('TILEOBJ_AS_EXP_RUNNER', 'N_HDU'))
        exp_basename = 'cat.exp.img'
        PSF_basename = 'galaxy_psf'

        exp_CCD_list = []
        PSF_list     = []
        for i in range(n_exp):
            exp_CCD_base = '{}{}_{}'.format(exp_basename, self._image_number, i)
            PSF_base = '{}{}_{}'.format(PSF_basename, self._image_number, i)

            for k_img in range(n_hdu):
                exp_CCD_name = '{}_{:02d}'.format(exp_CCD_base, k_img)
                exp_CCD_list.append(exp_CCD_name)

                PSF_name = '{}_{:02d}'.format(PSF_base, k_img)
                PSF_list.append(PSF_name)

        return exp_CCD_list, PSF_list
