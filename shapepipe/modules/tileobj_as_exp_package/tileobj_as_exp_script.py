# -*- coding: utf-8 -*-

"""TILEOBJ_AS_EXP_SCRIPT

This script contains a class to handle processing for the tileobj_as_exp_script module:
Writing objects selected on tiles to catalogue files in expoure-single-CCD format.

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



class TileObjAsExpError(Exception):
   """ TileObjAs_ExpError

   Generic error that is raised by this script.

   """

   pass



class tileobj_as_exp():

    """Class for writing selected objects on tiles to catalogues in
        exposure-single-CCD format

    Parameters
    ----------
    cat_tile_path: string
        path to tile catalogue file
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

    def __init__(self, cat_tile_path, output_dir, image_number, config, w_log):

        self._cat_tile_path = cat_tile_path
        self._output_dir    = output_dir
        self._image_number  = image_number
        self._config        = config
        self._w_log         = w_log



    def process(self):

        # List of exposures contributing to this tile
        exp_list_uniq = self.get_exposure_list()
        #print(exp_list_uniq)

        # List of exposure-single-CCD file names that were created by the
        # find_exposures modules
        exp_CCD_list = self.get_exp_CCD_list(len(exp_list_uniq))
        #print(exp_CCD_list)

        wcs_list = self.get_wcs_list(exp_CCD_list)

        self._w_log.info('Processed {} exposure-single-CCD catalogues'.format(len(exp_CCD_list)))



    def get_wcs_list(self, exp_CCD_list):
        """Return list of WCS corresponding to exposure-single-CCD images.

        Parameters
        ----------
        exp_CCD_list: list of strings
            list of exposure-single-CCD names

        Returns
        -------
        wcs_list: list of WCS objects
            list of WCS information
        """

        input_dir_exp_CCD = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'INPUT_DIR_EXP_CCD')
        ext_out           = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'OUTPUT_FILE_EXT')

        wcs_list = []
        for exp_CCD in exp_CCD_list:
            exp_CCD_path = '{}/{}.{}'.format(input_dir_exp_CCD, exp_CCD, ext_out)

            exp_CCD_img = io.FITSCatalog(exp_CCD_path, hdu_no=0)
            exp_CCD_img.open()
            header      = exp_img.get_header(hdu_no=0)
            exp_CCD_img.close()

            wcs_list.append(wcs.WCS(header))

        return wcs_list



    def get_exp_CCD_list(self, n_exp):
        """Return a list of exposure-single-CCD names.

        Parameters
        ----------
        n_exp: int
            number of exposures

        Returns
        -------
        exp_CCD_list: list of strings
            list of exposure-single-CCD names
        """

        n_hdu    = int(self._config.get('TILEOBJ_AS_EXP_RUNNER', 'N_HDU'))
        basename = 'exp'

        exp_CCD_list = []
        for i in range(n_exp):
            exp_CCD_base = '{}{}_{}'.format(basename, self._image_number, i)

            for k_img in range(n_hdu):
                exp_CCD_name = '{}_{:02d}'.format(exp_CCD_base, k_img)
                exp_CCD_list.append(exp_CCD_name)

        return exp_CCD_list


    def get_exposure_list(self):
        """Return list of exposure file used for the tile in process, from tiles FITS header

        Parameters
        ----------

        Returns
        -------
        exp_list_uniq: list of string
            list of exposure basenames
        """

        try:
            hdu   = fits.open(self._cat_tile_path)
            hist  = hdu[0].header['HISTORY']

        except:
            self._w_log.info('Error while reading tile catalogue FITS file {}, continuing...'.format(self._cat_tile_path))

        # MKDEBUG: This has to be checked, here we only want the input tile number.
        tile_num = self._image_number

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp  = h.split(' ')

            pattern = '(.*)\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise FindExposureError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile #{}'.format(len(exp_list_uniq), tile_num))
        self._w_log.info('{} duplicates were removed'.format(len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq


