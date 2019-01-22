# -*- coding: utf-8 -*-

"""FIND_EXPOSURES SCRIPT

This script contains a class to handle processing for the find_exposures modules: Identify exposures that are used in selected tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

:Version: 1.0

"""


import os
import re

import numpy as np
import astropy.io.fits as fits
from astropy import wcs
import sip_tpv as stp

import shapepipe.pipeline.file_io as sc



class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)



class find_exposures():

    """Class for identifying exposures that are used for a given tile.

    Parameters
    ----------
    cat_tile_path: string
        path to tile catalogue file
    config: config class
        config file content
    output_dir: string
        output directory
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

        # Input file basename
        ext_arr = config.get('FILE', 'FILE_EXT')
        ext     = ext_arr.split(',')[0]
        fbasen  = os.path.basename(cat_tile_path)
        m = re.search('(.*).{}'.format(ext), fbasen)
        if not m:
            raise MyError('Extension \'{}\' does not match tile cataog basename \'{}\''.format(ext, fbasen))
        self._cat_tile_basename = m.group(1)



    def process(self):

        exp_list_uniq = self.get_exposure_list()

        input_dir_exp = self._config.getlist('FIND_EXPOSURES', 'INPUT_DIR_EXP')
        num = []
        for i in range(len(input_dir_exp)):
            num_idx = self.create_exposures(exp_list_uniq, i)
            num.append(num_idx)

        for i in range(len(input_dir_exp)):
            self._w_log.info('Created {} exposure files of type {}'.format(num[i], i))



    def get_exposure_list(self):
        """Read list of exposure file used for the tile in process

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
                raise MyError('re match \'{}\' failed for filename \'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile #{}'.format(len(exp_list_uniq), tile_num))
        self._w_log.info('{} duplicates were removed'.format(len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq


    def create_hdus(self, num, exp_path, exp_base, transf_coord=True, transf_int=False):
        """Create FITS exposure files, one for each input HDU.

        Parameters
        ----------
        self: class find_exposure
            this class
        num: int
            current number of written files
        exp_path: string
            input exposure multi-HDU file path
        exp_base: string
            output file base name
        transf_coord: bool, optional, default=True
            if True, transform coordinates to astropy format
        transf_int: bool, optional, default=False
            if True, transform data to integer (e.g. flags)

        Returns
        -------
        dnum: int
            number of files written
        """

        img_file = sc.FITSCatalog(exp_path, hdu_no=1)
        img_file.open()

        ext_out = self._config.get('FIND_EXPOSURES', 'OUTPUT_FILE_EXT')

        hdu_max = len(img_file._cat_data)
        n_hdu   = int(self._config.get('FIND_EXPOSURES', 'N_HDU'))

        dnum = 0
        for k_img in range(1, n_hdu+1):

            h_img = img_file._cat_data[k_img].header
            coord_img = re.findall(r"[\w]+", h_img['DETSEC'])

            # TODO optional: Compare image, weight, and flag coordinates to be sure
            # the HDU's match

            # Change coordinates to astropy-readable format
            if transf_coord:
                stp.pv_to_sip(h_img)

            d = img_file._cat_data[k_img].data
            if transf_int:
                d = d.astype(np.int16)
            new_fits = fits.PrimaryHDU(data=d, header=h_img)

            out_name = '{}/{}-hdu{:02d}.{}'.format(self._output_dir, exp_base, dnum, ext_out)
            new_fits.writeto(out_name)

            dnum = dnum + 1

        img_file.close()

        return dnum


    def create_exposures(self, exp_list, idx_type):
        """Create exposure FITS files.

        Parameters
        ----------
        self: class find_exposure
            this class
        exp_list: list of string
            list of exposure basenames
        idx_type: int
            index for input image type, used for input dir and file pattern

        Returns
        -------
        num: int
            number of files written
        """

        input_dir_exp = self._config.getlist('FIND_EXPOSURES', 'INPUT_DIR_EXP')[idx_type]
        input_pattern = self._config.getlist('FIND_EXPOSURES', 'INPUT_FILE_PATTERN')[idx_type]
        ext_in = self._config.get('FIND_EXPOSURES', 'INPUT_FILE_EXT')

        num = 0
        for exp_name in exp_list:

            exp_path = '{}/{}{}.{}'.format(input_dir_exp, exp_name, input_pattern, ext_in)

            exp_base_new = '{}{}-{}'.format(exp_name, input_pattern, self._cat_tile_basename)

            # Data-specific options

            # Transform coordinates of image header
            transf_coord = (idx_type==0)

            # Transform flag to integer
            transf_int   = (idx_type==2)

            dnum = self.create_hdus(num, exp_path, exp_base_new, transf_coord=transf_coord, transf_int=transf_int)
            num  = num + dnum

        return num

