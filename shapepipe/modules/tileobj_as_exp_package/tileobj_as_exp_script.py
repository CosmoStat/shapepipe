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

        wcs_list, nx_list, ny_list = self.get_wcs_list(exp_CCD_list)

        dat_tile = self.get_tile_cat_data()

        n_tileobj = self.write_tileobj_to_exp_CCDs(dat_tile, exp_CCD_list, wcs_list, nx_list, ny_list)

        self._w_log.info('Number of exposure-single-CCD images on input: {}'.format(len(exp_CCD_list)))
        self._w_log.info('Number of headers with WCS information: {}'.format(len(wcs_list)))
        self._w_log.info('Number of tile-object catalogues in expoure-single-CCD format written: {}'.format(n_tileobj))



    def write_tileobj_to_exp_CCDs(self, dat_tile, exp_CCD_list, wcs_list, nx_list, ny_list):
        """Write objects from tile to exposure-single-CCD catalogues

        Parameters
        ----------
        dat_tile: astropy.table.Table
            tile catalogue data
        exp_CCD_list: list of strings
            exposure-single-CCD file name base list
        wcs_list: dictionary of wcs.WCS
            WCS information
        nx_list: dictionary of int
            number of pixels in x
        ny_list: dictionary of int
            number of pixels in y

        Returns
        -------
        n_tileobj: int
            number of tile-object catalogues written
        """

        out_base = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'OUTPUT_FILE_PATTERN')
        ext_out  = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'OUTPUT_FILE_EXT')
        sex_cat_template = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'SEX_CAT_TEMPLATE')

        ### Set columns to be returned in catalogue data

        # Basic columns from FITS file

        # MKDEBUG TODO: Get from dat_tile
        cols = dat_tile.keys()

        # Assumes all but last column are float, and last (ID) is int.
        # I couldn't extract this information from dat_tile.
        dt = [float for c in cols[:-1]]
        dt.append(int)


        all_coord_tile_wcs = SkyCoord(ra=dat_tile['X_WORLD']*u.degree, dec=dat_tile['Y_WORLD']*u.degree)

        n_tileobj = 0
        for exp_CCD in exp_CCD_list:

            # Tile object coordinates in x, y using exposure-CCD WCS
            all_coord_tile_xy  = wcs_list[exp_CCD].all_world2pix(all_coord_tile_wcs.ra, all_coord_tile_wcs.dec, 0)

            # Find objects that are within this CCD
            ind_in_range       = ((all_coord_tile_xy[0] >= 0) & (all_coord_tile_xy[0] < nx_list[exp_CCD]) & \
                                  (all_coord_tile_xy[1] >= 0) & (all_coord_tile_xy[1] < ny_list[exp_CCD]))

            if ind_in_range.any():

                # Create temporary table with image coordinates from expoure, world coordinates from tile
                tileobj = Table([all_coord_tile_xy[0][ind_in_range],
                                 all_coord_tile_xy[1][ind_in_range],
                                 all_coord_tile_wcs.ra[ind_in_range],
                                 all_coord_tile_wcs.dec[ind_in_range],
                                 dat_tile['FWHM_IMAGE'][ind_in_range],
                                 dat_tile['ID'][ind_in_range]],
                                names=cols, dtype=dt)

                # Write to FITS file
                output_path  = '{}/{}{}.{}'.format(self._output_dir, out_base, exp_CCD, ext_out)
                exp_cat_file = io.FITSCatalog(output_path, open_mode=io.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=True)
                #import pdb
                #pdb.set_trace()
                exp_cat_file.save_as_fits(data=tileobj, names=cols, ext_name='LDAC_OBJECTS', sex_cat_path=sex_cat_template)

                n_tileobj = n_tileobj + 1

        return n_tileobj



    def get_tile_cat_data(self):
        """Return data from tiles catalogue.

        Parameters
        ----------
        None

        Returns
        -------
        dat_tile: astropy.table.Table
            catalogue data
        """

        ### Set columns to be returned in catalogue data

        # Basic columns from FITS file
        cols = ('X_IMAGE', 'Y_IMAGE', 'X_WORLD', 'Y_WORLD', 'FWHM_IMAGE')
        dt   = [float for c in cols]

        f_tile = io.FITSCatalog(self._cat_tile_path, SEx_catalog=True)
        f_tile.open()
        tmp = f_tile.get_data()

        # Use only columns given above
        dat_tile = Table([tmp[:][c] for c in cols], names=cols, dtype=dt)

        # Add column with galaxy IDs (running object number)
        dat_tile['ID'] = np.arange(0, len(dat_tile[cols[0]]))

        return dat_tile


    def get_wcs_list(self, exp_CCD_list):
        """Return list of WCS corresponding to exposure-single-CCD images.

        Parameters
        ----------
        exp_CCD_list: list of strings
            list of exposure-single-CCD names

        Returns
        -------
        wcs_list: dictionary of WCS objects
            list of WCS information
        """

        input_dir_exp_CCD = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'INPUT_DIR_EXP_CCD')
        ext               = self._config.get('TILEOBJ_AS_EXP_RUNNER', 'INPUT_EXP_CCD_EXT')

        nx_list = {}
        ny_list = {}
        wcs_list = {}
        for exp_CCD in exp_CCD_list:
            exp_CCD_path = '{}/{}.{}'.format(input_dir_exp_CCD, exp_CCD, ext)

            if os.path.isfile(exp_CCD_path):
                exp_CCD_img = io.FITSCatalog(exp_CCD_path, hdu_no=0)
                exp_CCD_img.open()
                header         = exp_CCD_img.get_header(hdu_no=0)
                ny_CCD, nx_CCD = exp_CCD_img.get_data().shape
                exp_CCD_img.close()

                wcs_list[exp_CCD] = wcs.WCS(header)
                nx_list[exp_CCD]  = nx_CCD
                ny_list[exp_CCD]  = ny_CCD

            else:
                msg = 'Exposure-single-CCD file \'{}\' not found, continuing'.format(exp_CCD_path)
                self._w_log.info(msg)

        return wcs_list, nx_list, ny_list



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
        basename = 'exp.img'

        exp_CCD_list = []
        for i in range(n_exp):
            exp_CCD_base = '{}{}_{}'.format(basename, self._image_number, i)

            for k_img in range(n_hdu):
                exp_CCD_name = '{}_{:02d}'.format(exp_CCD_base, k_img)
                exp_CCD_list.append(exp_CCD_name)

        return exp_CCD_list


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
            hdu   = fits.open(self._img_tile_path)
            hist  = hdu[0].header['HISTORY']

        except:
            self._w_log.info('Error while reading tile catalogue FITS file {}, continuing...'.format(self._img_tile_path))

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


