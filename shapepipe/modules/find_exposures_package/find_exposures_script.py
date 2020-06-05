# -*- coding: utf-8 -*-

"""FIND_EXPOSURES SCRIPT

This script contains a class to handle processing for the find_exposures
module: Identify exposures that are used in selected tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: June 2020

:Version: 1.1

"""


import re

import astropy.io.fits as fits

import shapepipe.pipeline.file_io as io


class find_exposures():

    """Class for identifying exposures that are used for a given tile.

    Parameters
    ----------
    img_tile_path: string
        path to tile image file
    config: config class
        config file content
    output_dir: string
        output directory
    image_number: string
        image number according to numbering scheme from config file
    w_log: log file class
        log file

    Returns
    -------
    None
    """

    def __init__(self, img_tile_path, output_dir, image_number, config, w_log):

        self._img_tile_path = img_tile_path
        self._output_dir = output_dir
        self._image_number = image_number
        self._config = config
        self._w_log = w_log

    def process(self):

        types = ['image', 'weight', 'flag']

        exp_list_uniq = self.get_exposure_list()

        if len(exp_list_uniq) > 0:
            for exp in exp_list_uniq:
                print(exp)

        return

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
            hdu = fits.open(self._img_tile_path)
            hist = hdu[0].header['HISTORY']

        except Exception:
            self._w_log.info('Error while reading tile image FITS file '
                             '{}, continuing...'.format(self._img_tile_path))

        tile_num = self._image_number

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp = h.split(' ')

            pattern = r'(.*)\.{1}.*'
            m = re.search(pattern, temp[3])
            if not m:
                raise IndexError('re match \'{}\' failed for filename '
                                 '\'{}\''.format(pattern, temp[3]))

            exp_name = m.group(1)

            exp_list.append(exp_name)

        exp_list_uniq = list(set(exp_list))

        self._w_log.info('Found {} exposures used in tile #{}'.format(
                         len(exp_list_uniq), tile_num))
        self._w_log.info('{} duplicates were removed'.format(
                         len(exp_list) - len(exp_list_uniq)))

        return exp_list_uniq
