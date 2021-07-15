# -*- coding: utf-8 -*-
"""SPLIT EXP SCRIPT

This module contains a class to identify single exposures that were used
to create tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

:Package: ShapePipe

"""


import re

import astropy.io.fits as fits

import shapepipe.pipeline.file_io as io


class FindExposures():
    """Find Exposures

    This class finds exposures that are used for a given tile.

    Parameters
    ----------
    img_tile_path : string
        path to tile image file
    output_path : string
        output file path
    w_log :
        log file
    """

    def __init__(self, img_tile_path, output_path, w_log):

        self._img_tile_path = img_tile_path
        self._output_path = output_path
        self._w_log = w_log

    def process(self):
        """Process

        Main function to identify exposures.
        """

        # Get list of exposures
        exp_list_uniq = self.get_exposure_list()

        # Write list to output ascii file
        f_out = open(self._output_path, 'w')
        if len(exp_list_uniq) > 0:
            for exp in exp_list_uniq:
                print(exp, file=f_out)
        f_out.close()

    def get_exposure_list(self):
        """Get Exposure List

        Return list of exposure file used for the tile in process, from tiles
        FITS header

        Parameters
        ----------

        Returns
        -------
        exp_list_uniq: list of string
            list of exposure basenames
        """

        try:
            # Get history from tiles FITS header
            hdu = fits.open(self._img_tile_path)
            hist = hdu[0].header['HISTORY']

        except Exception:
            # Key word not found -> raise error
            self._w_log.info(
                'Error while reading tile image FITS file '
                + f'{self._img_tile_path}, continuing...'
            )

        exp_list = []

        # Get exposure file names
        for h in hist:
            temp = h.split(' ')
