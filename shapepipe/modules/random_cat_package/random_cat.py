# -*- coding: utf-8 -*-
"""RANDOM CAT SCRIPT

This module contains a class to create a random catalogue, and to compute
the tile area accounding for overlapping and masked regions.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: December 2021

"""


import re
import os

import numpy as np
import astropy.io.fits as fits
from astropy import wcs

import shapepipe.pipeline.file_io as io
from shapepipe.utilities import cfis


class RandomCat():
    """Random Catalogue

    This class creates a random catalogue given a mask FITS file

    Parameters
    ----------
    input_image_path : str
        path to input image file
    input_mask_path : str
        path to input mask file
    output_path : str
        output file path for random catalogue
    n_rand : float
        number of random objects on output
    density : bool
        n_rand is interpreted per square degrees if True
    w_log : logging.Logger
        log file
    tile_list_path : str, optional, default=None
        List to all tile IDs, to remove objects in
        overlapping tile areas
    """

    def __init__(
        self,
        input_image_path,
        input_mask_path,
        output_path,
        n_rand,
        density,
        w_log,
        tile_list_path=None
    ):

        self._input_image_path = input_image_path
        self._input_mask_path = input_mask_path
        self._output_path = output_path
        self._n_rand = n_rand
        self._density = density
        self._w_log = w_log
        self._tile_list_path = tile_list_path

    def process(self):
        """Process

        Main function to identify exposures.
        """

        # Read image FITS file header
        try:
            img = fits.open(self._input_image_path)
            header = img[0].header
        except OSError:
        except IOError:
            # FITS file might contain only header.
            # Try as ascii file
            try:
                fin = open(self._input_image_path)
                header = fits.Header.fromtextfile(fin)
                fin.close()
            except:
                raise

        # Get WCS
        WCS = wcs.WCS(header)

        # Read mask FITS file
        hdu_mask = fits.open(self._input_mask_path)
        mask = hdu_mask[0].data

        # Number of pixels
        n_pix_x = mask.data.shape[0]
        n_pix_y = mask.data.shape[1]
        n_pix = n_pix_x * n_pix_y

        # Number of non-masked pixels
        n_unmasked = len(np.where(hdu_mask[0].data == 0)[0])

        # Compute various areas

        # Pixel area in deg^2
        area_pix = wcs.utils.proj_plane_pixel_area(WCS)

        # Tile area
        area_deg2 = area_pix * n_pix

        # Area of unmasked region
        area_deg2_eff = area_pix * n_unmasked

        # Compute number of requested objects
        if not self._density:
            # Use value from config file
            n_obj = self._n_rand
        else:
            # Compute number of objects from density
            n_obj = int(self._n_rand / area_deg2 * area_deg2_eff / area_deg2)

        self._w_log.info(f'Creating {n_obj} random objects')

        # Check that a reasonably large number of pixels is not masked
        if n_unmasked < n_obj:
            raise ValueError(
                f'Number of un-masked pixels {n_unmasked} is '
                + f'smaller than number of random objects requested {n_obj}'
            )

        # Draw points until n are in mask
        n_found = 0
        xy_rand = []
        while n_found < n_obj:
            idx_x = np.random.randint(n_pix_x)
            idx_y = np.random.randint(n_pix_y)

            # Add points with additional random sub-pixel value
            if mask[idx_x, idx_y] == 0:
                d = np.random.random(2)
                # MKDEBUG: the following seems to work, x and y interchanged
                xy_rand.append([idx_y + d[1], idx_x + d[0]])
                n_found = n_found + 1
        xy_rand = np.array(xy_rand)

        # Transform to WCS
        res = WCS.all_pix2world(xy_rand, 1)
        ra_rand = res[:, 0]
        dec_rand = res[:, 1]

        # Tile ID
        file_name = os.path.split(self._output_path)[1]
        file_base = os.path.splitext(file_name)[0]
        tile_ID_str = re.split('-', file_base)[1:]
        tile_id = float('.'.join(tile_ID_str))
        tile_id_array = np.ones(n_obj) * tile_id

        # Write to output
        cat_out = [
            ra_rand,
            dec_rand,
            xy_rand[:, 0],
            xy_rand[:, 1],
            tile_id_array
        ]
        column_names = ['RA', 'DEC', 'x', 'y', 'TILE_ID']

        # TODO: Add units to header
        output = io.FITSCatalog(
            self._output_path,
            open_mode=io.BaseCatalog.OpenMode.ReadWrite
        )
        output.save_as_fits(cat_out, names=column_names)

        # Remove overlapping regions
        if self._tile_list_path:
            self._w_log.info('Flag overlapping objects')
            ratio_non_overl_tot = cfis.remove_common_elements(
                output,
                self._tile_list_path,
                pos_param=['RA', 'DEC']
            )
        else:
            ratio_non_overl_tot = 1

        # Compute area without overlapping regions (approximation)
        area_deg2_non_overl = area_deg2 * ratio_non_overl_tot
        area_deg2_eff_non_overl = area_deg2_eff * ratio_non_overl_tot

        # Write area information to log file
        self._w_log.info(f'Total area = {area_deg2:.4f} deg^2')
        self._w_log.info(f'Unmasked area = {area_deg2_eff:.4f} deg^2')
        self._w_log.info(
            f'Ratio masked to total pixels = {n_unmasked / n_pix:.3f}'
        )

        self._w_log.info(
            f'Total area without overlap = {area_deg2_non_overl:.4f} deg^2'
        )
        self._w_log.info(
            'Unmaskewd area without overlap = '
            + f'{area_deg2_eff_non_overl:.4f} deg^2'
        )
        self._w_log.info(
            f'Ratio of non-overlap to total area = {ratio_non_overl_tot:.3f}'
        )
