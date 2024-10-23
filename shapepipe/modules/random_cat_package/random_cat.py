"""RANDOM CATALOGUE.

This module contains a class to create a random catalogue, and to compute
the tile area accounting for overlapping and masked regions.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import os
import re

import numpy as np

import astropy.io.fits as fits
from astropy import wcs
from astropy.table import Table

from reproject import reproject_to_healpix

from shapepipe.pipeline import file_io
from shapepipe.utilities import cfis


class RandomCat:
    """Random Catalogue.

    This class creates a random catalogue given a mask FITS file.

    Parameters
    ----------
    input_image_path : str
        Path to input image file
    input_mask_path : str
        Path to input mask file
    output_dir : str
        Output directory
    file_number_pattern : str
        ShapePipe image ID string
    output_file_pattern : str
        Output file pattern (base name) for random catalogue
    n_rand : float
        Number of random objects on output
    density : bool
        ``n_rand`` is interpreted per square degrees if ``True``
    w_log : logging.Logger
        Logging instance
    healpix_options : dict
        Parameters for HEALPix output mask file
    """

    def __init__(
        self,
        input_image_path,
        input_mask_path,
        output_dir,
        file_number_string,
        output_file_pattern,
        n_rand,
        density,
        w_log,
        healpix_options,
    ):

        self._input_image_path = input_image_path
        self._input_mask_path = input_mask_path
        self._output_dir = output_dir
        self._file_number_string = file_number_string
        self._output_file_pattern = output_file_pattern
        self._n_rand = n_rand
        self._density = density
        self._w_log = w_log
        self._healpix_options = healpix_options

    def save_as_healpix(self, hdu_mask, header):
        """Save As Healpix.

        Save mask as healpix FITS file.

        Parameters
        ----------
        hdu_mask : class HDUList
            HDU with 2D pixel mask image
        header : class Header
            Image header with WCS information

        """
        if not self._healpix_options:
            return

        mask_1d, footprint = reproject_to_healpix(
            (hdu_mask, header),
            'galactic',
            nside=self._healpix_options['OUT_NSIDE']
        )

        t = Table()
        t['flux'] = mask_1d
        t.meta['ORDERING'] = 'RING'
        t.meta['COORDSYS'] = 'G'
        t.meta['NSIDE'] = self._healpix_options['OUT_NSIDE']
        t.meta['INDXSCHM'] = 'IMPLICIT'

        output_path = (
            f'{output_dir}/{self._healpix_options["FILE_BASE"]}-'
            + f'{file_number_string}.fits'
        )
        t.write(output_path)

    def process(self):
        """Process.

        Main function to identify exposures.

        """
        # Read image FITS file header
        try:
            img = fits.open(self._input_image_path)
            header = img[0].header
        except (OSError, IOError) as error:
            # FITS file might contain only header.
            # Try as ascii file
            try:
                fin = open(self._input_image_path)
                header = fits.Header.fromtextfile(fin)
                fin.close()
            except Exception:
                raise

        # Get WCS
        WCS = wcs.WCS(header)

        # Read mask FITS file
        hdu_mask = fits.open(self._input_mask_path)
        mask = hdu_mask[0].data

        # Save mask in healpix format (if option is set)
        self._save_as_healpix(hdu_mask, header)

        # Number of pixels
        n_pix_x = mask.data.shape[0]
        n_pix_y = mask.data.shape[1]
        n_pix = n_pix_x * n_pix_y

        # Number of non-masked pixels
        n_unmasked = len(np.where(mask == 0)[0])

        # Compute various areas

        # Pixel area in deg^2
        area_pix = wcs.utils.proj_plane_pixel_area(WCS)

        # Tile area
        area_deg2 = area_pix * n_pix

        # Area of unmasked region
        area_deg2_eff = area_pix * n_unmasked

        # Compute number of requested objects
        if n_unmasked > 0:
            if not self._density:
                # Use value from config file
                n_obj = self._n_rand
            else:
                # Compute number of objects from density
                n_obj = int(
                    self._n_rand / area_deg2 * area_deg2_eff / area_deg2
                )

            # Check that a reasonably large number of pixels is not masked
            if n_unmasked < n_obj:
                raise ValueError(
                    f"Number of un-masked pixels {n_unmasked} is smaller "
                    + f"than number of random objects requested {n_obj}"
                )

        else:
            n_obj = 0

        self._w_log.info(f"Creating {n_obj} random objects")

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
        if n_unmasked > 0:
            ra_rand = res[:, 0]
            dec_rand = res[:, 1]
            x_rand = xy_rand[:, 0]
            y_rand = xy_rand[:, 1]
        else:
            ra_rand = []
            dec_rand = []
            x_rand = []
            y_rand = []

        # Tile ID
        output_path = (
            f"{self._output_dir}/{self._output_file_pattern}-"
            + f"{self._file_number_string}.fits"
        file_base = os.path.splitext(file_name)[0]
        tile_ID_str = re.split("-", file_base)[1:]
        tile_id = float(".".join(tile_ID_str))
        tile_id_array = np.ones(n_obj) * tile_id

        # Write to output
        cat_out = [ra_rand, dec_rand, x_rand, y_rand, tile_id_array]
        column_names = ["RA", "DEC", "x", "y", "TILE_ID"]

        # TODO: Add units to header
        output = file_io.FITSCatalogue(
            output_path, open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
        )
        output.save_as_fits(cat_out, names=column_names)

        # Write area information to log file
        self._w_log.info(f"Total area = {area_deg2:.4f} deg^2")
        self._w_log.info(f"Unmasked area = {area_deg2_eff:.4f} deg^2")
        self._w_log.info(
            f"Ratio masked to total pixels = {n_unmasked / n_pix:.3f}"
        )
