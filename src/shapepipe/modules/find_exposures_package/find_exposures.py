"""FIND EXPOSURES SCRIPT.

This module contains a class to identify single exposures that were used
to create tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import re

import astropy.io.fits as fits


class FindExposures():
    """Find Exposures.

    This class finds exposures that are used for a given tile.

    Parameters
    ----------
    img_tile_path : str
        Path to tile image file
    output_path : str
        Output file path
    w_log : logging.Logger
        Log file
    colnum: int
        column number for exposure name in fits header
    prefix: str
        prefix for exposures
    """

    def __init__(self, img_tile_path, output_path, w_log, colnum, prefix):

        self._img_tile_path = img_tile_path
        self._output_path = output_path
        self._w_log = w_log
        self._colnum = colnum
        self.prefix = prefix

    def process(self):
        """Process.

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
        """Get Exposure List.

        Return list of exposure file used for the tile in process, from tiles
        FITS header.

        Returns
        -------
        list
            List of exposure basenames

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
        # History entries have format as the following example:
        # "input image 2243881p.fits 6 extension(s)"
        for _hist in hist:
            temp = _hist.split(' ')

            pattern = r'(.*)\.{1}.*'
            pattern_match = re.search(pattern, temp[self._colnum])
            if not pattern_match:
                raise IndexError(
                    f're match \'{pattern}\' failed for filename'
                    + f' \'{temp[self._colnum]}\''
                )

            exp_name = pattern_match.group(1)
            exp_name = exp_name.removeprefix(self.prefix)
            # LSB exposure names have 's', header still says 'p'
            # exp_name = re.sub(r'p', 's', exp_name)

            exp_list.append(exp_name)

        # Remove potential duplicates
        exp_list_uniq = list(set(exp_list))

        # For log output
        n_exp_uniq = len(exp_list_uniq)
        n_duplicate = len(exp_list) - n_exp_uniq
        self._w_log.info(f'Found {n_exp_uniq} exposures used in tile')
        self._w_log.info(f'{n_duplicate} duplicates were removed')

        return exp_list_uniq
