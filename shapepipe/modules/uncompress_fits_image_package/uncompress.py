# -*- coding: utf-8 -*-

"""UNCOMRESS FITS

This module uncompresses fits images and saves them as a single-hdu FITS file.


:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2020

"""

from astropy.io import fits


class Uncompress(object):
    """Uncompress

    This class handles the uncompress process of compressed FITS files.

    Parameters
    ----------
    input_file_list : array of string
        input files
    output_pattern_list : array of string
        output file pattern
    output_dir : string
        output directory
    file_number_string : string
        image numbering (pipeline ID string)
    data_hdu : int
        input HDU number
    """

    def __init__(
        self,
        input_file_list,
        output_pattern_list,
        output_dir,
        file_number_string,
        data_hdu
    ):

        self._input_file_list = input_file_list
        self._output_pattern_list = output_pattern_list
        self._output_dir = output_dir
        self._file_number_string = file_number_string
        self._data_hdu = data_hdu

    def process(self):
        """Process

        Main function to process uncompress.
        """

        # Go through input list
        for idx in range(len(self._input_file_list)):

            # Get data and header
            data = fits.getdata(self._input_file_list[idx], self._data_hdu)
            header = fits.getheader(self._input_file_list[idx], self._data_hdu)

            # Create and write new FITS file with that HDU only
            hdu = fits.PrimaryHDU(data, header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(
                f'{self._output_dir}/'
                + f'{self._output_pattern_list[idx]}'
                + f'{self._file_number_string}.fits'
            )
