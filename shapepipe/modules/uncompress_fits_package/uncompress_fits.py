"""UNCOMPRESS FITS.

This module uncompresses FITS images and saves them as a single-HDU FITS file.

:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from astropy.io import fits


class Uncompress(object):
    """Uncompress.

    This class handles the uncompress process of compressed FITS files.

    Parameters
    ----------
    input_file_list : numpy.ndarray
        Input files
    output_pattern_list : numpy.ndarray
        Output file pattern
    output_dir : str
        Output directory
    file_number_string : str
        Image numbering (pipeline ID string)
    data_hdu : int
        Input HDU number

    """

    def __init__(
        self,
        input_file_list,
        output_pattern_list,
        output_dir,
        file_number_string,
        data_hdu,
    ):

        self._input_file_list = input_file_list
        self._output_pattern_list = output_pattern_list
        self._output_dir = output_dir
        self._file_number_string = file_number_string
        self._data_hdu = data_hdu

    def process(self):
        """Process.

        Main function to process uncompressing.

        """
        # Go through input list
        for idx in range(len(self._input_file_list)):

            # Get data and header
            try:
                data = fits.getdata(self._input_file_list[idx], self._data_hdu)
            except Exception as e:
                print(
                    "Error while reading data from FITS image "
                    + f"{self._input_file_list[idx]} at HDU {self._data_hdu}"
                )
                raise
            header = fits.getheader(self._input_file_list[idx], self._data_hdu)

            # Create and write new FITS file with that HDU only
            hdu = fits.PrimaryHDU(data, header)
            hdul = fits.HDUList([hdu])
            hdul.writeto(
                f"{self._output_dir}/"
                + f"{self._output_pattern_list[idx]}"
                + f"{self._file_number_string}.fits"
            )
