"""SPLIT EXPOSURE.

Class to split single-exposure multi-CCD mosaic images into single-exposure
single-CCD files, one HDU per CCD.

This module splits the different CCDs (HDUs in FITS files) of a
single exposure into separate files.

:Author: Axel Guinot

"""

import numpy as np
import sip_tpv as stp
from astropy.io import fits
from astropy.wcs import WCS

from shapepipe.pipeline import file_io


class SplitExposures(object):
    """Split Exposures.

    Parameters
    ----------
    input_file_list : list
        Input file paths, typically image, weight, and flag
    output_dir : str
        Output directory
    file_number_string : str
        Input file identified
    output_suffix : str
        Output file suffix
    n_hdu : int
        Number of HDUs (CCDs) of input files

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        file_number_string,
        output_suffix,
        n_hdu
    ):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._file_number_string = file_number_string
        self._output_suffix = output_suffix
        self._n_hdu = n_hdu

    def process(self):
        """Process.

        Process the splitting of single-exposure images.

        """
        for exp_path, output_suffix in zip(
            self._input_file_list,
            self._output_suffix
        ):

            transf_int = 'flag' in output_suffix
            transf_coord = 'image' in output_suffix
            save_header = 'image' in output_suffix

            self.create_hdus(
                exp_path,
                output_suffix,
                transf_coord,
                transf_int,
                save_header
            )

    def create_hdus(
        self,
        exp_path,
        output_suffix,
        transf_coord,
        transf_int,
        save_header
    ):
        """Create HDUs.

        Split a single exposures CCDs into separate files.

        Parameters
        ----------
        exp_path : str
            Path to the single exposure
        output_suffix : str
            Suffix for the output file
        transf_coord : bool
            Transform the WCS (``pv`` to ``sip``) if ``True``
        transf_int : bool
            Set data types to int if ``True``
        save_header : bool
            Save WCS information if ``True``

        """
        # Initialise array to store headers of all HDUs
        header_file = np.zeros(self._n_hdu, dtype='O')

        hdu_list = fits.open(exp_path)

        # Check completeness of HDUs
        if len(hdu_list) <= self._n_hdu:
            raise TypeError(
                f'Only {len(hdu_list)} out of requested {self._n_hdu} HDUs '
                + f'found in image {exp_path}. File might be corrupt or not '
                'completely downloaded.'
            )

        # Loop over HDUs = CCDs
        for idx in range(1, self._n_hdu + 1):

            # Get HDU header
            h = hdu_list[idx].header

            # Transform between astrometry systems if requested
            if transf_coord:
                stp.pv_to_sip(h)

            # Get CCD image data
            d = hdu_list[idx].data

            # Transform (bool) to int type if requested
            if transf_int:
                d = d.astype(np.int16)

            # Set output file name
            file_name = (
                f'{self._output_dir}/{output_suffix}'
                + f'{self._file_number_string}-{str(idx-1)}.fits'
            )

            # Write HDU to output file
            new_file = file_io.FITSCatalogue(
                file_name,
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
            )
            new_file.save_as_fits(data=d, image=True, image_header=h)

            # Save header to array if requested
            if save_header:
                try:
                    w = WCS(h)
                except Exception:
                    print(f'WCS error for file {exp_path}')
                    raise
                header_file[idx - 1] = {'WCS': w, 'header': h.tostring()}

        # Save header array if requested
        if save_header:
            file_name = (
                f'{self._output_dir}/headers{self._file_number_string}.npy'
            )
            np.save(file_name, header_file)
