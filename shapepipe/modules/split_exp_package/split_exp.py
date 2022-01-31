"""SPLIT EXPOSURE.

Class to split single-exposure multi-CCD mosaic images into single-exposure
single-CCD files, one HDU per CCD.

This module splits the different CCDs (= hdus in FITS files) of a
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
    input_file_list : list of str
        input file paths, typically image, weight, and flag
    output_dir : str
        output directory
    file_number_string : str
        input file identified
    output_suffix : str
        output file suffix
    n_hdu : int
        number of HDUs (= CCDs) of input files

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
            Transform the WCS (pv to sip) if True
        transf_int : bool
            Set data types to int if True
        save_header : bool
            Save WCS information if True

        """
        header_file = np.zeros(self._n_hdu, dtype='O')

        for idx in range(1, self._n_hdu + 1):

            h = fits.getheader(exp_path, idx)
            if transf_coord:
                stp.pv_to_sip(h)

            d = fits.getdata(exp_path, idx)
            if transf_int:
                d = d.astype(np.int16)

            file_name = (
                f'{self._output_dir}/{output_suffix}'
                + f'{self._file_number_string}-{str(idx-1)}.fits'
            )

            new_file = file_io.FITSCatalogue(
                file_name,
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
            )
            new_file.save_as_fits(data=d, image=True, image_header=h)

            if save_header:
                try:
                    w = WCS(h)
                except Exception:
                    print(f'WCS error for file {exp_path}')
                    raise
                header_file[idx - 1] = {'WCS': w, 'header': h.tostring()}

        if save_header:
            file_name = (
                f'{self._output_dir}/headers{self._file_number_string}.npy'
            )
            np.save(file_name, header_file)
