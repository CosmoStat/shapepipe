"""PASTE CATALOGUES.

Class to paste different (SExtractor) catalogs of objects.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Axel Guinot

"""

import os
import re

import numpy as np
from astropy import coordinates as coords
from astropy import units as u
from astropy.wcs import WCS

from shapepipe.pipeline import file_io


class PasteCat(object):
    """Paste Catalogue.

    Parameters
    ----------
    input_file_list : list of strings
        list of input catalogue paths to be pasted.
    output_path : str
        output file path of pasted catalog
    w_log : logging.Logger
        log file
    ext_name : list of str, optional, default = None
        HDU extension names, if None use input file names
    check_col_name : str, optional, default=None:
        if not None, use column with this key to check equal number
        of rows in each input catalog
    hdu_no : numpy.ndarray of int, optional, default = None
        hdu numbers of input catalog; by default set to 2 for all
        input files

    """

    def __init__(
        self,
        input_file_list,
        output_path,
        w_log,
        ext_name=None,
        check_col_name=None,
        hdu_no=None
    ):

        self._input_file_list = input_file_list
        self._output_path = output_path
        self._w_log = w_log
        self._ext_name = ext_name
        self._check_col_name = check_col_name
        if hdu_no is None:
            self._hdu_no = [2] * len(input_file_list)
        else:
            self._hdu_no = hdu_no

    def process(self):
        """Process.

        Process the pasting of catalogues.

        """
        # Create output catalog
        pasted_cat = file_io.FITSCatalogue(
            self._output_path,
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
        )

        for i, input_file in enumerate(self._input_file_list):

            self._w_log.info(f'Pasting catalogue \'{input_file}\'')

            # Read input data
            cat = file_io.FITSCatalogue(input_file)
            cat.open()
            data = np.copy(cat.get_data(self._hdu_no[i]))
            col_names = cat.get_col_names(self._hdu_no[i])
            cat.close()

            # Check equality if required by user
            if self._check_col_name:
                if i > 0:
                    if self._check_col_name not in col_names:
                        raise KeyError(
                            f'CHECK_COL_NAME key \'{self._check_col_name}\''
                            + 'not found in input catalog'
                        )
                    if not (
                        data[self._check_col_name]
                        == data_prev[self._check_col_name]
                    ).all():
                        raise Exception(
                            f'Column check using key'
                            + '\'{self._check_col_name}\''
                            + 'failed for input catalogs '
                            + f'#{i-1} and #{i}'
                        )
                data_prev = data

            # Add to output cat
            if self._ext_name:
                ext_name = self._ext_name[i]
            else:
                ext_name = input_file
            pasted_cat.save_as_fits(data, ext_name=ext_name)
