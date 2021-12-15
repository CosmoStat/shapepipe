# -*- coding: utf-8 -*-
"""MERGE SCRIPT

Class to merge separate catalogues.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2020

"""

import os
import re
import warnings

import numpy as np
from astropy.io import fits

from shapepipe.pipeline import file_io


class MergeSep(object):
    """Merge Sep
    Merge Separate Catalogues

    Parameters
    ----------
    input_file_list : list of str
        input file paths
    file_number_string : str
        file number following ShapePipe numbering scheme
    file_pattern : str
        file base name
    file_ext : str
        file extension
    output_dir : string
        output directory
    n_split_max : int
        number of separate input catalogues
    warning : str
        action when warning occurs, one in 'error', 'warning'
    w_log : logging.Logger
        log file
    """

    def __init__(
        self,
        input_file_list,
        file_number_string,
        file_pattern,
        file_ext,
        output_dir,
        n_split_max,
        warning,
        w_log
    ):

        self._input_file_list = input_file_list
        self._file_number_string = file_number_string
        self._file_pattern = file_pattern
        self._file_ext = file_ext
        self._output_dir = output_dir
        self._n_split_max = n_split_max
        self._warning = warning
        self._w_log = w_log

    def process(self):
        """Process

        Process merging of separate catalogues.
        """

        # Set warning action
        warnings.simplefilter(self._warning, UserWarning)

        # Loop over input files = outputs from different modules
        for idx, input_file in enumerate(self._input_file_list):

            # Get all input directories
            input_path_n = []
            input_path_n.append(input_file)
            for n in range(2, self._n_split_max + 1):
                res = re.sub('1', str(n), input_file, 1)
                input_path_n.append(res)

            # Open first catalogue, read number of extensions and columns
            cat0 = file_io.FITSCatalogue(input_file, SEx_catalogue=True)
            cat0.open()
            list_ext_name = cat0.get_ext_name()
            list_col_name = cat0.get_col_names()
            cat0.close()

            # Create empty dictionary
            # data dimension = n_extension x n_column x n_obj
            data = {}
            for hdu_ind, ext_name in enumerate(list_ext_name):
                if ext_name == 'PRIMARY':
                    continue
                data[ext_name] = {}
                for col_name in list_col_name:
                    data[ext_name][col_name] = []

            # Read and append all data, including first catalogue
            for n in range(self._n_split_max):
                cat_path = input_path_n[n]
                if os.path.exists(cat_path):
                    cat = file_io.FITSCatalogue(cat_path, SEx_catalogue=True)
                    cat.open()

                    for hdu_ind, ext_name in enumerate(list_ext_name):
                        if ext_name == 'PRIMARY':
                            continue
                        for col_name in list_col_name:
                            data[ext_name][col_name] += list(
                                cat.get_data(hdu_ind)[col_name]
                            )

                    cat.close()
                else:
                    msg = f'Input catalogue \'{cat_path}\' not found'
                    warnings.warn(msg)
                    wmsg = f'Warning: {msg}'
                    self._w_log.info(wmsg)
                    print(wmsg)

            # Save combined catalogue
            output_name = (
                f'{self._output_dir}/{self._file_pattern[idx]}'
                + f'{self._file_number_string}{self._file_ext[idx]}'
            )
            output = file_io.FITSCatalogue(
                output_name,
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
            )
            for hdu_ind, ext_name in enumerate(list_ext_name):
                if ext_name == 'PRIMARY':
                    continue
                output.save_as_fits(
                    data[ext_name],
                    names=list_col_name,
                    ext_name=ext_name
                )

        return None, None
