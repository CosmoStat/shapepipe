# -*- coding: utf-8 -*-

"""PASTE CAT RUNNER

This module pastes different (SExtractor) catalogs of objects with identical IDs.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Axel Guinot

:Date: 10/2020

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io

import os
import re

import numpy as np

from astropy import units as u
from astropy import coordinates as coords
from astropy.wcs import WCS


class PasteCat(object):
    """PasteCat initialisation.

    Parameters
    ----------
    input_file_list : list of strings
        list of input catalogue paths to be pasted.
    output_path : string
        output file path of pasted catalog
    w_log :
        log file
    ext_name : list of strings, optional, default=None
        HDU extension names, if None use input file names
    check_col_name : string, optional, default=None:
        if not None, use column with this key to check equal number
        of rows in each input catalog
    hdu_no : array of int, optional, default=None
        hdu numbers of input catalog; by default set to 2 for all
        input files
    """

    def __init__(self, input_file_list, output_path, w_log,
                 ext_name=None, check_col_name=None, hdu_no=None):

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

        # Create output catalog
        pasted_cat = io.FITSCatalog(self._output_path,
                                    open_mode=io.BaseCatalog.OpenMode.ReadWrite)

        for i, input_file in enumerate(self._input_file_list):
            self._w_log.info('Pasting catalog \'{}\''.format(input_file))

            # Read input data
            cat = io.FITSCatalog(input_file)
            cat.open()
            data = np.copy(cat.get_data(hdu_no=self._hdu_no[i]))
            col_names = cat.get_col_names()
            cat.close()

            # Check equality
            if self._check_col_name:
                if i > 0:
                    if self._check_col_name not in col_names:
                        raise KeyError('CHECK_COL_NAME key \'{}\' not found in '
                                       'input catalog'
                                       ''.format(self._check_col_name))
                    if not (data[self._check_col_name] == data_prev[self._check_col_name]).all():
                        raise Exception('Column check using key \'{}\' failed for input catalogs '
                                        '#{} and #{}'
                                        ''.format(self._check_col_name, i-1, i))
                data_prev = data

            # Add to output cat
            if self._ext_name:
                ext_name = self._ext_name[i]
            else:
                ext_name = input_file
            pasted_cat.save_as_fits(data, ext_name=ext_name)


@module_runner(version='1.0',
               input_module='sextractor_runner',
               file_pattern='tile_sexcat',
               file_ext='.fits',
               depends='numpy',
               run_method='parallel')
def paste_cat_runner(input_file_list, run_dirs, file_number_string,
                     config, w_log):

    if config.has_option('PASTE_CAT_RUNNER', 'CHECK_COL_NAME'):
        check_col_name = config.get('PASTE_CAT_RUNNER', 'CHECK_COL_NAME')
    else:
        check_col_name = None

    if config.has_option('PASTE_CAT_RUNNER', 'HDU'):
        tmp = config.getlist('PASTE_CAT_RUNNER', 'HDU')
        hdu_no = [int(i) for i in tmp]
        if len(hdu_no) != len(input_file_list):
            raise IndexError('Different lengths for input file list ({}) and'
                             'HDU ({})'
                             ''.format(len(input_file_list), len(hdu_no)))
    else:
        hdu_no = None

    if config.has_option('PASTE_CAT_RUNNER', 'OUTPUT_FILE_PATTERN'):
        output_file_pattern = config.get('PASTE_CAT_RUNNER', 'OUTPUT_FILE_PATTERN')
    else:
        output_file_pattern = 'cat_pasted'

    if config.has_option('PASTE_CAT_RUNNER', 'EXT_NAME'):
        ext_name_list = config.getlist('PASTE_CAT_RUNNER', 'EXT_NAME')
        if len(ext_name_list) != len(input_file_list):
            raise ValueError('Input file list length ({}) and EXT_NAME list ({}) '
                             'need to be equal'
                             ''.format(len(input_file_list), len(ext_name_list)))
    else:
        ext_name_list = None

    file_ext = 'fits'

    output_path = '{}/{}{}.{}'.format(run_dirs['output'],
                                      output_file_pattern,
                                      file_number_string,
                                      file_ext)

    inst = PasteCat(input_file_list, output_path, w_log, ext_name=ext_name_list,
                    check_col_name=check_col_name, hdu_no=hdu_no)

    inst.process()

    return None, None
