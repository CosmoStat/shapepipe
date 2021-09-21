# -*- coding: utf-8 -*-

"""MERGE SEP CATS RUNNER

This module merges output catalogues that have been created by separate (parallel)
calls to ShapePipe with the respective modules. Example: ngmix.

:Author: Morgan Schmitz, Axel Guinot, Martin Kilbinger

"""


import numpy as np
from astropy.io import fits

import os
import re
import warnings

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
import shapepipe.pipeline.file_io as io


@module_runner(input_module='ngmix_runner', version='1.0',
               file_pattern=['ngmix'],
               file_ext=['.fits'], depends=['numpy'])
def merge_sep_cats_runner(input_file_list, run_dirs, file_number_string,
                          config, module_config_sec, w_log):

    n_split_max = config.getint(module_config_sec, 'N_SPLIT_MAX')

    file_pattern = config.getlist(module_config_sec, 'FILE_PATTERN')
    file_ext = config.getlist(module_config_sec, 'FILE_EXT')

    if config.has_option(module_config_sec, 'WARNING'):
        warning = config.get(module_config_sec, 'WARNING')
    else:
        warning = 'error'
    warnings.simplefilter(warning, UserWarning)

    # Loop over input files = outputs from different modules
    for idx in range(len(input_file_list)):

        # Get all input directories
        input_path_n = []
        input_path_n.append(input_file_list[idx])
        for n in range(2, n_split_max + 1):
            res = re.sub('1', str(n), input_file_list[idx], 1)
            input_path_n.append(res)

        # Open first catalogue, read number of extensions and columns
        cat0 = io.FITSCatalog(input_file_list[idx], SEx_catalog=True)
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
        for n in range(n_split_max):
            cat_path = input_path_n[n]
            if os.path.exists(cat_path):
                cat = io.FITSCatalog(cat_path, SEx_catalog=True)
                cat.open()

                for hdu_ind, ext_name in enumerate(list_ext_name):
                    if ext_name == 'PRIMARY':
                        continue
                    for col_name in list_col_name:
                        data[ext_name][col_name] += list(cat.get_data(hdu_ind)[col_name])

                cat.close()
            else:
                msg = 'Input catalogue \'{}\' not found'.format(cat_path)
                warnings.warn(msg)
                wmsg = 'Warning: {}'.format(msg)
                w_log.info(wmsg)
                print(wmsg)

        # Save combined catalogue
        output_name = '{}/{}{}{}'.format(run_dirs['output'], file_pattern[idx],
                                         file_number_string, file_ext[idx])
        output = io.FITSCatalog(output_name,
                                open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        for hdu_ind, ext_name in enumerate(list_ext_name):
            if ext_name == 'PRIMARY':
                continue
            output.save_as_fits(data[ext_name], names=list_col_name, ext_name=ext_name)

    return None, None
