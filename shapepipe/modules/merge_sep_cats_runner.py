# -*- coding: utf-8 -*-

"""CREATE LOG EXP HEADER

This module merges output catalogues that have been created by separate (parallel)
calls to ShapePipe with the respective modules. Example: ngmix.

:Author: Morgan Schmitz, Axel Guinot, Martin Kilbinger

"""


import numpy as np
from astropy.io import fits

import os
import re

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
import shapepipe.pipeline.file_io as io


@module_runner(input_module='ngmix_runner', version='1.0',
               file_pattern=['ngmix'],
               file_ext=['.fits'], depends=['numpy'])
def merge_sep_cats_runner(input_file_list, run_dirs, file_number_string,
                          config, w_log):

    n_split_max = config.getint('MERGE_SEP_CATS_RUNNER', 'N_SPLIT_MAX')

    # Get all input directories
    input_path_n = []
    input_path_n.append(input_file_list[0])
    for n in range(2, n_split_max + 1):
        res = re.sub('1', str(n), input_file_list[0])
        input_path_n.append(res)

    # Open first catalogue, read number of extensions and columns
    cat0 = io.FITSCatalog(input_file_list[0], SEx_catalog=True)
    cat0.open()
    list_ext_name = cat0.get_ext_name()
    list_col_name = cat0.get_col_names()
    cat0.close()

    # Create empty dictionary
    # data: n_extension x n_column x n_obj
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
        cat = io.FITSCatalog(cat_path, SEx_catalog=True)
        cat.open()

        for hdu_ind, ext_name in enumerate(list_ext_name):
            if ext_name == 'PRIMARY':
                continue
            for col_name in list_col_name:
                data[ext_name][col_name] += list(cat.get_data(hdu_ind)[col_name])

        cat.close()


    # Save combined catalogue
    output_name = '{}/ngmix{}.fits'.format(run_dirs['output'], file_number_string)
    output = io.FITSCatalog(output_name,
                            open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    for hdu_ind, ext_name in enumerate(list_ext_name):
        if ext_name == 'PRIMARY':
            continue
        output.save_as_fits(data[ext_name], names=list_col_name, ext_name=ext_name)

    return None, None
