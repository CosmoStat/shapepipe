# -*- coding: utf-8 -*-

"""CREATE LOG EXP HEADER

This module merge the "header" files output of the split_exp_runner.py module.
It create a binnary file that contain the wcs of each CCDs for each exposures.

:Author: Axel Guinot

"""


import os
import re

import numpy as np
from sqlitedict import SqliteDict

from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='split_exp_runner', version='1.0',
               file_pattern=['header'],
               file_ext=['.npy'], depends=['numpy', 'sqlitedict'],
               run_method='serial')
def merge_headers_runner(input_file_list, run_dirs, file_number_string,
                         config, w_log):

    output_dir = run_dirs['output']
    if config.has_option('MERGE_HEADER_RUNNER', 'OUTPUT_PATH'):
        output_dir = config.getexpanded('MERGE_HEADER_RUNNER', 'OUTPUT_PATH')

    final_file = SqliteDict(output_dir + '/log_exp_headers.sqlite')
    for file_path in input_file_list:
        key = re.split('headers-', os.path.splitext(os.path.split(file_path[0])[1])[0])[1]
        final_file[key] = np.load(file_path[0], allow_pickle=True)

    final_file.commit()
    final_file.close()

    return None, None
