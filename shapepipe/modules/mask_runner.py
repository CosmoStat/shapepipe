# -*- coding: utf-8 -*-

"""MASK RUNNER

This file is the pipeline runner for the mask package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.mask_package.mask_script import mask


@module_runner(version='1.0', file_pattern=['image', 'weight', 'flag'],
               file_ext=['.fits','.fits','.fits'], depends=['numpy','astropy'],
               executes=['ww','findgsc2.2'])
def mask_runner(input_file_list, output_dir, file_number_string,
                   config, w_log):

    if len(input_file_list) == 2:
        ext_flag_name = None
    elif len(input_file_list) == 3:
        ext_flag_name = input_file_list[2]
    else:
        raise ValueError("Input files must be 'image', 'weight' and 'ext_flags' (optionnal)")
    
    config_file = config.get('MASK_RUNNER', 'MASK_CONFIG_PATH')
    suffix = config.get('MASK_RUNNER', 'SUFFIX')

    inst = mask(*input_file_list[:2], suffix.replace(" ",""), file_number_string, config_file, output_dir, path_external_flag=ext_flag_name)
    stdout, stderr = inst.make_mask()

    return stdout, stderr