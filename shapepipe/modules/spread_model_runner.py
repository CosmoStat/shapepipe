# -*- coding: utf-8 -*-

"""SPREAD MODEL RUNNER

This module computes the spread model.

:Author: Axel Guinot

:Date: 2019, 2020

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.SpreadModel_package import SpreadModel_script as sm



@module_runner(
    version='1.0',
    input_module=['sextractor_runner', 'psfex_interp_runner_me',
                  'vignetmaker_runner'],
    file_pattern=['sexcat', 'galaxy_psf', 'weight_vign'],
    file_ext=['.fits', '.sqlite', '.fits'],
    depends=['numpy', 'galsim'],
    run_method='parallel'
)
def spread_model_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):

    sex_cat_path, psf_cat_path, weight_cat_path = input_file_list

    if config.has_option('SPREAD_MODEL_RUNNER', 'SUFFIX'):
        suffix = config.get('SPREAD_MODEL_RUNNER', 'SUFFIX')
        if (suffix.lower() != 'none') & (suffix != ''):
            suffix = suffix + '_'
        else:
            suffix = ''
    else:
        suffix = ''

    pixel_scale = config.getfloat('SPREAD_MODEL_RUNNER', 'PIXEL_SCALE')
    output_mode = config.get('SPREAD_MODEL_RUNNER', 'OUTPUT_MODE')

    file_name = f'{suffix}sexcat_sm{file_number_string}.fits'
    output_path = f'{run_dirs["output"]}/{file_name}'

    inst = sm.SpreadModel(sex_cat_path, psf_cat_path, weight_cat_path,
                          output_path, pixel_scale, output_mode)

    inst.process()

    return None, None
