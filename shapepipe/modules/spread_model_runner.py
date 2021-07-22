# -*- coding: utf-8 -*-

"""SPREAD MODEL RUNNER

Module runner for ``spread_model``

:Author: Axel Guinot

:Date: 2019, 2020

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.spread_model_package import spread_model as sm


@module_runner(
    version='1.0',
    input_module=[
        'sextractor_runner',
        'psfex_interp_runner_me',
        'vignetmaker_runner'
    ],
    file_pattern=['sexcat', 'galaxy_psf', 'weight_vign'],
    file_ext=['.fits', '.sqlite', '.fits'],
    depends=['numpy', 'galsim'],
    run_method='parallel'
)
def spread_model_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    w_log
):

    # Get input files
    sex_cat_path, psf_cat_path, weight_cat_path = input_file_list

    # Get file suffix (optional)
    if config.has_option('SPREAD_MODEL_RUNNER', 'SUFFIX'):
        suffix = config.get('SPREAD_MODEL_RUNNER', 'SUFFIX')
        if (suffix.lower() != 'none') & (suffix != ''):
            suffix = suffix + '_'
        else:
            suffix = ''
    else:
        suffix = ''

    # Get pixel scale and output mode
    pixel_scale = config.getfloat('SPREAD_MODEL_RUNNER', 'PIXEL_SCALE')
    output_mode = config.get('SPREAD_MODEL_RUNNER', 'OUTPUT_MODE')

    # Set output file path
    file_name = f'{suffix}sexcat_sm{file_number_string}.fits'
    output_path = f'{run_dirs["output"]}/{file_name}'

    # Create spread model class instance
    inst = sm.SpreadModel(
        sex_cat_path,
        psf_cat_path,
        weight_cat_path,
        output_path,
        pixel_scale,
        output_mode
    )

    # Process spread model computation
    inst.process()

    return None, None
