# -*- coding: utf-8 -*-

"""MERGE STARCAT RUNNER

This module is used to merge the validation stars resulting from
validation runner or fit_validation runner.

:Author: Tobias Liaudat, Morgan Schmitz, Axel Guinot, Martin Kilbinger

"""

from shapepipe.modules.module_decorator import module_runner
import shapepipe.modules.merge_starcat_package.merge as merge


@module_runner(
    input_module=['mccd_fit_val_runner'],
    version='1.0',
    file_pattern=['validation_psf'],
    file_ext=['.fits'],
    numbering_scheme='-0000000',
    depends=['numpy', 'astropy'],
    run_method='serial'
)
def merge_starcat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    # Read config file options
    psf_model = config.get(module_config_sec, 'PSF_MODEL')
    allowed_psf_models = ['psfex', 'mccd']
    if psf_model not in allowed_psf_models:
        raise ValueError(
            f'Invalid config entry PSF_MODEL={psf_model} found,'
            + f'needs to be one of {allowed_psf_models}'
        )

    output_dir = run_dirs['output']

    if psf_model == 'mccd':
        MSC = merge.MergeStarCatMCCD
    else:
        MSC = merge.MergeStarCatPSFEX

    ms = MSC(input_file_list, output_dir, w_log)

    ms.process()

    return None, None
