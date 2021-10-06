# -*- coding: utf-8 -*-

"""MCCD MERGE STARCAT RUNNER

This module is used to merge the validation stars resulting from the MCCD
validation runner or fit_validation runner.

:Author: Tobias Liaudat

"""

from shapepipe.modules.module_decorator import module_runner
import shapepipe.modules.mccd_merge_starcat_package.merge as merge


@module_runner(
    input_module=['mccd_fit_val_runner', 'mccd_val_runner'],
    version='1.0',
    file_pattern=['validation_psf'],
    file_ext=['.fits'],
    numbering_scheme='-0000000',
    depends=['numpy', 'astropy'],
    run_method='serial'
)
def mccd_merge_starcat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    w_log.info('Merging validation results.')

    output_dir = run_dirs['output']

    stamp_size = 51
    rad = 10
    ms = merge.MergeStarCat(input_file_list, output_dir, w_log, stamp_size, rad)

    ms.process()

    return None, None
