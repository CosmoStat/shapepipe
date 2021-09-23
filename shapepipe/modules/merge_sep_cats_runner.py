# -*- coding: utf-8 -*-

"""MERGE SEP CATS RUNNER

This module merges output catalogues that have been created by separate (parallel)
calls to ShapePipe with the respective modules. Example: ngmix.

:Author: Morgan Schmitz, Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2020, 2021

:Package: ShapePipe

"""


from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.merge_sep_cats_package import merge_script as merge


@module_runner(
    input_module='ngmix_runner',
    version='1.0',
    file_pattern=['ngmix'],
    file_ext=['.fits'],
    depends=['numpy'])
def merge_sep_cats_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    # Get config entries
    n_split_max = config.getint(module_config_sec, 'N_SPLIT_MAX')

    file_pattern = config.getlist(module_config_sec, 'FILE_PATTERN')
    file_ext = config.getlist(module_config_sec, 'FILE_EXT')

    if config.has_option(module_config_sec, 'WARNING'):
        warning = config.get(module_config_sec, 'WARNING')
    else:
        warning = 'error'

    inst = merge.MergeSep(
        input_file_list,
        file_number_string,
        file_pattern,
        file_ext,
        run_dirs['output'],
        n_split_max,
        warning,
        w_log
    )

    inst.process()

    return None, None
