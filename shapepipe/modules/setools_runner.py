# -*- coding: utf-8 -*-

"""SETOOLS RUNNER

This file is the pipeline runner for the SETools package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.SETools_package import SETools_script as setools


@module_runner(input_module='sextractor_runner_exp', version='1.0',
               file_pattern=['sexcat'], file_ext=['.fits'],
               depends=['numpy', 'matplotlib'])
def setools_runner(input_file_list, run_dirs, file_number_string,
                   config, w_log):

    config_file = config.getexpanded('SETOOLS_RUNNER', 'SETOOLS_CONFIG_PATH')

    inst = setools.SETools(input_file_list[0], run_dirs['output'],
                           file_number_string, config_file)
    inst.process(w_log)

    return None, None
