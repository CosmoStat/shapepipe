# -*- coding: utf-8 -*-

"""SETOOLS RUNNER

This file is the pipeline runner for the SETools package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.SETools_package import SETools_script as setools


@module_runner(input_module='sextractor_runner', version='1.0',
               file_pattern=['sexcat'], file_ext=['.fits'],
               depends=['numpy', 'matplotlib'])
def setools_runner(input_file_list, output_dir, file_number_string,
                   config, w_log):

    config_file = config.get('SETOOLS_RUNNER', 'SETOOLS_CONFIG_PATH')

    inst = setools.SETools(input_file_list[0], config_file, output_dir)
    inst.process()

    return None, None
