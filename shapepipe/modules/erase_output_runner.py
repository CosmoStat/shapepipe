# -*- coding: utf-8 -*-

"""ERASE OUTPUT RUNNER

This module erase all the output of a specific module.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline.execute import execute


@module_runner(version='1.0',
               file_pattern=['output_pattern'],
               file_ext=['.output_ext'],
               numbering_scheme='_0')
def erase_output_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):

    stdout, stderr = execute('rm -f {0}'.format(*input_file_list))

    return stdout, stderr
