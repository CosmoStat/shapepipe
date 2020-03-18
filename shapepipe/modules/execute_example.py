# -*- coding: utf-8 -*-

"""EXECUTE MODULE EXAMPLE

This module defines methods for an example command line execution module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='python_example', version='1.0',
               file_pattern='pyex_output', file_ext='.cat', executes='head',
               run_method='parallel')
def execute_example(input_file_list, run_dirs, file_number_string, *args):

    command_line = 'head {}'.format(input_file_list[0])
    output_file_name = '{}/head_output{}.txt'.format(run_dirs['output'],
                                                     file_number_string)

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, 'w')
    text_file.write(stdout)

    return stdout, stderr
