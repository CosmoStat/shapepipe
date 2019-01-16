# -*- coding: utf-8 -*-

"""EXECUTE MODULE EXAMPLE

This module defines methods for an example command line execution module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='python_example', version='1.0',
               file_pattern='process', file_ext='.cat', executes='head')
def execute_example(input_file_list, output_dir, job_name, *args):

    command_line = 'head {}'.format(input_file_list[0])
    output_file_name = '{}/{}_head_out.txt'.format(output_dir, job_name)

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, 'w')
    text_file.write(stdout)

    return stdout, stderr
