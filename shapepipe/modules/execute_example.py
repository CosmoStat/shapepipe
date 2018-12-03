# -*- coding: utf-8 -*-

"""EXECUTE MODULE EXAMPLE

This module defines methods for an example command line execution module.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module='python_example', file_pattern='process',
               file_ext='.cat')
def execute_example(worker_dict, filehd, config, w_log):

    command_line = 'head {}'.format(worker_dict['process'])
    output_file_name = (filehd.format(filehd.output_dir,
                        '{}_head_out.txt'.format(worker_dict['job_name'])))

    stdout, stderr = execute(command_line)

    text_file = open(output_file_name, 'w')
    text_file.write(stdout)

    return stdout, stderr
