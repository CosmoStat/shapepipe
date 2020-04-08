# -*- coding: utf-8 -*-

"""ERASE OUTPUT RUNNER

This module erase all the output of a specific module.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline.execute import execute


<<<<<<< HEAD
def get_module_run_dir(output_dir):
    """
    """

    return '/' + '/'.join(re.split('/', output_dir)[1:-2]) + '/psfex_runner/output'


@module_runner(input_module='sextractor_runner',
               version='1.0',
               file_pattern=['output_pattern'],
               file_ext=['.output_ext'],
               numbering_scheme='_0')
def erase_output_runner(input_file_list, output_dir, file_number_string,
                        config, w_log):

    # module_name = config.get("ERASE_OUTPUT_RUNNER", "MODULE_NAME")



=======
@module_runner(version='1.0',
               file_pattern=['output_pattern'],
               file_ext=['.output_ext'],
               numbering_scheme='_0')
def erase_output_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):

>>>>>>> 6f38177730a70145384f44c8990da7fc70ee1f9d
    stdout, stderr = execute('rm -f {0}'.format(*input_file_list))

    return stdout, stderr
