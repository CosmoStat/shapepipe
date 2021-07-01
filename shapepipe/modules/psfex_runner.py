# -*- coding: utf-8 -*-

"""PSFEX RUNNER

This module run PSFEx.

:Author: Axel Guinot

"""

import re
import os
from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.psfex_package.psfex import PSFex_caller


@module_runner(input_module='setools_runner',
        version='1.0',
        file_pattern=['star_selection'],
        file_ext=['.fits'],
        executes='psfex'
)

def psfex_runner(input_file_list,
        run_dirs,
        file_number_string,
        config,
        w_log
    ):
    """
    Runs the psfex wrapper package.

    Parameters:
    ----------
    input_file_list: list
        List with the full path to the input file
    run_dirs: dict
        Dictionary with directories
    file_number_string: str
        Specific number of the data file being processed
    config:
        Config object
    w_log:
        Log object


    Returns:

    stdout, stderr: str
        Strings with the output and error output of execution.

    """
    inputs = [input_file_list,
                     run_dirs,
                     file_number_string,
                     config,
                     w_log]

    # extract psfex  run configurations
    psfex_executable_path = config.getexpanded("PSFEX_RUNNER",
                                "EXEC_PATH")
    output_dir = run_dirs['output']

    outcatalog_name = f'psfex_cat{file_number_string}.cat'
    psfex_config_file = config.getexpanded("PSFEX_RUNNER",
                            "DOT_PSFEX_FILE")
    input_file_path = input_file_list[0]

    # check image options
    if config.has_option('PSFEX_RUNNER', "CHECKIMAGE"):
        check_image_list = config.getlist("PSFEX_RUNNER", "CHECKIMAGE")
    else:
        check_image_list = ['']

    # prepare the psfex command line
    PSFex_call  = PSFex_caller(psfex_executable_path,
            input_file_path,
            psfex_config_file,
            output_dir,
            outcatalog_name,
            check_image_list)

    # generates the psfex command
    command_line = PSFex_call.generate_command()

    w_log.info(f'Running command \'{command_line}\'')
    stderr, stdout = execute(command_line)

    # move psfex errors reported as stdout to stderr
    check_error = re.findall('error', stdout.lower())
    check_error2 = re.findall('all done', stdout.lower())

    if check_error == []:
        stderr2 = ''
    else:
        stderr2 = stdout

    if check_error2 == []:
        stderr2 = stdout

    return stdout, stderr2
