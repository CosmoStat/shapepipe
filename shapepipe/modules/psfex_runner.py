# -*- coding: utf-8 -*-

"""PSFEX RUNNER

This module run PSFEx.

:Author: Axel Guinot

"""

import re
import os
from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.psfex_package.psfex_script import PSFEx_caller


@module_runner(
    input_module='setools_runner',
    version='1.0',
    file_pattern=['star_selection'],
    file_ext=['.fits'],
    executes='psfex'
)
def psfex_runner(
        input_file_list,
        run_dirs,
        file_number_string,
        config,
        module_config_sec,
        w_log
):

    # extract psfex  run configurations
    psfex_executable_path = config.getexpanded(
        module_config_sec,
        "EXEC_PATH"
    )
    output_dir = run_dirs['output']

    outcatalog_name = f"{output_dir}/psfex_cat{file_number_string}.cat"

    psfex_config_file = config.getexpanded(
        module_config_sec,
        "DOT_PSFEX_FILE"
    )

    input_file_path = input_file_list[0]

    # check image options
    if config.has_option(module_config_sec, "CHECKIMAGE"):
        check_image_list = config.getlist(module_config_sec, "CHECKIMAGE")
    else:
        check_image_list = ['']

    # prepare the psfex command line
    PSFEx_call = PSFEx_caller(
        psfex_executable_path,
        input_file_path,
        psfex_config_file,
        output_dir,
        outcatalog_name,
        check_image_list
    )

    # generates the psfex command
    command_line = PSFEx_call.generate_command()

    w_log.info(f'Running command \'{command_line}\'')
    stderr, stdout = execute(command_line)

    # move psfex errors reported as stdout to stderr
    stdout, stderr = PSFEx_call.parse_errors(stderr, stdout)

    return stdout, stderr
