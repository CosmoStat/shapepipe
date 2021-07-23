# -*- coding: utf-8 -*-

"""SETOOLS RUNNER

Module runner for ``setools``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.setools_package import setools


@module_runner(
    input_module='sextractor_runner_exp',
    version='1.0',
    file_pattern=['sexcat'],
    file_ext=['.fits'],
    depends=['numpy', 'matplotlib']
)
def setools_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    w_log,
):

    # Get path to setools configuration file
    config_file = config.getexpanded('SETOOLS_RUNNER', 'SETOOLS_CONFIG_PATH')

    # Create instance of SETools
    se_inst = setools.SETools(
        input_file_list[0],
        run_dirs['output'],
        file_number_string,
        config_file,
    )

    # Process inputs
    se_inst.process(w_log)

    # No return objects
    return None, None
