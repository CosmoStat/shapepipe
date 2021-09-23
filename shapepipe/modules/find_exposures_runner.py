# -*- coding: utf-8 -*-

"""FIND_EXPOSURES RUNNER

Module runner for ``find_exposures``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.find_exposures_package import find_exposures_script as fe


@module_runner(
    version='1.1',
    file_pattern=['image'],
    file_ext='.fits',
    depends=['numpy', 'astropy', 'sip_tpv'],
    numbering_scheme='_0'
)
def find_exposures_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    # Get input file name of tile
    input_file_name = input_file_list[0]

    # Create output ascii file name
    output_path = f'{run_dirs["output"]}/exp_numbers{file_number_string}.txt'

    # Create class
    inst = fe.FindExposures(
        input_file_name,
        output_path,
        w_log
    )

    # Run processing
    inst.process()

    # No return objects
    return None, None
