# -*- coding: utf-8 -*-

"""FIND_EXPOSURES RUNNER

This module runs find_exposures: Identify exposures that are used in selected tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.find_exposures_package import find_exposures_script


@module_runner(version='1.0', file_pattern=['image'], file_ext='.fits',
               depends=['numpy', 'astropy', 'sip_tpv'], numbering_scheme='_0')
def find_exposures_runner(input_file_list, output_dir, file_number_string,
                          config, w_log):

    input_file_name = input_file_list[0]

    inst = find_exposures_script.find_exposures(input_file_name, output_dir,
                                                file_number_string, config,
                                                w_log)
    inst.process()

    return None, None
