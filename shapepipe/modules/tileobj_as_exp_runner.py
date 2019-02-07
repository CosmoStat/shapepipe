# -*- coding: utf-8 -*-

"""TILEOBJ_AS_EXP_RUNNER

This module runs tileobj_as_exp_runner: Write objects selected on tiles to catalogues in exposure-single-CCD format.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: January 2019

"""


from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner

from shapepipe.modules.tileobj_as_exp_package import tileobj_as_exp_script as tileobj_as_exp


@module_runner(version='1.0', file_pattern=['image', 'cat'], file_ext=['.fits', '.fits'],
               depends=['numpy', 'astropy'], numbering_scheme='_0')
def tileobj_as_exp_runner(input_file_list, output_dir, file_number_string, config, w_log):

    input_file_names = input_file_list

    inst = tileobj_as_exp.tileobj_as_exp(input_file_names, output_dir, file_number_string, config, w_log)
    inst.process()

    return None, None

