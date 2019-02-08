# -*- coding: utf-8 -*-

"""COMBINE_MEXP_RUNNER

This module runs combine_mexp_runner: Combine information on objects and PSF
from multiple exposure-single-CCD catalogues.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: February 2019

"""


from shapepipe.pipeline.execute import execute
from shapepipe.modules.module_decorator import module_runner

from shapepipe.modules.combine_mexp_package import combine_mexp_script as combine_mexp


@module_runner(version='1.0', file_pattern=['image', 'cat.exp.img', 'galaxy_psf'],
               file_ext=['.fits', '.fits', '.fits'],
               depends=['numpy', 'astropy'], numbering_scheme='_0')
def combine_mexp_runner(input_file_list, output_dir, file_number_string, config, w_log):

    input_file_names = input_file_list

    inst = combine_mexp.combine_mexp(input_file_names, output_dir, file_number_string, config, w_log)
    inst.process()

    return None, None
