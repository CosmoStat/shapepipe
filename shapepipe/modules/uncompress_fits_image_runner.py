# -*- coding: utf-8 -*-

"""UNCOMPRESS FITS IMAGE RUNNER

Module runner for ``uncompress_fits_image_runner``
This module uncompress fits images and save them on a single hdu fits.

:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2020

:Package: ShapePipe

"""

from shapepipe.modules.module_decorator import module_runner
import shapepipe.modules.uncompress_fits_image_package.uncompress as uz


@module_runner(
    version='1.1',
    file_pattern=['image'],
    file_ext=['.fits'],
    numbering_scheme='_0'
)
def uncompress_fits_image_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    # Get output patterns
    output_pattern_list = config.getlist(module_config_sec, 'OUTPUT_PATTERN')

    # Get HDU number
    if config.has_option(module_config_sec, 'HDU_DATA'):
        data_hdu = config.getint(module_config_sec, 'HDU_DATA')
    else:
        data_hdu = 0

    # Check consistency of input and output list lengths
    if len(input_file_list) != len(output_pattern_list):
        raise ValueError(
            f'Lists INPUT_PATH ({len(input_file_list)})'
            + f' and OUTPUT_PATTERN ({len(output_pattern_list)})'
            + 'need to be of equal length.'
        )

    # Create instance of uncompress
    uz_inst = uz.Uncompress(
        input_file_list,
        output_pattern_list,
        run_dirs['output'],
        file_number_string,
        data_hdu)

    uz_inst.process()

    return None, None
