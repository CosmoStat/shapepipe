# -*- coding: utf-8 -*-

"""MASK RUNNER

Module runner for ``mask``.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.mask_package.mask import Mask


@module_runner(
    version='1.0',
    file_pattern=['image', 'weight', 'flag'],
    file_ext=['.fits', '.fits', '.fits'],
    depends=['numpy', 'astropy'],
    executes=['ww', 'findgsc2.2'],
    numbering_scheme='_0',
)
def mask_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):

    # Get number of input files
    n_inputs = len(input_file_list)

    # Set options for 2 inputs
    if n_inputs == 2:
        ext_flag_name = None
        ext_star_cat = None

    # Set options for 3 inputs
    elif n_inputs == 3:
        if config.getboolean('MASK_RUNNER', 'USE_EXT_FLAG'):
            ext_flag_name = input_file_list[2]
            ext_star_cat = None
        elif config.getboolean('MASK_RUNNER', 'USE_EXT_STAR'):
            ext_flag_name = None
            ext_star_cat = input_file_list[2]
        else:
            raise ValueError(
                f'Found {n_inputs} inputs but was expecting external flag or '
                + 'external star catalogue in the MASK_RUNNER section of the '
                + 'config file.'
            )

    # Set options for 4 inputs
    elif n_inputs == 4:
        if (
            config.getboolean('MASK_RUNNER', 'USE_EXT_FLAG')
            and config.getboolean('MASK_RUNNER', 'USE_EXT_STAR')
        ):
            ext_flag_name = input_file_list[2]
            ext_star_cat = input_file_list[3]
        else:
            raise ValueError(
                f'Found {n_inputs} inputs but was expecting external flag and '
                + 'external star catalogue in the MASK_RUNNER section of the '
                + 'config file.'
            )

    # Raise error for invalid settings
    else:
        raise ValueError(
            f'Found {n_inputs} inputs and these must be "image", "weight" and '
            + '"ext_flags", "ext_star_cat" (optional). Check the MASK_RUNNER '
            + 'section of the config file to make sure you have the '
            + 'appropriate settings.'
        )

    # Get path to mask configuration options
    config_file = config.getexpanded('MASK_RUNNER', 'MASK_CONFIG_PATH')

    # Get mask HDU number
    if config.has_option('MASK_RUNNER', 'HDU'):
        hdu = config.getint('MASK_RUNNER', 'HDU')
    else:
        hdu = 0

    # Get mask mask file name suffix
    if config.has_option('MASK_RUNNER', 'SUFFIX'):
        suffix = config.get('MASK_RUNNER', 'SUFFIX')
    else:
        suffix = ''

    # Get mask base ouput file name
    if config.has_option('MASK_RUNNER', 'OUTNAME_BASE'):
        outname_base = config.get('MASK_RUNNER', 'OUTNAME_BASE')
    else:
        outname_base = 'flag'

    # Create instance of Mask
    inst = Mask(
        *input_file_list[:2],
        image_suffix=suffix.replace(" ", ""),
        image_num=file_number_string,
        config_filepath=config_file,
        output_dir=run_dirs['output'],
        path_external_flag=ext_flag_name,
        outname_base=outname_base,
        star_cat_path=ext_star_cat,
        hdu=hdu,
    )

    # Make the mask
    stdout, stderr = inst.make_mask()

    # Return stdout and stderr
    return stdout, stderr
