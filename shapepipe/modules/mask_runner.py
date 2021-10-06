# -*- coding: utf-8 -*-

"""MASK RUNNER

This file is the pipeline runner for the mask package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.mask_package.mask_script import mask


@module_runner(version='1.0',
               file_pattern=['image', 'weight', 'flag'],
               file_ext=['.fits', '.fits', '.fits'],
               depends=['numpy', 'astropy'], executes=['ww', 'findgsc2.2'],
               numbering_scheme='_0')
def mask_runner(input_file_list, run_dirs, file_number_string,
                config, module_config_sec, w_log):

    if len(input_file_list) == 2:
        ext_flag_name = None
        ext_star_cat = None
    elif len(input_file_list) == 3:
        if config.getboolean(module_config_sec, 'USE_EXT_FLAG'):
            ext_flag_name = input_file_list[2]
            ext_star_cat = None
        elif config.getboolean(module_config_sec, 'USE_EXT_STAR'):
            ext_flag_name = None
            ext_star_cat = input_file_list[2]
        else:
            raise ValueError('Expecting external flag or external star '
                             'catalog.')
    elif len(input_file_list) == 4:
        if (config.getboolean(module_config_sec, 'USE_EXT_FLAG') and
                config.getboolean(module_config_sec, 'USE_EXT_STAR')):
            ext_flag_name = input_file_list[2]
            ext_star_cat = input_file_list[3]
        else:
            raise ValueError('Expecting external flag and external star '
                             'catalog.')
    else:
        raise ValueError("Input file list of length {} found, must be "
                         "'image', 'weight' and 'ext_flags', 'ext_star_cat' "
                         "(optional)"
                         "".format(len(input_file_list)))

    config_file = config.getexpanded(module_config_sec, 'MASK_CONFIG_PATH')

    if config.has_option(module_config_sec, 'HDU'):
        hdu = config.getint(module_config_sec, 'HDU')
    else:
        hdu = 0

    if config.has_option(module_config_sec, 'SUFFIX'):
        suffix = config.get(module_config_sec, 'SUFFIX')
    else:
        suffix = ''

    if config.has_option(module_config_sec, 'OUTNAME_BASE'):
        outname_base = config.get(module_config_sec, 'OUTNAME_BASE')
    else:
        outname_base = 'flag'

    inst = mask(*input_file_list[:2], suffix.replace(" ", ""),
                file_number_string, config_file, run_dirs['output'],
                path_external_flag=ext_flag_name, outname_base=outname_base,
                star_cat_path=ext_star_cat, hdu=hdu)
    stdout, stderr = inst.make_mask()

    return stdout, stderr
