# -*- coding: utf-8 -*-

""" NGMIX RUNNER

This file contains methods to run ngmix for shape measurement.

:Author: Axel Guinot

"""

import re
import numpy as np

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.ngmix_package import ngmix


@module_runner(
    input_module=[
        'sextractor_runner',
        'psfex_interp_runner',
        'vignetmaker_runner'
    ],
    version='0.0.1',
    file_pattern=[
        'tile_sexcat',
        'image',
        'exp_background',
        'galaxy_psf',
        'weight',
        'flag'
    ],
    file_ext=['.fits', '.sqlite', '.sqlite', '.sqlite', '.sqlite', '.sqlite'],
    depends=['numpy', 'ngmix', 'galsim', 'sqlitedict', 'astropy']
)
def ngmix_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    # Read config file entries

    # Photometric zero point
    zero_point = config.getfloat(module_config_sec, 'MAG_ZP')

    # Pixel scale
    pixel_scale = config.getfloat(module_config_sec, 'PIXEL_SCALE')

    # Path to merged single-exposure single-HDU headers
    f_wcs_path = config.getexpanded(module_config_sec, 'LOG_WCS')

    # First and last galaxy ID to process
    id_obj_min = config.getint(module_config_sec, 'ID_OBJ_MIN')
    id_obj_max = config.getint(module_config_sec, 'ID_OBJ_MAX')

    # Initialise class instance
    inst = ngmix.Ngmix(
        input_file_list,
        run_dirs['output'],
        file_number_string,
        zero_point,
        pixel_scale,
        f_wcs_path,
        w_log,
        id_obj_min=id_obj_min,
        id_obj_max=id_obj_max
    )

    # Process ngmix shape measurement and metacalibration
    inst.process()

    return None, None
