# -*- coding: utf-8 -*-

"""PSFEXINTERP RUNNER

This file is the pipeline runner for the PSFExInterpolation package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.PSFExInterpolation_package import interpolation_script


@module_runner(input_module=['psfex_runner', 'setools_runner'], version='1.0',
               file_pattern=['star_selection', 'galaxy_selection'],
               file_ext=['.psf', '.fits'],
               depends=['numpy', 'astropy', 'galsim'])
def psfexinterp_runner(input_file_list, output_dir, file_number_string,
                       config, w_log):

    mode = config.get('PSFEXINTERP_RUNNER', 'MODE')

    pos_params = config.getlist('PSFEXINTERP_RUNNER', 'POSITION_PARAMS')
    get_shapes = config.getboolean('PSFEXINTERP_RUNNER', 'GET_SHAPES')
    star_thresh = config.getint('PSFEXINTERP_RUNNER', 'STAR_THRESH')

    if mode == 'C':
        psfcat_path, galcat_path = input_file_list

        inst = interpolation_script.PSFExInterpolator(psfcat_path, galcat_path,
                                                      output_dir, file_number_string,
                                                      pos_params, get_shapes, star_thresh)
        inst.process()

    elif mode == 'ME':
        dot_psf_dir = config.getexpanded('PSFEXINTERP_RUNNER', 'ME_DOT_PSF_DIR')
        dot_psf_pattern = config.get('PSFEXINTERP_RUNNER', 'ME_DOT_PSF_PATTERN')
        f_wcs_path = config.getexpanded('PSFEXINTERP_RUNNER', 'ME_LOG_WCS')

        galcat_path = input_file_list[0]

        inst = interpolation_script.PSFExInterpolator(None, galcat_path,
                                                      output_dir, file_number_string,
                                                      pos_params, get_shapes, star_thresh)

        inst.process_me(dot_psf_dir, dot_psf_pattern, f_wcs_path)

    else:
        ValueError('MODE has to be in : [C, ME]')

    return None, None
