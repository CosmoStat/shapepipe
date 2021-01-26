# -*- coding: utf-8 -*-

"""PSFEX_INTERP RUNNER

This file is the pipeline runner for the PSFExInterpolation package.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.PSFExInterpolation_package import interpolation_script


@module_runner(input_module=['psfex_runner', 'setools_runner'], version='1.0',
               file_pattern=['star_selection', 'galaxy_selection'],
               file_ext=['.psf', '.fits'],
               depends=['numpy', 'astropy', 'galsim', 'sqlitedict'])
def psfex_interp_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):

    mode = config.get('PSFEX_INTERP_RUNNER', 'MODE')

    pos_params = config.getlist('PSFEX_INTERP_RUNNER', 'POSITION_PARAMS')
    get_shapes = config.getboolean('PSFEX_INTERP_RUNNER', 'GET_SHAPES')
    star_thresh = config.getint('PSFEX_INTERP_RUNNER', 'STAR_THRESH')
    chi2_thresh = config.getint('PSFEX_INTERP_RUNNER', 'CHI2_THRESH')

    if mode == 'CLASSIC':
        psfcat_path, galcat_path = input_file_list

        inst = interpolation_script.PSFExInterpolator(psfcat_path, galcat_path,
                                                      run_dirs['output'],
                                                      file_number_string,
                                                      w_log, pos_params,
                                                      get_shapes, star_thresh,
                                                      chi2_thresh)
        inst.process()

    elif mode == 'MULTI-EPOCH':
        dot_psf_dir = config.getexpanded('PSFEX_INTERP_RUNNER',
                                         'ME_DOT_PSF_DIR')
        dot_psf_pattern = config.get('PSFEX_INTERP_RUNNER',
                                     'ME_DOT_PSF_PATTERN')
        f_wcs_path = config.getexpanded('PSFEX_INTERP_RUNNER', 'ME_LOG_WCS')

        galcat_path = input_file_list[0]

        inst = interpolation_script.PSFExInterpolator(None, galcat_path,
                                                      run_dirs['output'],
                                                      file_number_string,
                                                      w_log, pos_params,
                                                      get_shapes, star_thresh,
                                                      chi2_thresh)

        inst.process_me(dot_psf_dir, dot_psf_pattern, f_wcs_path)

    elif mode == 'VALIDATION':
        psfcat_path, galcat_path, psfex_cat_path = input_file_list

        inst = interpolation_script.PSFExInterpolator(psfcat_path, galcat_path,
                                                      run_dirs['output'],
                                                      file_number_string,
                                                      w_log, pos_params,
                                                      get_shapes, star_thresh,
                                                      chi2_thresh)

        inst.process_validation(psfex_cat_path)

    else:
        ValueError('MODE has to be in : [C, ME]')

    return None, None
