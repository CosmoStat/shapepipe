# -*- coding: utf-8 -*-

"""MCCD INTERPOLATION.

This file is the pipeline runner for the MCCD_interpolation package.

:Author: Tobias Liaudat based on Axel Guinot's code

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.MCCD_package import mccd_interpolation_script\
    as mccd_interp
from shapepipe.modules.MCCD_package import shapepipe_auxiliary_mccd as aux_mccd
import os


@module_runner(input_module=['setools_runner'],
               version='1.0',
               file_pattern=['galaxy_selection'],
               file_ext=['.npy'],
               depends=['numpy', 'astropy', 'galsim', 'sqlitedict', 'mccd'],
               run_method='parallel')
def mccd_interp_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):

    mode = config.getexpanded('MCCD_INTERP_RUNNER', 'MODE')
    pos_params = config.getlist('MCCD_INTERP_RUNNER', 'POSITION_PARAMS')
    get_shapes = config.getboolean('MCCD_INTERP_RUNNER', 'GET_SHAPES')
    mccd_model_extension = '.npy'
    output_dir = run_dirs['output']

    if mode == 'CLASSIC':
        psf_model_dir = config.getexpanded('MCCD_INTERP_RUNNER',
                                           'PSF_MODEL_DIR')
        psf_model_pattern = config.get('MCCD_INTERP_RUNNER',
                                       'PSF_MODEL_PATTERN')
        psf_separator = config.get('MCCD_INTERP_RUNNER',
                                   'PSF_MODEL_SEPARATOR')
        # psfcat_path, galcat_path = input_file_list
        galcat_path = input_file_list[0]

        # verify that the MCCD model exists
        ccd_id = file_number_string.split(psf_separator)[-1]
        exposure_id = file_number_string.split(psf_separator)[-2]

        psf_model_path = psf_model_dir + psf_model_pattern + psf_separator +\
            exposure_id + mccd_model_extension

        saving_path = output_dir + '/galaxy_psf' + file_number_string +\
            '.fits'

        if not os.path.exists(psf_model_path):
            error_msg = "The corresponding PSF model was not not found."
            w_log.info("Error message: {}.\n On catalog with id: {}.".format(
                error_msg, file_number_string))

            return None, None

        w_log.info('Interpolating catalog %s..' % file_number_string)
        output_msg = aux_mccd.mccd_interpolation_pipeline(
            mccd_model_path=psf_model_path,
            galcat_path=galcat_path,
            pos_params=pos_params,
            ccd_id=int(ccd_id),
            saving_path=saving_path,
            get_shapes=get_shapes)

        if output_msg is not None:
            w_log.info("Error message: {}.\n On catalog with id: {}.".format(
                output_msg, file_number_string))

    elif mode == 'MULTI-EPOCH':
        psf_model_dir = config.getexpanded('MCCD_INTERP_RUNNER',
                                           'PSF_MODEL_DIR')
        psf_model_pattern = config.get('MCCD_INTERP_RUNNER',
                                       'PSF_MODEL_PATTERN')
        f_wcs_path = config.getexpanded('MCCD_INTERP_RUNNER',
                                        'ME_LOG_WCS')

        galcat_path = input_file_list[0]

        inst = mccd_interp.MCCDinterpolator(None,
                                            galcat_path,
                                            output_dir,
                                            file_number_string,
                                            w_log,
                                            pos_params,
                                            get_shapes)

        inst.process_me(psf_model_dir, psf_model_pattern, f_wcs_path)

    elif mode == 'VALIDATION':
        ValueError('''MODE has to be in MULTI-EPOCH or CLASSIC.
        For validation use MCCD validation runner''')

    else:
        ValueError('MODE has to be in : [CLASSIC, MULTI-EPOCH]')

    return None, None
