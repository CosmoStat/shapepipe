# -*- coding: utf-8 -*-

"""MCCD FIT VAL RUNNER.

This file is the pipeline fit and validation runner for the MCCD package.

:Author: Tobias Liaudat

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.MCCD_package import shapepipe_auxiliary_mccd as aux_mccd
import mccd


@module_runner(input_module=['mccd_preprocessing_runner'], version='1.0',
               file_pattern=['train_star_selection', 'test_star_selection'],
               file_ext=['.fits', '.fits'], numbering_scheme='-0000000',
               depends=['numpy', 'mccd', 'galsim'], run_method='parallel')
def mccd_fit_val_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):
    # Recover the MCCD config file and its params
    config_file_path = config.getexpanded('MCCD', 'CONFIG_PATH')
    mccd_mode = config.get('MCCD', 'MODE')
    verbose = config.getboolean('MCCD', 'VERBOSE')

    # Parse MCCD config file
    mccd_parser = mccd.auxiliary_fun.MCCDParamsParser(config_file_path)
    mccd_parser.parse_document()

    # Prepare inputs to run the main fit function
    output_dir = run_dirs['output'] + '/'
    fit_saving_name = 'fitted_model'
    val_saving_name = 'validation_psf'
    # Extract the training star catalog
    trainstar_path = input_file_list[0]
    # Validation stars are in the second position of the list
    teststar_path = input_file_list[1]

    if mccd_mode == 'FIT_VALIDATION':

        aux_mccd.mccd_fit_pipeline(
            trainstar_path=trainstar_path,
            file_number_string=file_number_string,
            mccd_parser=mccd_parser,
            output_dir=output_dir,
            verbose=verbose,
            saving_name=fit_saving_name)

        # Fitted model is found in the output directory
        mccd_model_path = output_dir + fit_saving_name + file_number_string \
            + '.npy'

        aux_mccd.mccd_validation_pipeline(
            teststar_path=teststar_path,
            mccd_model_path=mccd_model_path,
            mccd_parser=mccd_parser,
            output_dir=output_dir,
            file_number_string=file_number_string,
            w_log=w_log,
            val_saving_name=val_saving_name)

    else:
        raise ValueError('''mccd_fit_val_runner should be called when the
        MODE is "FIT_VALIDATION".''')

    return None, None
