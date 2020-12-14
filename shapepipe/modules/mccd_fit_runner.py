# -*- coding: utf-8 -*-

"""MCCD FIT RUNNER.

This file is the pipeline fit runner for the MCCD package.

:Author: Tobias Liaudat

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.MCCD_package import shapepipe_auxiliary_mccd as aux_mccd
import mccd


@module_runner(input_module=['mccd_preprocessing_runner'], version='1.0',
               file_pattern=['train_star_selection'],
               file_ext=['.fits'], numbering_scheme='-0000000',
               depends=['numpy', 'mccd', 'galsim'], run_method='parallel')
def mccd_fit_runner(input_file_list, run_dirs, file_number_string,
                    config, w_log):
    # Recover the MCCD config file and its params
    config_file_path = config.getexpanded('MCCD', 'CONFIG_PATH')
    mccd_mode = config.get('MCCD', 'MODE')
    verbose = config.getboolean('MCCD', 'VERBOSE')

    # Parse MCCD config file
    mccd_parser = mccd.auxiliary_fun.MCCDParamsParser(config_file_path)
    mccd_parser.parse_document()

    # Prepare inputs to run the main fit function
    trainstar_path = input_file_list[0]
    output_dir = run_dirs['output'] + '/'
    saving_name = 'fitted_model'

    if mccd_mode == 'FIT':
        aux_mccd.mccd_fit_pipeline(
            trainstar_path, file_number_string, mccd_parser, output_dir,
            verbose, saving_name)

    else:
        raise ValueError('''mccd_fit_runner should be called when the
        MODE is "FIT".''')

    return None, None
