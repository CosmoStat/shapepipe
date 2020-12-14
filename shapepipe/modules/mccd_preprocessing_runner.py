# -*- coding: utf-8 -*-

"""MCCD PREPROCESSING RUNNER.

This file is the pipeline runner for the preprocessing of the inputs for the
MCCD algorithm.

:Author: Tobias Liaudat

"""
from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.MCCD_package import shapepipe_auxiliary_mccd as aux_mccd
import mccd


@module_runner(input_module=['setools_runner'], version='1.0',
               file_pattern=['star_split_ratio_80', 'star_split_ratio_20'],
               file_ext=['.fits', '.fits'],
               depends=['numpy', 'mccd', 'galsim', 'astropy'],
               run_method='serial')
def mccd_preprocessing_runner(input_file_list, run_dirs, file_number_string,
                              config, w_log):
    # Recover the MCCD config file and its params
    config_file_path = config.getexpanded('MCCD', 'CONFIG_PATH')
    mccd_mode = config.get('MCCD', 'MODE')
    verbose = config.getboolean('MCCD', 'VERBOSE')

    # Parse MCCD config file
    mccd_parser = mccd.auxiliary_fun.MCCDParamsParser(config_file_path)
    mccd_parser.parse_document()
    input_dict = mccd_parser.get_inputs_kw()

    # Extract useful parameters
    separator = input_dict['separator']
    min_n_stars = input_dict['min_n_stars']
    outlier_std_max = input_dict['outlier_std_max']

    if mccd_mode == 'FIT':
        input_file_pos_list = [0]
        save_name_list = ['train_star_selection']
        min_n_stars_list = [min_n_stars]

    elif mccd_mode == 'FIT_VALIDATION':
        if len(input_file_list[0]) == 1:
            input_file_pos_list = [0]
            save_name_list = ['train_star_selection']
            min_n_stars_list = [min_n_stars]
        else:
            input_file_pos_list = [0, 1]
            save_name_list = ['train_star_selection', 'test_star_selection']
            min_n_stars_list = [min_n_stars, 1]

    elif mccd_mode == 'VALIDATION':
        input_file_pos_list = [0]
        save_name_list = ['test_star_selection']
        min_n_stars_list = [1]

    else:
        raise ValueError('''MODE should be in ["FIT", "FIT_VALIDATION",
         "VALIDATION"].''')

    # Use the outfile from the pipeline and ignore the output directory from
    # the MCCD config file
    # Output paths for both newly generates datasets
    output_mccd_path = run_dirs['output'] + '/'

    [aux_mccd.mccd_preprocessing_pipeline(
     input_file_list=input_file_list,
     output_path=output_mccd_path,
     input_file_position=_input_pos,
     min_n_stars=_min_stars,
     separator=separator,
     CCD_id_filter_list=None,
     outlier_std_max=outlier_std_max,
     save_masks=False,
     save_name=_save_name,
     save_extension='.fits',
     verbose=verbose,
     print_fun=w_log.info)
     for _input_pos, _save_name, _min_stars in
     zip(input_file_pos_list, save_name_list, min_n_stars_list)]

    return None, None
