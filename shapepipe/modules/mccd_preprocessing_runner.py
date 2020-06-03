# -*- coding: utf-8 -*-

"""MCCD PREPROCESSING RUNNER

This file is the pipeline runner for the preprocessing of the inputs for the
MCCD algorithm.

:Author: Tobias Liaudat

"""
from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
from astropy.io import fits
import numpy as np
import mccd_rca.my_mccd_rca as rca
import mccd_rca.utils as utils
import mccd_rca.mccd_utils as mccd_utils
try:
    import galsim.hsm as hsm
    from galsim import Image
    import_fail = False
except ImportError:
    import_fail = True


@module_runner(input_module=['setools_runner'], version='1.0',
               file_pattern=['star_split_ratio_80','star_split_ratio_20'],
               file_ext=['.fits','.fits'],
               depends=['numpy', 'mccd_rca', 'galsim', 'astropy'],
               run_method='serial')
def mccd_preprocessing_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):

    save_masks = False
    numbering_scheme = config.get('FILE','NUMBERING_SCHEME')
    separator = numbering_scheme[0]

    # Output paths for both newly generates datasets
    output_mccd_path = run_dirs['output'] + '/'

    try:
        CCD_id_filter_list = config.get('MCCD_PREPROCESSING', 'CCDS')
        CCD_id_filter_list = np.array([int(elem) for elem in CCD_id_filter_list.split(',')])
    except Exception:
        CCD_id_filter_list = np.arange(40)
    min_train_n_stars = config.getint('MCCD_PREPROCESSING', 'MIN_TRAIN_STARS')
    min_val_n_stars = config.getint('MCCD_PREPROCESSING', 'MIN_VAL_STARS')

    # -------------------------- #
    # Preprocess
    mccd_star_nb = 0

    for it_tr_te in range(2):
        # Train and test catalog loop
        if it_tr_te == 0:
            print('Processing the training dataset..')
            mccd_inputs = mccd_utils.MccdInputs(separator=separator,save_masks=save_masks)
            # catalog_ids = mccd_inputs.preprocess_data(data_path, pipeline_train_pattern)
            catalog_ids = mccd_inputs.proprocess_pipeline_data(input_file_list, element_position=0)

            train_bool = True
            min_n_stars = min_train_n_stars

        elif it_tr_te == 1:
            print('Processing the validation dataset..')
            mccd_inputs = mccd_utils.MccdInputs(separator=separator,save_masks=save_masks)
            # catalog_ids = mccd_inputs.preprocess_data(data_path, pipeline_test_pattern)
            catalog_ids = mccd_inputs.proprocess_pipeline_data(input_file_list, element_position=1)
            train_bool = False
            min_n_stars = min_val_n_stars

        # Loop over the catalogs
        for it in range(catalog_ids.shape[0]):
            # For each observation position
            catalog_id = catalog_ids[it]
            star_list, pos_list, mask_list, ccd_list, SNR_list = mccd_inputs.get_inputs(catalog_id)
            print_fun = lambda x : w_log.info(x)
            star_list, pos_list, mask_list, ccd_list, SNR_list, _ = mccd_inputs.outlier_rejection(star_list, pos_list,
                                                                                                mask_list, ccd_list,
                                                                                                SNR_list, shape_std_max=5., print_fun=print_fun)

            mccd_star_list = []
            mccd_pos_list = []
            mccd_mask_list = []
            mccd_ccd_list = []
            mccd_SNR_list = []

            for j in range(len(star_list)):
                # For each CCD
                if ccd_list[j] in CCD_id_filter_list:
                    try:
                        # Separate the stars for train and test
                        n_stars = star_list[j].shape[2]

                        if n_stars >= min_n_stars:
                            mccd_star_list.append(star_list[j])
                            mccd_pos_list.append(pos_list[j])
                            mccd_mask_list.append(mask_list[j])
                            mccd_ccd_list.append(ccd_list[j] * np.ones(n_stars))
                            if SNR_list is not None:
                                mccd_SNR_list.append(SNR_list[j])
                        else:
                            msg = 'Not enough stars in catalog_id %s ,ccd %d. Total stars = %d' % (catalog_id, ccd_list[j],
                                                                                                   n_stars)
                            w_log.info(msg)
                            print(msg)

                    except Exception:
                        msg = 'Warning! Problem detected in catalog_id %s ,ccd %d' % (catalog_id, ccd_list[j])
                        print(msg)
                        w_log.info(msg)

            if mccd_pos_list:
                # If the list is not empty
                # Concatenate as fits can't handle list of numpy arrays and turn into reg format
                mccd_stars = utils.reg_format(np.concatenate(mccd_star_list, axis=2))
                mccd_poss = np.concatenate(mccd_pos_list, axis=0)
                mccd_ccds = np.concatenate(mccd_ccd_list, axis=0)

                if save_masks is True:
                    mccd_masks = utils.reg_format(np.concatenate(mccd_mask_list, axis=2))
                else:
                    # Send an array of False (None cannot be used in .fits)
                    mccd_masks = np.zeros((mccd_poss.shape[0]), dtype=bool)

                if SNR_list is not None:
                    mccd_SNRs = np.concatenate(mccd_SNR_list, axis=0)
                else:
                    # Send an array of False (None cannot be used in .fits)
                    mccd_SNRs = np.zeros((mccd_poss.shape[0]), dtype=bool)

                mccd_star_nb += mccd_stars.shape[0]

                # Save the fits file
                train_dic = {'VIGNET_LIST': mccd_stars, 'GLOB_POSITION_IMG_LIST': mccd_poss,
                             'MASK_LIST': mccd_masks, 'CCD_ID_LIST': mccd_ccds, 'SNR_WIN_LIST': mccd_SNRs}

                mccd_utils.save_fits(train_dic, train_bool, catalog_id, output_mccd_path,input_file_list[0][0])

        if it_tr_te == 0:
            print('Finished the training dataset processing.')
            print('Total training stars = %d' % (mccd_star_nb))
            mccd_star_nb = 0
        elif it_tr_te == 1:
            print('Finished the validation dataset processing.')
            print('Total validation stars = %d' % (mccd_star_nb))
            mccd_star_nb = 0


    return None, None
