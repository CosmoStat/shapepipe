# -*- coding: utf-8 -*-

"""MCCD auxiliary functions.

This module contains auxiliary functions for the MCCD package that are needed
by the MCCD runners.

:Author: Tobias Liaudat

"""

import os
import numpy as np
import mccd
from astropy.io import fits
import galsim
from shapepipe.pipeline import file_io as sc
import pprint

NOT_ENOUGH_STARS = 'Not enough stars to train the model.'


def mccd_preprocessing_pipeline(
    input_file_list,
    output_path,
    input_file_position=0,
    min_n_stars=20,
    separator='-',
    CCD_id_filter_list=None,
    outlier_std_max=100.,
    save_masks=True,
    save_name='train_star_selection',
    save_extension='.fits',
    verbose=True,
    print_fun=None
):
    r"""Preprocess input catalog.

    Parameters
    ----------
    input_file_list: list of list of str or list of str
        Input file list as taken from shapepipe's input.
    output_path: str
        Path to the folder where to save the preprocessed files.
    input_file_position: int
        Element position from the ``input_file_list`` to preprocess.
        Default is ``0``.
    min_n_stars: int
        Minimum number of stars in order to preprocess the CCD.
        Default is ``20``.
    separator: str
        Separator string that separates the catalog id and the CCD id.
        Default is ``'-'``.
    CCD_id_filter_list: list of int or None
        A list that correspond to the CCDs that should be preprocessed.
        If it is None, all the CCDs are preprocessed.
        (Current version: Hardcoded for the MegaCam scenario).
        Default is ``None``.
    outlier_std_max: float
        Parameter regulating the shape outlier removal. Default is very high
        so as it is not done at all. A decent number would be ``10``.
        Default is ``100.``.
    save_masks: bool
        If masks should be saved in the new file.
        Default is ``True``.
    save_name: str
        Name to save the preprocessed file.
        Default is ``'train_star_selection'``.
    save_extension: str
        Extension of the saved file.
        Default is ``.fits``.
    verbose: bool
        Verbose mode.
        Default is ``True``.
    print_fun: function or None


    Returns
    -------
    mccd_inputs: class
        An instance of ``MccdInputs`` class used for the input preprocessing.

    """
    mccd_star_nb = 0

    if CCD_id_filter_list is None:
        CCD_id_filter_list = np.arange(40)
    else:
        CCD_id_filter_list = np.array(CCD_id_filter_list)

    if verbose:
        if print_fun is None:
            def print_fun(x):
                print(x)
    else:
        def print_fun(x):
            pass

    print_fun('Processing dataset..')
    mccd_inputs = mccd.mccd_utils.MccdInputs(separator=separator)
    catalog_ids = mccd_inputs.proprocess_pipeline_data(
        input_file_list,
        element_position=input_file_position
    )

    # Loop over the catalogs
    for it in range(catalog_ids.shape[0]):
        # For each observation position
        catalog_id = catalog_ids[it]
        star_list, pos_list, mask_list, ccd_list, SNR_list, RA_list, \
            DEC_list = mccd_inputs.get_inputs(catalog_id)

        star_list, pos_list, mask_list, ccd_list, SNR_list, RA_list, \
            DEC_list, _ = mccd_inputs.outlier_rejection(
                star_list,
                pos_list,
                mask_list,
                ccd_list,
                SNR_list,
                RA_list,
                DEC_list,
                shape_std_max=outlier_std_max,
                print_fun=print_fun
            )

        mccd_star_list = []
        mccd_pos_list = []
        mccd_mask_list = []
        mccd_ccd_list = []
        mccd_SNR_list = []
        mccd_RA_list = []
        mccd_DEC_list = []

        for j in range(len(star_list)):
            # For each CCD
            if ccd_list[j] in CCD_id_filter_list:
                try:
                    n_stars = star_list[j].shape[2]

                    if n_stars >= min_n_stars:
                        mccd_star_list.append(star_list[j])
                        mccd_pos_list.append(pos_list[j])
                        mccd_mask_list.append(mask_list[j])
                        mccd_ccd_list.append(ccd_list[j] * np.ones(n_stars))
                        if SNR_list is not None:
                            mccd_SNR_list.append(SNR_list[j])
                        if RA_list is not None:
                            mccd_RA_list.append(RA_list[j])
                            mccd_DEC_list.append(DEC_list[j])
                    else:
                        msg = (
                            f"Not enough stars in catalog_id {catalog_id} "
                            + f",ccd {ccd_list[j]:d}. "
                            + f"Total stars = {n_stars:d}."
                        )
                        print_fun(msg)

                except Exception:
                    msg = (
                        f"Warning! Problem detected in catalog_id "
                        + f"{catalog_id} ,ccd {ccd_list[j]:d}"
                    )
                    print_fun(msg)

        if mccd_pos_list:
            # If the list is not empty
            # Concatenate, as fits can't handle list of numpy arrays and
            # turn into reg format
            mccd_stars = mccd.utils.reg_format(
                np.concatenate(mccd_star_list, axis=2))
            mccd_poss = np.concatenate(mccd_pos_list, axis=0)
            mccd_ccds = np.concatenate(mccd_ccd_list, axis=0)

            if save_masks is True:
                mccd_masks = mccd.utils.reg_format(
                    np.concatenate(mccd_mask_list, axis=2))
            else:
                # Send an array of False (None cannot be used in .fits)
                mccd_masks = np.zeros((mccd_poss.shape[0]), dtype=bool)

            if SNR_list is not None:
                mccd_SNRs = np.concatenate(mccd_SNR_list, axis=0)
            else:
                # Send an array of False (None cannot be used in .fits)
                mccd_SNRs = np.zeros((mccd_poss.shape[0]), dtype=bool)

            if RA_list is not None:
                mccd_RAs = np.concatenate(mccd_RA_list)
                mccd_DECs = np.concatenate(mccd_DEC_list)
            else:
                mccd_RAs = np.zeros((mccd_poss.shape[0]), dtype=bool)
                mccd_DECs = np.zeros((mccd_poss.shape[0]), dtype=bool)

            mccd_star_nb += mccd_stars.shape[0]

            # Save the fits file
            train_dic = {
                'VIGNET_LIST': mccd_stars,
                'GLOB_POSITION_IMG_LIST': mccd_poss,
                'MASK_LIST': mccd_masks,
                'CCD_ID_LIST': mccd_ccds,
                'SNR_WIN_LIST': mccd_SNRs,
                'RA_LIST': mccd_RAs,
                'DEC_LIST': mccd_DECs
            }

            saving_path = output_path + save_name + separator \
                + catalog_id + save_extension
            mccd.mccd_utils.save_to_fits(train_dic, saving_path)

    print_fun('Finished the training dataset processing.')
    print_fun(f"Total stars processed = {mccd_star_nb:d}")

    return mccd_inputs


def mccd_fit_pipeline(
    trainstar_path,
    file_number_string,
    mccd_parser,
    output_dir,
    verbose,
    saving_name='fitted_model',
    w_log=None
):
    r"""Fit the MCCD model to the observations."""
    # Extract the MCCD parameters from the parser
    mccd_inst_kw = mccd_parser.get_instance_kw()
    mccd_fit_kw = mccd_parser.get_fit_kw()
    use_SNR_weight = mccd_parser.get_extra_kw('use_SNR_weight')

    # Print the model configuration so that it is saved in log files
    w_log.info('MCCD configuration parameters:')
    w_log.info('[INPUTS]')
    inputs_dict_str = pprint.pformat({'use_SNR_weight': use_SNR_weight})
    w_log.info(inputs_dict_str)
    w_log.info('[INSTANCE]')
    inst_dict_str = pprint.pformat(mccd_inst_kw)
    w_log.info(inst_dict_str)
    w_log.info('[FIT]')
    fit_dict_str = pprint.pformat(mccd_fit_kw)
    w_log.info(fit_dict_str)
    w_log.info('End of MCCD configuration parameters.')

    # Open fits file
    starcat = fits.open(trainstar_path, memmap=False)

    mccd.auxiliary_fun.mccd_fit(
        starcat=starcat[1],
        mccd_inst_kw=mccd_inst_kw,
        mccd_fit_kw=mccd_fit_kw,
        output_dir=output_dir,
        catalog_id=file_number_string,
        sex_thresh=-1e5,
        use_SNR_weight=use_SNR_weight,
        verbose=verbose,
        saving_name=saving_name
    )

    starcat.close()


def mccd_validation_pipeline(
    teststar_path,
    mccd_model_path,
    mccd_parser,
    output_dir,
    file_number_string,
    w_log,
    val_saving_name
):
    r"""Validate the MCCD trained model against a set of observations."""
    w_log.info(f"Validating catalog {file_number_string}..")

    # Get MCCD parameters
    save_extension = '.fits'
    mccd_val_kw = mccd_parser.get_val_kw()
    testcat = fits.open(teststar_path, memmap=False)

    # Check if there is the fitted model
    if os.path.isfile(mccd_model_path):

        val_dict = mccd.auxiliary_fun.mccd_validation(
            mccd_model_path=mccd_model_path,
            testcat=testcat[1],
            **mccd_val_kw,
            sex_thresh=-1e5
        )

        testcat.close()

        val_saving_path = output_dir + val_saving_name + \
            file_number_string + save_extension

        # Save validation dictionary to fits file
        mccd.mccd_utils.save_to_fits(val_dict, val_saving_path)

        w_log.info(f"Validation catalog < {val_saving_path} > saved.")

    else:
        w_log.info(
            f"Fitted model corresponding to catalog"
            + f" {file_number_string} was not found."
        )


def mccd_interpolation_pipeline(
    mccd_model_path,
    galcat_path,
    pos_params,
    ccd_id,
    saving_path,
    get_shapes
):
    r"""Interpolate MCCD model."""
    # Import MCCD model
    mccd_model = mccd.mccd_quickload(mccd_model_path)
    # Open galaxy catalog
    galcat = fits.open(galcat_path, memmap=False)

    # Extract positions
    x_pos = galcat[2].data[pos_params[0]]
    y_pos = galcat[2].data[pos_params[1]]
    interp_pos = np.array([x_pos, y_pos]).T
    # Close catalog
    galcat.close()

    # Recover PSFs
    interp_PSFs = mccd_model.estimate_psf(test_pos=interp_pos, ccd_n=ccd_id)

    if interp_PSFs is not None:
        if get_shapes:
            PSF_moms = [
                galsim.hsm.FindAdaptiveMom(galsim.Image(psf), strict=False)
                for psf in interp_PSFs
            ]

            PSF_shapes = np.array([
                [
                    moms.observed_shape.g1,
                    moms.observed_shape.g2,
                    moms.moments_sigma,
                    int(bool(moms.error_message))
                ]
                for moms in PSF_moms
            ])

        shapepipe_write_output(
            saving_path=saving_path,
            example_fits_path=galcat_path,
            get_shapes=get_shapes,
            interp_PSFs=interp_PSFs,
            PSF_shapes=PSF_shapes
        )

        return None
    else:
        return NOT_ENOUGH_STARS


def shapepipe_write_output(
    saving_path,
    example_fits_path,
    get_shapes,
    interp_PSFs,
    PSF_shapes=None
):
    r"""Save interpolated PSFs dictionary to fits file.

    The saved files are compatible with the previous shapepipe's standard.
    """
    output = sc.FITSCatalog(
        saving_path,
        open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
        SEx_catalog=True
    )

    if get_shapes:
        data = {
            'VIGNET': interp_PSFs,
            'E1_PSF_HSM': PSF_shapes[:, 0],
            'E2_PSF_HSM': PSF_shapes[:, 1],
            'SIGMA_PSF_HSM': PSF_shapes[:, 2],
            'FLAG_PSF_HSM': PSF_shapes[:, 3].astype(int)
        }
    else:
        data = {'VIGNET': interp_PSFs}

    output.save_as_fits(data, sex_cat_path=example_fits_path)
