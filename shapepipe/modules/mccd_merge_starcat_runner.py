# -*- coding: utf-8 -*-

"""MCCD MERGE STARCAT RUNNER.

This module is used to merge the validation stars results of the MCCD
validation runner or fit_validation runner.

:Author: Tobias Liaudat

"""

import numpy as np
from astropy.io import fits
from shapepipe.pipeline import file_io as sc
from shapepipe.modules.module_decorator import module_runner


@module_runner(
    input_module=['mccd_fit_val_runner', 'mccd_val_runner'],
    version='1.0',
    file_pattern=['validation_psf'],
    file_ext=['.fits'],
    numbering_scheme='-0000000',
    depends=['numpy', 'astropy'],
    run_method='serial'
)
def mccd_merge_starcat_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log
):

    w_log.info('Merging validation results..')
    hdu_table = 1

    output_dir = run_dirs['output']

    x, y = [], []
    ra, dec = [], []
    g1_psf, g2_psf, size_psf = [], [], []
    g1, g2, size = [], [], []
    flag_psf, flag_star = [], []
    ccd_nb = []
    pixel_mse = []
    pixel_sum = []
    masked_pixel_mse = []
    masked_pixel_sum = []
    size_mse = []
    masked_size = []
    pix_norm_mse, size_norm_mse = [], []
    pix_filt_mse, size_filt_mse = [], []

    star_noise_var, star_noise_size = [], []
    model_var, model_var_size = [], []

    bad_catalogs = 0

    # Construction of the mask
    shap = np.array([51, 51])
    cent = np.array([25, 25])
    rad = 10

    my_mask = np.zeros((51, 51), dtype=bool)

    for i in range(0, shap[0]):
        for j in range(0, shap[1]):
            if np.sqrt((i - cent[0]) ** 2 + (j - cent[1]) ** 2) <= rad:
                my_mask[i, j] = True

    for name in input_file_list:
        starcat_j = fits.open(name[0], memmap=False)

        stars = np.copy(starcat_j[hdu_table].data['VIGNET_LIST'])
        stars[stars < -1e6] = 0
        psfs = np.copy(starcat_j[hdu_table].data['PSF_VIGNET_LIST'])

        # Pixel mse calculation
        pix_val = np.sum((stars - psfs) ** 2)
        pix_sum = np.sum((stars - psfs))
        masked_diffs = np.array(
            [(_star - _psf)[my_mask] for _star, _psf in zip(stars, psfs)]
        )
        masked_pix_val = np.sum(masked_diffs ** 2)
        masked_pix_sum = np.sum(masked_diffs)

        # Star noise variance (using masked stars)
        star_noise_var_val = np.array(
            [np.var(_star[~my_mask]) for _star in stars]
        )
        res_var_val = np.array(
            [np.var(_star - _psf) for _star, _psf in zip(stars, psfs)]
        )

        # Variance of the model
        # (Residual variance  - Star variance (using masked stars))
        model_var_val = res_var_val - star_noise_var_val
        model_var_val = model_var_val[model_var_val > 0]

        # if pix_val < 1e20:
        # Normalised pixel mse calculation
        stars_norm_vals = np.sqrt(np.sum(stars ** 2, axis=(1, 2)))
        psfs_norm_vals = np.sqrt(np.sum(psfs ** 2, axis=(1, 2)))
        # Select non zero stars & psfs
        non_zero_elems = np.logical_and(
            (psfs_norm_vals != 0),
            (stars_norm_vals != 0)
        )
        # Calculate the filtered mse calculation
        pix_filt_val = np.sum(
            (stars[non_zero_elems] - psfs[non_zero_elems]) ** 2
        )
        # Calculate the normalized (& filtered) mse calculation
        stars_norm_vals = stars_norm_vals[non_zero_elems].reshape(-1, 1, 1)
        psfs_norm_vals = psfs_norm_vals[non_zero_elems].reshape(-1, 1, 1)
        pix_norm_val = np.sum(
            (
                stars[non_zero_elems] / stars_norm_vals -
                psfs[non_zero_elems] / psfs_norm_vals
            ) ** 2
        )
        # Calculate sizes
        filt_size = stars[non_zero_elems].size
        regular_size = stars.size
        regular_masked_size = stars.shape[0] * (np.sum(my_mask))

        # Append the results to the lists
        pixel_mse.append(pix_val)
        pixel_sum.append(pix_sum)
        masked_pixel_mse.append(masked_pix_val)
        masked_pixel_sum.append(masked_pix_sum)
        size_mse.append(regular_size)
        masked_size.append(regular_masked_size)

        pix_norm_mse.append(pix_norm_val)
        size_norm_mse.append(filt_size)
        pix_filt_mse.append(pix_filt_val)
        size_filt_mse.append(filt_size)

        star_noise_var.append(star_noise_var_val)
        star_noise_size.append(star_noise_var_val.size)
        model_var.append(model_var_val)
        model_var_size.append(model_var_val.size)

        # positions
        x += list(starcat_j[hdu_table].data['GLOB_POSITION_IMG_LIST'][:, 0])
        y += list(starcat_j[hdu_table].data['GLOB_POSITION_IMG_LIST'][:, 1])

        # RA and DEC positions
        try:
            ra += list(starcat_j[hdu_table].data['RA_LIST'][:])
            dec += list(starcat_j[hdu_table].data['DEC_LIST'][:])
        except Exception:
            ra += list(np.zeros(
                starcat_j[hdu_table].data[
                    'GLOB_POSITION_IMG_LIST'
                ][:, 0].shape,
                dtype=int,
            ))
            dec += list(np.zeros(
                starcat_j[hdu_table].data[
                    'GLOB_POSITION_IMG_LIST'
                ][:, 0].shape,
                dtype=int,
            ))

        # shapes (convert sigmas to R^2)
        g1_psf += list(starcat_j[hdu_table].data['PSF_MOM_LIST'][:, 0])
        g2_psf += list(starcat_j[hdu_table].data['PSF_MOM_LIST'][:, 1])
        size_psf += list(starcat_j[hdu_table].data['PSF_MOM_LIST'][:, 2] ** 2)
        g1 += list(starcat_j[hdu_table].data['STAR_MOM_LIST'][:, 0])
        g2 += list(starcat_j[hdu_table].data['STAR_MOM_LIST'][:, 1])
        size += list(starcat_j[hdu_table].data['STAR_MOM_LIST'][:, 2] ** 2)

        # flags
        flag_psf += list(starcat_j[hdu_table].data['PSF_MOM_LIST'][:, 3])
        flag_star += list(starcat_j[hdu_table].data['STAR_MOM_LIST'][:, 3])

        # ccd id list
        ccd_nb += list(starcat_j[hdu_table].data['CCD_ID_LIST'])

        starcat_j.close()

    # ----- #
    # Stats about the merging
    def rmse_calc(values, sizes):
        return np.sqrt(
            np.nansum(np.array(values)) / np.nansum(np.array(sizes))
        )

    def rmse_calc_2(values, sizes):
        return np.sqrt(
            np.nansum(np.array(values) ** 2) / np.nansum(np.array(sizes))
        )

    def mean_calc(values, sizes):
        return (
            np.nansum(np.array(values)) / np.nansum(np.array(sizes))
        )

    def std_calc(values):
        return np.nanstd(np.array(values))

    # Regular pixel RMSE
    tot_pixel_rmse = rmse_calc(pixel_mse, size_mse)
    w_log.info(
        f"MCCD_merge_starcat: Regular Total pixel RMSE ="
        + f" {tot_pixel_rmse:.5e}\n"
    )

    # Regular Total pixel mean
    tot_pixel_mean = mean_calc(pixel_sum, size_mse)
    w_log.info(
        f"MCCD_merge_starcat: Regular Total pixel mean ="
        + f" {tot_pixel_mean:.5e}\n"
    )

    # Regular Total MASKED pixel RMSE
    tot_masked_pixel_rmse = rmse_calc(masked_pixel_mse, masked_size)
    w_log.info(
        f"MCCD_merge_starcat: Regular Total MASKED pixel RMSE ="
        + f" {tot_masked_pixel_rmse:.5e}\n"
    )

    # Regular Total MASKED pixel mean
    tot_masked_pixel_mean = mean_calc(masked_pixel_sum, masked_size)
    w_log.info(
        f"MCCD_merge_starcat: Regular Total MASKED pixel mean ="
        + f" {tot_masked_pixel_mean:.5e}\n"
    )

    # Normalized pixel RMSE
    tot_pix_norm_rmse = rmse_calc(pix_norm_mse, size_norm_mse)
    w_log.info(
        "MCCD_merge_starcat: Normalized Total pixel RMSE ="
        + f" {tot_pix_norm_rmse:.5e}\n"
    )

    # Normalized filtered pixel RMSE
    tot_pix_filt_rmse = rmse_calc(pix_filt_mse, size_filt_mse)
    w_log.info(
        "MCCD_merge_starcat: Filtered Total pixel RMSE ="
        + f" {tot_pix_filt_rmse:.5e}\n"
    )

    concat_model_var = np.concatenate(np.array(model_var))
    concat_star_noise_var = np.concatenate(np.array(star_noise_var))

    # Model variance
    mean_model_var = mean_calc(concat_model_var, model_var_size)
    std_model_var = std_calc(concat_model_var)
    rmse_model_var = rmse_calc_2(concat_model_var, model_var_size)
    w_log.info(
        f"MCCD-RCA variance:\n Mean Variance= {mean_model_var:.5e}\n"
        + f"Std Variance= {std_model_var:.5e}\n"
        + f"RMSE Variance= {rmse_model_var:.5e}\n"
    )

    # Star Noise Variance
    mean_star_var = mean_calc(concat_star_noise_var, star_noise_size)
    std_star_var = std_calc(concat_star_noise_var)
    rmse_star_var = rmse_calc_2(concat_star_noise_var, star_noise_size)
    w_log.info(
        f"Masked stars variance:\n Mean Variance= {mean_star_var:.5e}\n"
        + f"Std Variance= {std_star_var:.5e}\n"
        + f"RMSE Variance= {rmse_star_var:.5e}\n"
    )

    # Mask and transform to numpy arrays
    flagmask = np.abs(np.array(flag_star) - 1) * np.abs(np.array(flag_psf) - 1)
    psf_e1 = np.array(g1_psf)[flagmask.astype(bool)]
    psf_e2 = np.array(g2_psf)[flagmask.astype(bool)]
    psf_r2 = np.array(size_psf)[flagmask.astype(bool)]
    star_e1 = np.array(g1)[flagmask.astype(bool)]
    star_e2 = np.array(g2)[flagmask.astype(bool)]
    star_r2 = np.array(size)[flagmask.astype(bool)]

    rmse, mean, std_dev = stats_calculator(star_e1, psf_e1)
    w_log.info(
        f"Moment residual e1:\n Mean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
        + f" RMSE= {rmse:.5e}\n"
    )

    rmse, mean, std_dev = stats_calculator(star_e2, psf_e2)
    w_log.info(
        f"Moment residual e2:\n Mean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
        + f"RMSE= {rmse:.5e}\n"
    )

    rmse, mean, std_dev = stats_calculator(star_r2, psf_r2)
    w_log.info(
        f"Moment residual R2:\n Mean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
        + f"RMSE= {rmse:.5e}\n"
    )

    print(f"MCCD: Number of stars: {star_e1.shape[0]:d}")
    w_log.info(f"MCCD: Number of stars: {star_e1.shape[0]:d}")

    output = sc.FITSCatalog(
        output_dir + '/full_starcat-0000000.fits',
        open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
        SEx_catalog=True
    )

    # convert back to sigma for consistency
    data = {
        'X': x,
        'Y': y,
        'RA': ra,
        'DEC': dec,
        'E1_PSF_HSM': g1_psf,
        'E2_PSF_HSM': g2_psf,
        'SIGMA_PSF_HSM': np.sqrt(size_psf),
        'E1_STAR_HSM': g1,
        'E2_STAR_HSM': g2,
        'SIGMA_STAR_HSM': np.sqrt(size),
        'FLAG_PSF_HSM': flag_psf,
        'FLAG_STAR_HSM': flag_star,
        'CCD_NB': ccd_nb
    }
    print('Writing full catalog...')

    output.save_as_fits(data, sex_cat_path=input_file_list[0][0])
    print('... Done.')

    return None, None


def stats_calculator(val_ref, val_model):
    residual = val_ref - val_model

    rmse = np.sqrt(np.mean(residual ** 2))
    mean = np.mean(residual)
    std_dev = np.std(residual)

    return rmse, mean, std_dev
