# -*- coding: utf-8 -*-
"""MCCD MERGE STARCAT SCRIPT

This module contains a class to identify single exposures that were used
to create tiles.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>, Tobias Liaudat

:Date: 2021

:Package: ShapePipe

"""


import numpy as np
from astropy.io import fits
from shapepipe.pipeline import file_io as sc


class MergeStarCat(object):
    """Merge Star Catalogue

    Merge MCCD star catalogues

    Parameters
    ----------
    input_file_list : list of str
        input files
    output_dir : str
        output directory
    w_log :
        log file
    stamp_size : int
        stamp size, in pixels
    rad : int
        radius for mask, in pixels, required for
        some statistics computations
    hdu_table : int, optional, default=1
        HDU number
    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        w_log,
        stamp_size,
        rad,
        hdu_table=1
    ):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log
        self._stamp_size = stamp_size
        self._rad = rad
        self._hdu_table = hdu_table

    @staticmethod
    def rmse_calc(values, sizes):
        """Rmse Calculation

        Calculate square root of mean over input values. If values are (x - <x>),
        this function computes the RMSE.

        Parameters
        ----------
        values : list of float
            sums of pixel values for a list of input images
        sizes : list of int
            number of pixels for a list of input images

        Returns
        -------
        rmse : float
            root mean square error
        """

        rmse = np.sqrt(MergeStarCat.mean_calc(values, sizes))

        return rmse

    @staticmethod
    def rmse_calc_2(values, sizes):
        """Rmse calculation 2

        TBD

        Parameters
        ----------
        values : list of float
            sums of pixel values for a list of input images
        sizes : list of int
            number of pixels for a list of input images

        Returns
        -------
        rmse : float
            root mean square error
        """

        rmse = np.sqrt(
            np.nansum(np.array(values) ** 2) / np.nansum(np.array(sizes))
        )

        return rmse

    @staticmethod
    def mean_calc(values, sizes):
        """"Mean calculation

        Calculate pixel mean over all input images

        Parameters
        ----------
        values : list of float
            sums of pixel values for a list of input images
        sizes : list of int
            number of pixels for a list of input images

        Returns
        -------
        mean : float
            mean
        """

        mean = (
            np.nansum(np.array(values)) / np.nansum(np.array(sizes))
        )

        return mean

    @staticmethod
    def std_calc(values):
        """"Standard deviation calculation

        Calculate pixel standard deviation over all input images

        Parameters
        ----------
        values : list of float
            sums of pixel values for a list of input images
        sizes : list of int
            number of pixels for a list of input images

        Returns
        -------
        std : float
            standard deviation
        """

        std = np.nanstd(np.array(values))

        return std

    @staticmethod
    def stats_calculator(val_ref, val_model):
        """Stats Calculator

        Calculate RMSE, mean, and standard deviation of residuals

        Parameters
        ----------
        val_ref : list of float
            reference values
        val_model : list of flost
            model values

        Returns
        -------
        rmse : float
            root mean square error
        mean : float
            mean
        std_dev : float
            standard deviation
        """

        residual = val_ref - val_model

        rmse = np.sqrt(np.mean(residual ** 2))
        mean = np.mean(residual)
        std_dev = np.std(residual)

        return rmse, mean, std_dev

    def process(self):

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
        shap = np.array([self._stamp_size, self._stamp_size])
        stamp_size_half = int(self._stamp_size / 2)
        cent = np.array([stamp_size_half, stamp_size_half])

        my_mask = np.zeros((self._stamp_size, self._stamp_size), dtype=bool)

        for i in range(0, shap[0]):
            for j in range(0, shap[1]):
                if np.sqrt((i - cent[0]) ** 2 + (j - cent[1]) ** 2) <= self._rad:
                    my_mask[i, j] = True

        for name in self._input_file_list:
            starcat_j = fits.open(name[0], memmap=False)

            stars = np.copy(starcat_j[self._hdu_table].data['VIGNET_LIST'])
            stars[stars < -1e6] = 0
            psfs = np.copy(starcat_j[self._hdu_table].data['PSF_VIGNET_LIST'])

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
            x += list(
                starcat_j[self._hdu_table].data['GLOB_POSITION_IMG_LIST'][:, 0]
            )
            y += list(
                starcat_j[self._hdu_table].data['GLOB_POSITION_IMG_LIST'][:, 1]
            )

            # RA and DEC positions
            try:
                ra += list(starcat_j[self._hdu_table].data['RA_LIST'][:])
                dec += list(starcat_j[self._hdu_table].data['DEC_LIST'][:])
            except Exception:
                ra += list(np.zeros(
                    starcat_j[self._hdu_table].data[
                        'GLOB_POSITION_IMG_LIST'
                    ][:, 0].shape,
                    dtype=int,
                ))
                dec += list(np.zeros(
                    starcat_j[self._hdu_table].data[
                        'GLOB_POSITION_IMG_LIST'
                    ][:, 0].shape,
                    dtype=int,
                ))

            # shapes (convert sigmas to R^2)
            g1_psf += list(
                starcat_j[self._hdu_table].data['PSF_MOM_LIST'][:, 0]
            )
            g2_psf += list(
                starcat_j[self._hdu_table].data['PSF_MOM_LIST'][:, 1]
            )
            size_psf += list(
                starcat_j[self._hdu_table].data['PSF_MOM_LIST'][:, 2] ** 2
            )
            g1 += list(starcat_j[self._hdu_table].data['STAR_MOM_LIST'][:, 0])
            g2 += list(starcat_j[self._hdu_table].data['STAR_MOM_LIST'][:, 1])
            size += list(
                starcat_j[self._hdu_table].data['STAR_MOM_LIST'][:, 2] ** 2
            )

            # flags
            flag_psf += list(
                starcat_j[self._hdu_table].data['PSF_MOM_LIST'][:, 3]
            )
            flag_star += list(
                starcat_j[self._hdu_table].data['STAR_MOM_LIST'][:, 3]
            )

            # ccd id list
            ccd_nb += list(starcat_j[self._hdu_table].data['CCD_ID_LIST'])

            starcat_j.close()

        # Regular pixel RMSE
        tot_pixel_rmse = MergeStarCat.rmse_calc(pixel_mse, size_mse)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total pixel RMSE ='
            + f' {tot_pixel_rmse:.5e}\n'
        )

        # Regular Total pixel mean
        tot_pixel_mean = MergeStarCat.mean_calc(pixel_sum, size_mse)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total pixel mean ='
            + f' {tot_pixel_mean:.5e}\n'
        )

        # Regular Total MASKED pixel RMSE
        tot_masked_pixel_rmse = MergeStarCat.rmse_calc(
            masked_pixel_mse, masked_size
        )
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total MASKED pixel RMSE ='
            + f' {tot_masked_pixel_rmse:.5e}\n'
        )

        # Regular Total MASKED pixel mean
        tot_masked_pixel_mean = MergeStarCat.mean_calc(
            masked_pixel_sum, masked_size
        )
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total MASKED pixel mean ='
            + f' {tot_masked_pixel_mean:.5e}\n'
        )

        # Normalized pixel RMSE
        tot_pix_norm_rmse = MergeStarCat.rmse_calc(pix_norm_mse, size_norm_mse)
        self._w_log.info(
            'MCCD_merge_starcat: Normalized Total pixel RMSE ='
            + f' {tot_pix_norm_rmse:.5e}\n'
        )

        # Normalized filtered pixel RMSE
        tot_pix_filt_rmse = MergeStarCat.rmse_calc(pix_filt_mse, size_filt_mse)
        self._w_log.info(
            'MCCD_merge_starcat: Filtered Total pixel RMSE ='
            + f' {tot_pix_filt_rmse:.5e}\n'
        )

        concat_model_var = np.concatenate(np.array(model_var))
        concat_star_noise_var = np.concatenate(np.array(star_noise_var))

        # Model variance
        mean_model_var = MergeStarCat.mean_calc(
            concat_model_var, model_var_size
        )
        std_model_var = MergeStarCat.std_calc(concat_model_var)
        rmse_model_var = MergeStarCat.rmse_calc_2(
            concat_model_var, model_var_size
        )
        self._w_log.info(
            f'MCCD-RCA variance:\nMean Variance= {mean_model_var:.5e}\n'
            + f'Std Variance= {std_model_var:.5e}\n'
            + f'RMSE Variance= {rmse_model_var:.5e}\n'
        )

        # Star Noise Variance
        mean_star_var = MergeStarCat.mean_calc(
            concat_star_noise_var, star_noise_size
        )
        std_star_var = MergeStarCat.std_calc(concat_star_noise_var)
        rmse_star_var = MergeStarCat.rmse_calc_2(
            concat_star_noise_var, star_noise_size
        )
        self._w_log.info(
            f'Masked stars variance:\nMean Variance= {mean_star_var:.5e}\n'
            + f'Std Variance= {std_star_var:.5e}\n'
            + f'RMSE Variance= {rmse_star_var:.5e}\n'
        )

        # Mask and transform to numpy arrays
        flagmask = (
            np.abs(np.array(flag_star) - 1) * np.abs(np.array(flag_psf) - 1)
        )
        psf_e1 = np.array(g1_psf)[flagmask.astype(bool)]
        psf_e2 = np.array(g2_psf)[flagmask.astype(bool)]
        psf_r2 = np.array(size_psf)[flagmask.astype(bool)]
        star_e1 = np.array(g1)[flagmask.astype(bool)]
        star_e2 = np.array(g2)[flagmask.astype(bool)]
        star_r2 = np.array(size)[flagmask.astype(bool)]

        rmse, mean, std_dev = MergeStarCat.stats_calculator(star_e1, psf_e1)
        self._w_log.info(
            f'Moment residual e1:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        rmse, mean, std_dev = MergeStarCat.stats_calculator(star_e2, psf_e2)
        self._w_log.info(
            f'Moment residual e2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        rmse, mean, std_dev = MergeStarCat.stats_calculator(star_r2, psf_r2)
        self._w_log.info(
            f'Moment residual R2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        self._w_log.info(f'MCCD: Number of stars: {star_e1.shape[0]:d}')

        output = sc.FITSCatalog(
            self._output_dir + '/full_starcat-0000000.fits',
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

        output.save_as_fits(data, sex_cat_path=self._input_file_list[0][0])
