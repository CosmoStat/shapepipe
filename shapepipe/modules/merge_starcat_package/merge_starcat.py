"""MERGE STAR CATALOGUES.

This module contains a class to identify single exposures that were used
to create tiles.

:Authors: Martin Kilbinger <martin.kilbinger@cea.fr>, Tobias Liaudat,
    Morgan Schmitz, Axel Guinot

"""

import re

import numpy as np
from astropy.io import fits

from shapepipe.pipeline import file_io


class MergeStarCatMCCD(object):
    """Merge Star Catalogue MCCD.

    Merge star catalogues of MCCD PSF model output.

    Parameters
    ----------
    input_file_list : list
        Input files
    output_dir : str
        Output directory
    w_log : logging.Logger
        Logging instance
    stamp_size : int, optional
        Stamp size, in pixels; default is ``51``
    rad : int, optional
        Radius for mask, in pixels; default is ``10``
    hdu_table : int, optional
        HDU number; default is ``1``

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        w_log,
        stamp_size=51,
        rad=10,
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
        r"""Calculate RMSE.

        Calculate square root of mean over input values.
        If ``values`` is an array with element :math:`j` being
        :math:`\sum_j^{N_j}x_{i, j}^2`, where :math:`x_{ij}`
        is the residual (ground truth - estimation), and
        sizes is the array :math:`N_j`, then
        this function computes the RMSE.

        Parameters
        ----------
        values : list
            Sums of pixel values for a list of input images
        sizes : list
            Number of pixels for a list of input images

        Returns
        -------
        rmse : float
            Root mean square error

        See Also
        --------
        MergeStarCatMCCD.mean_calc

        """
        rmse = np.sqrt(MergeStarCatMCCD.mean_calc(values, sizes))

        return rmse

    @staticmethod
    def rmse_calc_2(values, sizes):
        r"""Calculate RMSE 2.

        Calculate square root of mean over squared input valuess.
        If ``values`` is an array with element :math:`j` being
        :math:`\sum_j^{N_j}x_{i, j}`, where :math:`x_{ij}`
        is the residual (ground truth - estimation), and
        sizes is the array :math:`N_j`, then
        this function computes the RMSE.

        Parameters
        ----------
        values : list
            Sums of pixel values for a list of input images
        sizes : list
            Number of pixels for a list of input images

        Returns
        -------
        float
            Root mean square error

        """
        rmse = np.sqrt(
            np.nansum(np.array(values) ** 2) / np.nansum(np.array(sizes))
        )

        return rmse

    @staticmethod
    def mean_calc(values, sizes):
        """Calculate Mean.

        Calculate pixel mean over all input images.

        Parameters
        ----------
        values : list
            Sums of pixel values for a list of input images
        sizes : list
            Number of pixels for a list of input images

        Returns
        -------
        float
            Mean

        """
        mean = (
            np.nansum(np.array(values)) / np.nansum(np.array(sizes))
        )

        return mean

    @staticmethod
    def std_calc(values):
        """Calculate Standard Deviation.

        Calculate pixel standard deviation over all input images.

        Parameters
        ----------
        values : list
            Sums of pixel values for a list of input images
        sizes : list
            Number of pixels for a list of input images

        Returns
        -------
        float
            Standard deviation

        """
        std = np.nanstd(np.array(values))

        return std

    @staticmethod
    def stats_calculator(val_ref, val_model):
        """Calculate Stats.

        Calculate RMSE, mean, and standard deviation of residuals.

        Parameters
        ----------
        val_ref : list
            Reference values
        val_model : list
            Model values

        Returns
        -------
        tuple
            Root mean square error, mean and standard deviation

        """
        residual = val_ref - val_model

        rmse = np.sqrt(np.mean(residual ** 2))
        mean = np.mean(residual)
        std_dev = np.std(residual)

        return rmse, mean, std_dev

    def process(self):
        """Process.

        Process merging.

        """
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

        idx = np.arange(0, shap[0])
        jdx = np.arange(0, shap[1])
        inside_circle = np.sqrt(
            (idx[np.newaxis, :] - cent[0]) ** 2
            + (jdx[:, np.newaxis] - cent[1]) ** 2
        ) <= self._rad
        my_mask[inside_circle] = True

        for name in self._input_file_list:
            starcat_j = fits.open(name[0], memmap=False)

            try:
                stars = np.copy(starcat_j[self._hdu_table].data['VIGNET_LIST'])
            except ValueError:
                print(f'Error for file {name[0]}, check FITS file integrity')
                raise
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
                [np.var(_star[np.invert(my_mask)]) for _star in stars]
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
                    stars[non_zero_elems] / stars_norm_vals
                    - psfs[non_zero_elems] / psfs_norm_vals
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

        # Shortcut name
        MSC = MergeStarCatMCCD

        # Regular pixel RMSE
        tot_pixel_rmse = MSC.rmse_calc(pixel_mse, size_mse)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total pixel RMSE ='
            + f' {tot_pixel_rmse:.5e}\n'
        )

        # Regular Total pixel mean
        tot_pixel_mean = MSC.mean_calc(pixel_sum, size_mse)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total pixel mean ='
            + f' {tot_pixel_mean:.5e}\n'
        )

        # Regular Total MASKED pixel RMSE
        tot_masked_pixel_rmse = MSC.rmse_calc(masked_pixel_mse, masked_size)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total MASKED pixel RMSE ='
            + f' {tot_masked_pixel_rmse:.5e}\n'
        )

        # Regular Total MASKED pixel mean
        tot_masked_pixel_mean = MSC.mean_calc(masked_pixel_sum, masked_size)
        self._w_log.info(
            f'MCCD_merge_starcat: Regular Total MASKED pixel mean ='
            + f' {tot_masked_pixel_mean:.5e}\n'
        )

        # Normalized pixel RMSE
        tot_pix_norm_rmse = MSC.rmse_calc(pix_norm_mse, size_norm_mse)
        self._w_log.info(
            'MCCD_merge_starcat: Normalized Total pixel RMSE ='
            + f' {tot_pix_norm_rmse:.5e}\n'
        )

        # Normalized filtered pixel RMSE
        tot_pix_filt_rmse = MSC.rmse_calc(pix_filt_mse, size_filt_mse)
        self._w_log.info(
            'MCCD_merge_starcat: Filtered Total pixel RMSE ='
            + f' {tot_pix_filt_rmse:.5e}\n'
        )

        concat_model_var = np.concatenate(np.array(model_var))
        concat_star_noise_var = np.concatenate(np.array(star_noise_var))

        # Model variance
        mean_model_var = MSC.mean_calc(concat_model_var, model_var_size)
        std_model_var = MSC.std_calc(concat_model_var)
        rmse_model_var = MSC.rmse_calc_2(concat_model_var, model_var_size)
        self._w_log.info(
            f'MCCD-RCA variance:\nMean Variance= {mean_model_var:.5e}\n'
            + f'Std Variance= {std_model_var:.5e}\n'
            + f'RMSE Variance= {rmse_model_var:.5e}\n'
        )

        # Star Noise Variance
        mean_star_var = MSC.mean_calc(concat_star_noise_var, star_noise_size)
        std_star_var = MSC.std_calc(concat_star_noise_var)
        rmse_star_var = MSC.rmse_calc_2(concat_star_noise_var, star_noise_size)
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

        rmse, mean, std_dev = MSC.stats_calculator(star_e1, psf_e1)
        self._w_log.info(
            f'Moment residual e1:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        rmse, mean, std_dev = MSC.stats_calculator(star_e2, psf_e2)
        self._w_log.info(
            f'Moment residual e2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        rmse, mean, std_dev = MSC.stats_calculator(star_r2, psf_r2)
        self._w_log.info(
            f'Moment residual R2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n'
            + f'RMSE= {rmse:.5e}\n'
        )

        self._w_log.info(f'MCCD: Number of stars: {star_e1.shape[0]:d}')

        # Prepare output FITS catalogue
        output = file_io.FITSCatalogue(
            f'{self._output_dir}/full_starcat-0000000.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True
        )

        # Collect columns
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

        # Write file
        output.save_as_fits(data, sex_cat_path=self._input_file_list[0][0])


class MergeStarCatPSFEX(object):
    """Merge Star Catalogue PSFEx.

    Merge star catalogues of PSFEx PSF model output.

    Parameters
    ----------
    input_file_list : list
        Input files
    output_dir : str
        Output directory
    w_log : logging.Logger
        Logging instance
    hdu_table : int, optional
        HDU number; default is ``2``

    """

    def __init__(self, input_file_list, output_dir, w_log, hdu_table=2):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log
        self._hdu_table = hdu_table

    def process(self):
        """Process.

        Process merging.

        """
        x, y, ra, dec = [], [], [], []
        g1_psf, g2_psf, size_psf = [], [], []
        g1, g2, size = [], [], []
        flag_psf, flag_star = [], []
        mag, snr, psfex_acc = [], [], []
        ccd_nb = []

        self._w_log.info(
            f'Merging {len(self._input_file_list)} star catalogues'
        )

        for name in self._input_file_list:
            starcat_j = fits.open(name[0], memmap=False)

            data_j = starcat_j[self._hdu_table].data

            # positions
            x += list(data_j['X'])
            y += list(data_j['Y'])
            ra += list(data_j['RA'])
            dec += list(data_j['DEC'])

            # shapes (convert sigmas to R^2)
            g1_psf += list(data_j['E1_PSF_HSM'])
            g2_psf += list(data_j['E2_PSF_HSM'])
            size_psf += list(data_j['SIGMA_PSF_HSM']**2)
            g1 += list(data_j['E1_STAR_HSM'])
            g2 += list(data_j['E2_STAR_HSM'])
            size += list(data_j['SIGMA_STAR_HSM']**2)

            # flags
            flag_psf += list(data_j['FLAG_PSF_HSM'])
            flag_star += list(data_j['FLAG_STAR_HSM'])

            # misc
            mag += list(data_j['MAG'])
            snr += list(data_j['SNR'])
            psfex_acc += list(data_j['ACCEPTED'])

            # CCD number
            ccd_nb += [
                re.split(r"\-([0-9]*)\-([0-9]+)\.", name[0])[-2]
            ] * len(data_j['RA'])

        # Prepare output FITS catalogue
        output = file_io.FITSCatalogue(
            f'{self._output_dir}/full_starcat-0000000.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True
        )

        # Collect columns
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
            'MAG': mag,
            'SNR': snr,
            'ACCEPTED': psfex_acc,
            'CCD_NB': ccd_nb
        }

        # Write file
        output.save_as_fits(
            data,
            overwrite=True,
            sex_cat_path=self._input_file_list[0][0],
        )


class MergeStarCatSetools(object):
    """Merge Star Catalogue Setools.

    Merge star catalogues of Setools output.

    Parameters
    ----------
    input_file_list : list
        Input files
    output_dir : str
        Output directory
    w_log : logging.Logger
        Logging instance
    hdu_table : int, optional
        HDU number; default is ``2``

    """

    def __init__(self, input_file_list, output_dir, w_log, hdu_table=2):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log
        self._hdu_table = hdu_table

    @classmethod
    def get_moments(cls, data):
        """Get Moments.

        Return second-order moments.

        Parameters
        ----------
        data : dict
            input data

        Returns
        -------
        m11 : float
            second-order moment along xy
        m20 : float
            second-order moment along x
        m02 : float
            second-order moment along y

        """
        # SExtractor output. First and second moments are normalised.
        # Second moments are centred.
        q11 = 'X2WIN_IMAGE'
        q22 = 'Y2WIN_IMAGE'
        q12 = 'XYWIN_IMAGE'

        # Second moments
        m11 = data[q12]
        m20 = data[q11]
        m02 = data[q22]

        return m11, m20, m02

    @classmethod
    def get_ellipticity(cls, m11, m20, m02, typ):
        """Get Ellipticity.

        Compute ellipticity from second-order moments.

        Parameters
        ----------
        m11 : float
            second-order moment along xy
        m20 : float
            second-order moment along x
        m02 : float
            second-order moment along y
        typ : str
            ellipticity type, allowed are 'epsilon', 'chi'

        Returns
        -------
        list
            ellipticity components

        """
        if typ == 'epsilon':
            # Determinant = (Q_11 Q_22 - Q_12^2)^(1/2)
            det = np.sqrt(m20 * m02 - m11 * m11)
        elif typ == 'chi':
            det = 0
        else:
            raise ValueError(f'Invalid ellipticity type {type}')

        # Denominator = Q_11 + Q_22 [ + 2 * det]
        den = m20 + m02 + 2 * det

        # Ellipticity = (Q_11 - Q_22 + 2 i Q_12) / den
        ell = (m20 - m02 + 1j * 2 * m11) / den

        if type == 'chi':
            # chi estimates 2*g, so to get g we have to divide by 2
            ell = ell / 2

        return ell.real, ell.imag

    def process(self):
        """Process.

        Process merging.

        """
        x, y, ra, dec = [], [], [], []
        eps1, ep2, chi1, chi2, size = [], [], [], [], []
        flags, flags_ext = [], []
        mag, snr = [], []
        ccd_nb = []

        self._w_log.info(
            f'Merging {len(self._input_file_list)} star catalogues'
        )

        for name in self._input_file_list:
            starcat_j = fits.open(name[0], memmap=False)

            data_j = starcat_j[self._hdu_table].data

            # positions
            x += list(data_j['XWIN_IMAGE'])
            y += list(data_j['YWIN_IMAGE'])
            ra += list(data_j['XWIN_WORLD'])
            dec += list(data_j['YWIN_WORLD'])

            m11, m20, m02 = self.get_moments(data_j)
            eps1, eps2 = self.get_ellipticity(m11, m20, m02, 'epsilon')
            chi1, chi2 = self.get_ellipticity(m11, m20, m02, 'chi')

            size += list(data_j['FLUX_RADIUS'])

            # flags
            flags += list(data_j['FLAGS_WIN'])
            flags_ext += list(data_j['IMAFLAGS_ISO'])

            # misc
            mag += list(data_j['MAG_WIN'])
            snr += list(data_j['SNR_WIN'])

            # CCD number
            ccd_nb += [
                re.split(r"\-([0-9]*)\-([0-9]+)\.", name[0])[-2]
            ] * len(data_j['XWIN_IMAGE'])

        # Prepare output FITS catalogue
        output = file_io.FITSCatalogue(
            f'{self._output_dir}/full_starcat-0000000.fits',
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True
        )

        # Collect columns
        # convert back to sigma for consistency
        data = {
            'X': x,
            'Y': y,
            'RA': ra,
            'DEC': dec,
            'EPS1': eps1,
            'EPS2': eps2,
            'CHI1': chi1,
            'CHI2': chi2,
            'SIZE': size,
            'FLAGS': flags,
            'FLAGS_EXT': flags_ext,
            'MAG': mag,
            'SNR': snr,
            'CCD_NB': ccd_nb,
        }

        # Write file
        output.save_as_fits(
            data,
            overwrite=True,
            sex_cat_path=self._input_file_list[0][0],
        )
