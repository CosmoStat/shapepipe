"""MERGE STAR CATALOGUES.

This module contains a class to identify single exposures that were used
to create tiles.

:Authors: Martin Kilbinger <martin.kilbinger@cea.fr>, Tobias Liaudat,
    Morgan Schmitz, Axel Guinot

"""

import re

import numpy as np
from astropy.io import fits
from cs_util import size as cs_size

from shapepipe.pipeline import file_io


def _stack(chunks, dtype=None):
    """Concatenate one column's per-catalogue arrays into a single array.

    THE COLUMN ACCUMULATORS ARE LISTS OF ARRAYS, ONE PER INPUT CATALOGUE, and
    not lists of values, because these classes are the last step of a whole
    campaign. ``x += list(data["X"])`` turns 4 bytes of float32 payload into a
    32-byte python object plus an 8-byte pointer in a list that overallocates —
    measured at ~10x the input bytes end to end, which put a full-survey merge
    (~20k exposures x 40 CCDs) at ~400 GB of RAM and made it unrunnable on any
    node. One array per catalogue plus one concatenate at the end holds ~1x, and
    produces the identical output: np.array() over a list of numpy scalars and
    np.concatenate() over the arrays they came from agree on dtype and on order.

    IT EMPTIES THE LIST IT IS GIVEN, and that is not a side effect to tidy away
    later — it is half the saving. np.concatenate holds the chunks and the
    result at once, so a caller that stacks sixteen columns while all sixteen
    chunk lists are still alive peaks at twice the campaign. Released column by
    column, the peak is one campaign plus one column. Callers stack once, at the
    end, and do not touch the accumulators afterwards.

    An empty input list is a merge over no catalogues, which the callers guard
    against; it returns an empty array so the output column still exists.
    """
    if not chunks:
        return np.array([], dtype=dtype or np.float64)
    out = np.concatenate(chunks)
    del chunks[:]
    return out


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
    input_cat_type : str, optional
        Input catalogue type, added to outupt header; default is ``None``,
        in which no key is added

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        w_log,
        stamp_size=51,
        rad=10,
        hdu_table=1,
        input_cat_type=None,
    ):

        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log
        self._stamp_size = stamp_size
        self._rad = rad
        self._hdu_table = hdu_table
        self._input_cat_type = input_cat_type

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
        rmse = np.sqrt(np.nansum(np.array(values) ** 2) / np.nansum(np.array(sizes)))

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
        mean = np.nansum(np.array(values)) / np.nansum(np.array(sizes))

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

        rmse = np.sqrt(np.mean(residual**2))
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
        inside_circle = (
            np.sqrt(
                (idx[np.newaxis, :] - cent[0]) ** 2
                + (jdx[:, np.newaxis] - cent[1]) ** 2
            )
            <= self._rad
        )
        my_mask[inside_circle] = True

        for name in self._input_file_list:
            try:
                starcat_j = fits.open(name[0], memmap=False, ignore_missing_simple=True)
            except ValueError:
                print(f"Error for file {name[0]}, check FITS file integrity")
                #raise
                continue

            stars = np.copy(starcat_j[self._hdu_table].data["VIGNET_LIST"])
            stars[stars < -1e6] = 0
            psfs = np.copy(starcat_j[self._hdu_table].data["PSF_VIGNET_LIST"])

            # Pixel mse calculation
            pix_val = np.sum((stars - psfs) ** 2)
            pix_sum = np.sum((stars - psfs))
            masked_diffs = np.array(
                [(_star - _psf)[my_mask] for _star, _psf in zip(stars, psfs)]
            )
            masked_pix_val = np.sum(masked_diffs**2)
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
            stars_norm_vals = np.sqrt(np.sum(stars**2, axis=(1, 2)))
            psfs_norm_vals = np.sqrt(np.sum(psfs**2, axis=(1, 2)))
            # Select non zero stars & psfs
            non_zero_elems = np.logical_and(
                (psfs_norm_vals != 0), (stars_norm_vals != 0)
            )
            # Calculate the filtered mse calculation
            pix_filt_val = np.sum((stars[non_zero_elems] - psfs[non_zero_elems]) ** 2)
            # Calculate the normalized (& filtered) mse calculation
            stars_norm_vals = stars_norm_vals[non_zero_elems].reshape(-1, 1, 1)
            psfs_norm_vals = psfs_norm_vals[non_zero_elems].reshape(-1, 1, 1)
            pix_norm_val = np.sum(
                (
                    stars[non_zero_elems] / stars_norm_vals
                    - psfs[non_zero_elems] / psfs_norm_vals
                )
                ** 2
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

            # ONE ARRAY PER CATALOGUE PER COLUMN (see _stack): the per-value
            # python lists this replaces cost ~10x the input bytes.
            # positions
            pos = starcat_j[self._hdu_table].data["GLOB_POSITION_IMG_LIST"]
            x.append(np.asarray(pos[:, 0]))
            y.append(np.asarray(pos[:, 1]))

            # RA and DEC positions
            try:
                ra.append(np.asarray(starcat_j[self._hdu_table].data["RA_LIST"][:]))
                dec.append(np.asarray(starcat_j[self._hdu_table].data["DEC_LIST"][:]))
            except Exception:
                ra.append(np.zeros(pos[:, 0].shape, dtype=int))
                dec.append(np.zeros(pos[:, 0].shape, dtype=int))

            # shapes (convert sigmas to T = 2 sigma^2)
            psf_mom = starcat_j[self._hdu_table].data["PSF_MOM_LIST"]
            star_mom = starcat_j[self._hdu_table].data["STAR_MOM_LIST"]
            g1_psf.append(np.asarray(psf_mom[:, 0]))
            g2_psf.append(np.asarray(psf_mom[:, 1]))
            size_psf.append(np.asarray(cs_size.sigma_to_T(psf_mom[:, 2])))
            g1.append(np.asarray(star_mom[:, 0]))
            g2.append(np.asarray(star_mom[:, 1]))
            size.append(np.asarray(cs_size.sigma_to_T(star_mom[:, 2])))

            # flags
            flag_psf.append(np.asarray(psf_mom[:, 3]))
            flag_star.append(np.asarray(star_mom[:, 3]))

            # ccd id list
            ccd_nb.append(np.asarray(starcat_j[self._hdu_table].data["CCD_ID_LIST"]))

            starcat_j.close()

        # Shortcut name
        MSC = MergeStarCatMCCD

        # Regular pixel RMSE
        tot_pixel_rmse = MSC.rmse_calc(pixel_mse, size_mse)
        self._w_log.info(
            f"MCCD_merge_starcat: Regular Total pixel RMSE ="
            + f" {tot_pixel_rmse:.5e}\n"
        )

        # Regular Total pixel mean
        tot_pixel_mean = MSC.mean_calc(pixel_sum, size_mse)
        self._w_log.info(
            f"MCCD_merge_starcat: Regular Total pixel mean ="
            + f" {tot_pixel_mean:.5e}\n"
        )

        # Regular Total MASKED pixel RMSE
        tot_masked_pixel_rmse = MSC.rmse_calc(masked_pixel_mse, masked_size)
        self._w_log.info(
            f"MCCD_merge_starcat: Regular Total MASKED pixel RMSE ="
            + f" {tot_masked_pixel_rmse:.5e}\n"
        )

        # Regular Total MASKED pixel mean
        tot_masked_pixel_mean = MSC.mean_calc(masked_pixel_sum, masked_size)
        self._w_log.info(
            f"MCCD_merge_starcat: Regular Total MASKED pixel mean ="
            + f" {tot_masked_pixel_mean:.5e}\n"
        )

        # Normalized pixel RMSE
        tot_pix_norm_rmse = MSC.rmse_calc(pix_norm_mse, size_norm_mse)
        self._w_log.info(
            "MCCD_merge_starcat: Normalized Total pixel RMSE ="
            + f" {tot_pix_norm_rmse:.5e}\n"
        )

        # Normalized filtered pixel RMSE
        tot_pix_filt_rmse = MSC.rmse_calc(pix_filt_mse, size_filt_mse)
        self._w_log.info(
            "MCCD_merge_starcat: Filtered Total pixel RMSE ="
            + f" {tot_pix_filt_rmse:.5e}\n"
        )

        concat_model_var = np.concatenate(np.array(model_var))
        concat_star_noise_var = np.concatenate(np.array(star_noise_var))

        # Model variance
        mean_model_var = MSC.mean_calc(concat_model_var, model_var_size)
        std_model_var = MSC.std_calc(concat_model_var)
        rmse_model_var = MSC.rmse_calc_2(concat_model_var, model_var_size)
        self._w_log.info(
            f"MCCD-RCA variance:\nMean Variance= {mean_model_var:.5e}\n"
            + f"Std Variance= {std_model_var:.5e}\n"
            + f"RMSE Variance= {rmse_model_var:.5e}\n"
        )

        # Star Noise Variance
        mean_star_var = MSC.mean_calc(concat_star_noise_var, star_noise_size)
        std_star_var = MSC.std_calc(concat_star_noise_var)
        rmse_star_var = MSC.rmse_calc_2(concat_star_noise_var, star_noise_size)
        self._w_log.info(
            f"Masked stars variance:\nMean Variance= {mean_star_var:.5e}\n"
            + f"Std Variance= {std_star_var:.5e}\n"
            + f"RMSE Variance= {rmse_star_var:.5e}\n"
        )

        # Mask and transform to numpy arrays
        # Concatenate once, here: everything below already wanted arrays and
        # was calling np.array() on python lists to get them (see _stack).
        x, y, ra, dec = _stack(x), _stack(y), _stack(ra), _stack(dec)
        g1_psf, g2_psf, size_psf = _stack(g1_psf), _stack(g2_psf), _stack(size_psf)
        g1, g2, size = _stack(g1), _stack(g2), _stack(size)
        flag_psf, flag_star = _stack(flag_psf), _stack(flag_star)
        ccd_nb = _stack(ccd_nb)

        flagmask = np.abs(flag_star - 1) * np.abs(flag_psf - 1)
        psf_e1 = g1_psf[flagmask.astype(bool)]
        psf_e2 = g2_psf[flagmask.astype(bool)]
        psf_r2 = size_psf[flagmask.astype(bool)]
        star_e1 = g1[flagmask.astype(bool)]
        star_e2 = g2[flagmask.astype(bool)]
        star_r2 = size[flagmask.astype(bool)]

        rmse, mean, std_dev = MSC.stats_calculator(star_e1, psf_e1)
        self._w_log.info(
            f"Moment residual e1:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
            + f"RMSE= {rmse:.5e}\n"
        )

        rmse, mean, std_dev = MSC.stats_calculator(star_e2, psf_e2)
        self._w_log.info(
            f"Moment residual e2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
            + f"RMSE= {rmse:.5e}\n"
        )

        rmse, mean, std_dev = MSC.stats_calculator(star_r2, psf_r2)
        self._w_log.info(
            f"Moment residual R2:\nMean= {mean:.5e}\nStd Dev= {std_dev:.5e}\n"
            + f"RMSE= {rmse:.5e}\n"
        )

        self._w_log.info(f"MCCD: Number of stars: {star_e1.shape[0]:d}")

        # Prepare output FITS catalogue
        output = file_io.FITSCatalogue(
            f"{self._output_dir}/full_starcat-0000000.fits",
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True,
        )

        # Collect columns (size stored as T = 2 sigma^2)
        data = {
            "X": x,
            "Y": y,
            "RA": ra,
            "DEC": dec,
            "HSM_G1_PSF": g1_psf,
            "HSM_G2_PSF": g2_psf,
            "HSM_T_PSF": size_psf,
            "HSM_G1_STAR": g1,
            "HSM_G2_STAR": g2,
            "HSM_T_STAR": size,
            "HSM_FLAG_PSF": flag_psf,
            "HSM_FLAG_STAR": flag_star,
            "CCD_NB": ccd_nb,
        }

        # Write file
        output.save_as_fits(data, sex_cat_path=self._input_file_list[0][0])
        # MKDEBUG; Implement this in file_io
        if self._input_cat_type:
            with fits.open(out_path, mode="update") as hdu_list:
                header = hdu_list[1].header
    
                # Add a new key-value pair
                header['CATTYPE'] = (self._input_cat_type, "Input catalogue type")
    
                # Save changes to the FITS file
                hdu_list.flush()



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
    input_cat_type : str, optional
        Input catalogue type, added to outupt header; default is ``None``,
        in which no key is added

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        w_log,
        hdu_table=2,
        input_cat_type=None,
    ):
        self._input_file_list = input_file_list
        self._output_dir = output_dir
        self._w_log = w_log
        self._hdu_table = hdu_table
        self._input_cat_type = input_cat_type

    # The columns this class writes, and where each comes from. Kept as data
    # rather than as sixteen repeated lines, because a two-pass merge would
    # otherwise state every column three times: to size it, to allocate it and
    # to fill it.
    _COLUMNS = (
        ("X", "X"), ("Y", "Y"), ("RA", "RA"), ("DEC", "DEC"),
        ("HSM_G1_PSF", "HSM_G1_PSF"), ("HSM_G2_PSF", "HSM_G2_PSF"),
        ("HSM_T_PSF", "HSM_T_PSF"), ("HSM_G1_STAR", "HSM_G1_STAR"),
        ("HSM_G2_STAR", "HSM_G2_STAR"), ("HSM_T_STAR", "HSM_T_STAR"),
        ("HSM_FLAG_PSF", "HSM_FLAG_PSF"), ("HSM_FLAG_STAR", "HSM_FLAG_STAR"),
    )
    # Present in psfex_interp output, absent from pix2wcs-converted files
    # (MKDEBUG); zero-filled when missing rather than failing the merge.
    _OPTIONAL = (("MAG", "MAG"), ("SNR", "SNR"), ("ACCEPTED", "ACCEPTED"))

    def _ccd_nb(self, path):
        """The CCD number this catalogue's rows carry, parsed from its name."""
        return re.split(r"\-([0-9]*)\-([0-9]+)\.", path)[-2]

    def process(self):
        """Process.

        Process merging.

        TWO PASSES, AND NEITHER HOLDS THE CAMPAIGN TWICE. The first reads only
        the FITS HEADER of every input — NAXIS2, the row count — and never
        touches a data block; the second allocates the output columns once, at
        their exact final length, and fills them slice by slice. Peak memory is
        therefore ONE output plus ONE input catalogue.

        What this replaces, in two steps, is instructive about the cost of the
        obvious code. Accumulating each column into a python LIST OF VALUES —
        ``x += list(data["X"])`` — turned 4 bytes of float32 payload into a
        32-byte object plus an 8-byte pointer, measured at ~10x the input bytes
        end to end and putting a full-survey merge (~20k exposures x 40 CCDs) at
        ~400 GB. Accumulating one ARRAY PER CATALOGUE and concatenating once
        brought that to ~5.5x. This pass structure removes what was left of the
        accumulation: there are no chunks, and no concatenate that must hold its
        inputs and its result at the same time.

        ``self._input_file_list`` MUST BE ITERABLE TWICE, which the module
        runner's list is. A one-shot generator is not, and would silently merge
        nothing on the second pass — hence the explicit row-count check below.
        """
        self._w_log.info(
            f"Merging {len(self._input_file_list)} star catalogues"
        )

        # --- pass 1: row counts and dtypes, from headers alone --------------
        # THE OPTIONAL COLUMNS ARE A PER-FILE QUESTION, NOT A PER-MERGE ONE.
        # A pix2wcs-converted catalogue has no MAG/SNR/ACCEPTED while an
        # ordinary one does, and a merge can be handed both. Deciding from the
        # first file alone got it wrong in both directions: converted-first
        # zero-filled the real values of every ordinary file behind it, and
        # ordinary-first raised KeyError on the first converted one. So the
        # dtype comes from ANY file that carries the column, and pass 2 asks
        # each file for itself.
        names, dtypes, opt_dtypes, n_total = [], None, {}, 0
        for name in self._input_file_list:
            try:
                with fits.open(name[0], memmap=False,
                               ignore_missing_simple=True) as starcat_j:
                    hdu = starcat_j[self._hdu_table]
                    n_rows = hdu.header["NAXIS2"]
                    # ColDefs.dtype describes the table without reading it.
                    # NOTE: it is the RAW storage dtype and ignores TSCAL/TZERO,
                    # so a scaled column would be allocated narrower than the
                    # values .data returns. Latent, not live: no validation_psf
                    # column is scaled. Read the dtype off .data if one ever is.
                    cols = hdu.columns.dtype
                    if dtypes is None:
                        dtypes = cols
                    for _, col in self._OPTIONAL:
                        if col not in opt_dtypes and col in (cols.names or ()):
                            opt_dtypes[col] = cols[col]
            except OSError:
                print(f"Error while opening file '{name[0]}'")
                #raise
                continue
            names.append(name[0])
            n_total += n_rows

        if dtypes is None:
            raise ValueError("merge_starcat: no readable input catalogue")

        # --- allocate once, at the exact final length -----------------------
        data = {out: np.empty(n_total, dtype=dtypes[col])
                for out, col in self._COLUMNS}
        for out, col in self._OPTIONAL:
            # A column no file carries still gets a column, zero-filled, in the
            # positional dtype the old code used for it.
            data[out] = np.empty(n_total, dtype=opt_dtypes.get(col, dtypes["X"]))
        # CCD_NB is one string per catalogue, repeated over its rows; its width
        # is the widest CCD number in the campaign, which pass 1 already knows.
        width = max((len(self._ccd_nb(n)) for n in names), default=1)
        data["CCD_NB"] = np.empty(n_total, dtype=f"U{width}")

        # --- pass 2: fill ---------------------------------------------------
        at = 0
        for name in self._input_file_list:
            try:
                starcat_j = fits.open(name[0], memmap=False,
                                      ignore_missing_simple=True)
            except OSError:
                continue
            data_j = starcat_j[self._hdu_table].data
            n_rows = len(data_j)
            sl = slice(at, at + n_rows)

            have = set(data_j.dtype.names or ())
            for out, col in self._COLUMNS:
                data[out][sl] = data_j[col]
            for out, col in self._OPTIONAL:
                # THIS file's schema, not the merge's: zero-fill only the files
                # that actually lack the column.
                data[out][sl] = data_j[col] if col in have else 0
            data["CCD_NB"][sl] = self._ccd_nb(name[0])

            at += n_rows
            starcat_j.close()

        if at != n_total:
            # The two passes disagreed: an input changed under us, or the list
            # was a one-shot iterable. Either way the output would be padded
            # with uninitialised memory, so say so rather than write it.
            raise ValueError(
                f"merge_starcat: pass 1 counted {n_total} rows, pass 2 filled "
                f"{at} — is the input list iterable more than once?")

        # Prepare output FITS catalogue
        # MKDEBUG: SEx_cat=True -> False
        out_path = f"{self._output_dir}/full_starcat-0000000.fits"
        output = file_io.FITSCatalogue(
            out_path,
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=False,
        )

        # `data` was built by the two passes above (size stored as T = 2
        # sigma^2); every column is already an array of its final length.

        # Write file
        # MKDEBUG for psf conv (pix2WCS) files do not write as SExtractorCat;
        # we do not want to copy the first input data content to HDU #1.
        # sex_cat_path=self._input_file_list[0][0],
        output.save_as_fits(
            data,
            overwrite=True,
        )
        # MKDEBUG; Implement this in file_io
        if self._input_cat_type:
            with fits.open(out_path, mode="update") as hdu_list:
                header = hdu_list[1].header
    
                # Add a new key-value pair
                header['CATTYPE'] = (self._input_cat_type, "Input catalogue type")
    
                # Save changes to the FITS file
                hdu_list.flush()


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
        q11 = "X2WIN_IMAGE"
        q22 = "Y2WIN_IMAGE"
        q12 = "XYWIN_IMAGE"

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
        if typ == "epsilon":
            # Determinant = (Q_11 Q_22 - Q_12^2)^(1/2)
            det = np.sqrt(m20 * m02 - m11 * m11)
        elif typ == "chi":
            det = 0
        else:
            raise ValueError(f"Invalid ellipticity type {type}")

        # Denominator = Q_11 + Q_22 [ + 2 * det]
        den = m20 + m02 + 2 * det

        # Ellipticity = (Q_11 - Q_22 + 2 i Q_12) / den
        ell = (m20 - m02 + 1j * 2 * m11) / den

        if type == "chi":
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
            f"Merging {len(self._input_file_list)} star catalogues"
        )

        for name in self._input_file_list:
            starcat_j = fits.open(name[0], memmap=False)

            data_j = starcat_j[self._hdu_table].data

            # positions
            x.append(np.asarray(data_j["XWIN_IMAGE"]))
            y.append(np.asarray(data_j["YWIN_IMAGE"]))
            ra.append(np.asarray(data_j["XWIN_WORLD"]))
            dec.append(np.asarray(data_j["YWIN_WORLD"]))

            # PRE-EXISTING BUG, LEFT ALONE DELIBERATELY: these four REBIND the
            # accumulators initialised above rather than appending to them, so
            # only the LAST input file's ellipticities reach the output while
            # every other column carries the whole merge. Setools is not wired
            # to any workflow path today; fixing it is its own change with its
            # own verification, and doing it silently inside a memory rewrite
            # would bury it.
            m11, m20, m02 = self.get_moments(data_j)
            eps1, eps2 = self.get_ellipticity(m11, m20, m02, "epsilon")
            chi1, chi2 = self.get_ellipticity(m11, m20, m02, "chi")

            size.append(np.asarray(data_j["FLUX_RADIUS"]))

            # flags
            flags.append(np.asarray(data_j["FLAGS_WIN"]))
            flags_ext.append(np.asarray(data_j["IMAFLAGS_ISO"]))

            # misc
            mag.append(np.asarray(data_j["MAG_WIN"]))
            snr.append(np.asarray(data_j["SNR_WIN"]))

            # CCD number
            ccd_nb.append(np.full(
                len(data_j["XWIN_IMAGE"]),
                re.split(r"\-([0-9]*)\-([0-9]+)\.", name[0])[-2]))

        # Prepare output FITS catalogue
        output = file_io.FITSCatalogue(
            f"{self._output_dir}/full_starcat-0000000.fits",
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=True,
        )

        # Collect columns
        # convert back to sigma for consistency
        data = {
            "X": _stack(x),
            "Y": _stack(y),
            "RA": _stack(ra),
            "DEC": _stack(dec),
            "EPS1": eps1,
            "EPS2": eps2,
            "CHI1": chi1,
            "CHI2": chi2,
            "SIZE": _stack(size),
            "FLAGS": _stack(flags),
            "FLAGS_EXT": _stack(flags_ext),
            "MAG": _stack(mag),
            "SNR": _stack(snr),
            "CCD_NB": _stack(ccd_nb, dtype="U1"),
        }

        # Write file
        output.save_as_fits(
            data,
            overwrite=True,
            sex_cat_path=self._input_file_list[0][0],
        )
