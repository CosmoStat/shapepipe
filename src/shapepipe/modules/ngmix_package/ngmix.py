"""NGMIX.

This module contains a class for ngmix shape measurement.

:Authors: Lucie Baumont, Axel Guinot

"""

import os
import re
import ngmix
import galsim
import numpy as np
from astropy.io import fits
from modopt.math.stats import sigma_mad
from ngmix.observation import Observation, ObsList
from scipy.spatial import cKDTree
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io

# Neighbour treatments selectable with the BLEND_HANDLING option.
BLEND_HANDLINGS = ("noisefill", "uberseg")

# Gaussian half-light radius per unit sigma: r50 = sqrt(2 ln 2) * sigma
# = 1.17741 * sigma.  With the ngmix area parameter T = 2 sigma^2 this
# gives r50 = SIGMA_TO_R50 * sqrt(T / 2) = sqrt(ln 2 * T).
SIGMA_TO_R50 = np.sqrt(2.0 * np.log(2.0))

METACAL_TYPES = ('noshear', '1p', '1m', '2p', '2m')


def get_mcal_flags(res):
    """Get Metacal Flags.

    Bitwise OR of the per-type metacal fit flags, the v1 contract for the
    downstream NGMIX_MCAL_FLAGS column: nonzero whenever any metacal
    type's galaxy fit failed.

    Parameters
    ----------
    res : dict
        MetacalBootstrapper result dict with one entry per metacal type.

    Returns
    -------
    int
        OR of all per-type ``flags``.
    """
    return int(np.bitwise_or.reduce(
        [res.get(name, {}).get('flags', 0) for name in METACAL_TYPES]
    ))


def get_prior(pixel_scale, rng, T_range=None, F_range=None):
    """Build ngmix joint prior for a 6-parameter galaxy model.

    Parameters
    ----------
    pixel_scale : float
        Pixel scale in arcsec (sets centroid prior width).
    rng : numpy.random.RandomState
        Random state for all priors.
    T_range : tuple of float, optional
        (min, max) for flat size prior; default (-1, 1e3).
    F_range : tuple of float, optional
        (min, max) for flat flux prior; default (-100, 1e9).

    Returns
    -------
    ngmix.joint_prior.PriorSimpleSep
    """
    if T_range is None:
        T_range = [-1.0, 1.0e3]
    if F_range is None:
        F_range = [-100.0, 1.0e9]

    cen_prior = ngmix.priors.CenPrior(
        cen1=0.0, cen2=0.0,
        sigma1=pixel_scale, sigma2=pixel_scale,
        rng=rng,
    )
    g_prior = ngmix.priors.GPriorBA(sigma=0.4, rng=rng)
    T_prior = ngmix.priors.FlatPrior(minval=T_range[0], maxval=T_range[1], rng=rng)
    F_prior = ngmix.priors.FlatPrior(minval=F_range[0], maxval=F_range[1], rng=rng)

    return ngmix.joint_prior.PriorSimpleSep(
        cen_prior=cen_prior,
        g_prior=g_prior,
        T_prior=T_prior,
        F_prior=F_prior,
    )


class Tile_cat():
    """Tile_cat.

    catalog measured on a tile

    Parameters
    ----------
    cat_path

    """
    def __init__(
        self,
        cat_path
    ):
        self.cat_path = cat_path
        if cat_path:
            self.get_data(cat_path)

    def get_data(self, cat_path):
        tile_cat = file_io.FITSCatalogue(
            cat_path,
            SEx_catalogue=True,
        )
        tile_cat.open()
        data = tile_cat.get_data()
        cols = data.dtype.names

        self.obj_id = np.copy(data['NUMBER'])
        self.ra = np.copy(data['XWIN_WORLD'])
        self.dec = np.copy(data['YWIN_WORLD'])

        # Optional columns — may be absent in external (non-SExtractor) catalogs
        self.flux = np.copy(data['FLUX_AUTO']) if 'FLUX_AUTO' in cols else None
        self.vign = np.copy(data['VIGNET']) if 'VIGNET' in cols else None

        tile_cat.close()

class Postage_stamp():
    """Galaxy Postage Stamp.

    Class to hold catalog of postage stamps for a single galaxy

    Parameters
    ----------
    bkg_sub: bool

    megacam_flip: bool
    We probably want to put weight and flag options here too

    """
    def __init__(
        self,
        bkg_sub=True,
        megacam_flip=True

    ):
        self.gals = []
        self.psfs = []
        self.weights = []
        self.flags = []
        self.bkg_rms = []
        # Per-epoch segmentation stamps (object NUMBERs), used only by the
        # "uberseg" blend handling; empty for the default noise-fill path.
        self.segs = []
        self.jacobs = []
        # Per-epoch full WCS and the object's sky position, used only by the
        # "wcs" centroid source (skipped for the default "hsm" path).
        self.wcs = []
        self.ra = []
        self.dec = []
        self.bkg_sub = bkg_sub
        self.megacam_flip = megacam_flip

class Vignet():
    """Vignet.

    Class to hold catalog of postage stamps

    Parameters
    ----------
    gal_vignet_path
    bkg_vignet_path
    psf_vignet_path
    weight_vignet_path
    flag_vignet_path
    f_wcs_path
    bkg_rms_vignet_path
    """
    def __init__(
        self,
        gal_vignet_path,
        bkg_vignet_path,
        psf_vignet_path,
        weight_vignet_path,
        flag_vignet_path,
        f_wcs_path,
        bkg_rms_vignet_path=None,

    ):
        self.f_wcs_file = SqliteDict(f_wcs_path)
        self.gal_vign_cat = SqliteDict(gal_vignet_path)
        self.bkg_vign_cat = SqliteDict(bkg_vignet_path)
        self.psf_vign_cat = SqliteDict(psf_vignet_path)
        self.weight_vign_cat = SqliteDict(weight_vignet_path)
        self.flag_vign_cat = SqliteDict(flag_vignet_path)
        self.bkg_rms_vign_cat = (
            SqliteDict(bkg_rms_vignet_path)
            if bkg_rms_vignet_path is not None
            else None
        )

    def close(self):
        self.f_wcs_file.close()
        self.gal_vign_cat.close()
        self.bkg_vign_cat.close()
        self.flag_vign_cat.close()
        self.weight_vign_cat.close()
        self.psf_vign_cat.close()
        if self.bkg_rms_vign_cat is not None:
            self.bkg_rms_vign_cat.close()

class Ngmix(object):
    """Ngmix.

    Class to handle NGMIX shapepe measurement.

    Parameters
    ----------
    input_file_list : list
        Input files
    output_dir : str
        Output directory
    file_number_string : str
        File numbering scheme
    zero_point : float
        Photometric zero point
    pixel_scale : float
        Pixel scale in arcsec
    f_wcs_path : str
        Path to merged single-exposure single-HDU headers
    w_log : logging.Logger
        Logging instance
    save_batch : int, optional
        Save output catalogue in batches of this size; detaul is ``-1`` (no
        batch save)
    id_obj_min : int, optional
        First galaxy ID to process, not used if the value is set to ``-1``;
        the default is ``-1``
    id_obj_max : int, optional
        Last galaxy ID to process, not used if the value is set to ``-1``;
        the default is ``-1``
    centroid_source : {"hsm", "wcs"}, optional
        How to place the galaxy Jacobian origin for the centroid prior. The
        default ``"hsm"`` re-centers on the HSM adaptive-moment centroid
        (robust for galaxies); ``"wcs"`` uses the catalog sky position
        projected through the WCS (better for stars, whose HSM moments are
        noisy). See :func:`make_ngmix_observation`.

    Raises
    ------
    IndexError
        If the length of the input file list is incorrect

    """

    def __init__(
        self,
        input_file_list,
        output_dir,
        file_number_string,
        zero_point,
        pixel_scale,
        f_wcs_path,
        w_log,
        save_batch=-1,
        id_obj_min=-1,
        id_obj_max=-1,
        centroid_source="hsm",
        blend_handling="noisefill",
    ):

        if len(input_file_list) not in {6, 7}:
            raise IndexError(
                f"Input file list has length {len(input_file_list)},"
                + " required is 6 or 7"
            )

        if blend_handling not in BLEND_HANDLINGS:
            raise ValueError(
                f"Unknown BLEND_HANDLING '{blend_handling}'; expected one of"
                + f" {BLEND_HANDLINGS}"
            )

        self._tile_cat_path = input_file_list[0]
        self._vignet_cat = Vignet(
            input_file_list[1],
            input_file_list[2],
            input_file_list[3],
            input_file_list[4],
            input_file_list[5],
            f_wcs_path,
            input_file_list[6] if len(input_file_list) == 7 else None,
        )

        self._output_dir = output_dir
        self._file_number_string = file_number_string

        self._zero_point = zero_point
        self._pixel_scale = pixel_scale

        self._f_wcs_path = f_wcs_path

        self._save_batch = save_batch
        self._id_obj_min = id_obj_min
        self._id_obj_max = id_obj_max
        self._centroid_source = centroid_source
        self._blend_handling = blend_handling

        self._w_log = w_log

        # Initiatlise random generator
        seed = int(''.join(re.findall(r'\d+', self._file_number_string)))
        self._rng = np.random.RandomState(seed)
        self._w_log.info(f'Random generator initialisation seed = {seed}')

    @classmethod
    def MegaCamFlip(self, vign, ccd_nb):
        """Flip for MegaCam.

        MegaPipe has CCDs that are upside down. This function flips the
        postage stamps in these CCDs. TO DO: This will give incorrect results
        when used with THELI ccds.  Fix this.

        Parameters
        ----------
        vign : numpy.ndarray
            Array containing the postage stamp to flip
        ccd_nb : int
            ID of the CCD containing the postage stamp

        Returns
        -------
        numpy.ndarray
            The flipped postage stamp

        """
        if ccd_nb < 18 or ccd_nb in [36, 37]:
            # swap x axis so origin is on top-right
            return np.rot90(vign, k=2)
        else:
            # swap y axis so origin is on bottom-left
            return vign

    def get_prior(self, T_range=None, F_range=None):
        """Get Prior.

        Returns
        -------
        ngmix.joint_prior.PriorSimpleSep
        """
        return get_prior(
            self._pixel_scale, self._rng,
            T_range=T_range, F_range=F_range,
        )

    def compile_results(self, results):
        """Compile Results.

        Prepare the results of NGMIX before saving. TO DO: add snr_r and T_r
        This needs to be updated
        Parameters
        ----------
        results : dict
            Results of NGMIX metacal

        Returns
        -------
        dict
            Compiled results ready to be written to a file
            note: psfo is the original image psf from psfex or mccd

        Raises
        ------
        KeyError
            If SNR key not found

        """
        # Output HDU order. Same set as METACAL_TYPES, but kept in this
        # fixed order so output catalogues stay byte-reproducible; the check
        # below guards against the two lists silently diverging.
        names = ["1m", "1p", "2m", "2p", "noshear"]
        if set(names) != set(METACAL_TYPES):
            raise ValueError(
                "compile_results metacal type list is out of sync with"
                + " METACAL_TYPES"
            )
        names2 = [
            'id',
            'n_epoch_model',
            'mcal_types_fail',
            'nfev_fit',
            'g1_psfo_ngmix',
            'g2_psfo_ngmix',
            'g1_err_psfo_ngmix',
            'g2_err_psfo_ngmix',
            'g1',
            'g1_err',
            'g2',
            'g2_err',
            'T',
            'T_err',
            'Tpsf',
            'Tpsf_err',
            'r50',
            'r50_err',
            'r50psf',
            'r50psf_err',
            'g1_psf',
            'g2_psf',
            'flux',
            'flux_err',
            's2n',
            'mag',
            'mag_err',
            'flags',
            'mcal_flags'
        ]
        output_dict = {k: {kk: [] for kk in names2} for k in names}
        for idx in range(len(results)):
            for name in names:
                fit = results[idx][name]

                # ngmix 2.x does not raise on fit failure: after ntry the
                # result keeps flags != 0 and carries none of the
                # measurement keys (g, g_cov, T, T_err, flux, flux_err,
                # s2n).  NaN-fill those so failed types are recorded with
                # their flags instead of crashing the tile on a KeyError.
                flux = fit.get("flux", np.nan)
                flux_err = fit.get("flux_err", np.nan)
                g = np.asarray(fit.get("g", (np.nan, np.nan)))
                g_cov = np.asarray(
                    fit.get("g_cov", np.full((2, 2), np.nan))
                )
                T_gal = fit.get("T", np.nan)
                T_gal_err = fit.get("T_err", np.nan)

                mag = -2.5 * np.log10(flux) + self._zero_point
                mag_err = np.abs(-2.5 * flux_err / (flux * np.log(10)))

                output_dict[name]["id"].append(results[idx]["obj_id"])
                output_dict[name]["n_epoch_model"].append(
                    results[idx]["n_epoch_model"]
                )
                output_dict[name]["mcal_types_fail"].append(
                    results[idx]["mcal_types_fail"]
                )
                # ngmix 2.x reports the solver's function-evaluation count
                # (nfev, ~tens-hundreds; -1 on some failures), not the v1
                # 1-5 retry count, so the column is named accordingly.
                output_dict[name]["nfev_fit"].append(
                    fit.get("nfev", np.nan)
                )
                output_dict[name]["g1_psfo_ngmix"].append(
                    results[idx]["g_PSFo"][0]
                )
                output_dict[name]["g2_psfo_ngmix"].append(
                    results[idx]["g_PSFo"][1]
                )
                output_dict[name]["g1_err_psfo_ngmix"].append(
                    results[idx]["g_err_PSFo"][0]
                )
                output_dict[name]["g2_err_psfo_ngmix"].append(
                    results[idx]["g_err_PSFo"][1]
                )
                output_dict[name]["T"].append(T_gal)
                output_dict[name]["T_err"].append(T_gal_err)
                output_dict[name]["Tpsf"].append(results[idx]["T_PSFo"])
                output_dict[name]["Tpsf_err"].append(results[idx]["T_err_PSFo"])
                output_dict[name]["g1_psf"].append(results[idx]["g_PSFo"][0])
                output_dict[name]["g2_psf"].append(results[idx]["g_PSFo"][1])

                # Galaxy half-light radius from the fitted area T = 2 sigma^2:
                # r50 = sqrt(ln 2 * T), with d r50 / d T = r50 / (2 T)
                if T_gal > 0:
                    r50_gal = SIGMA_TO_R50 * np.sqrt(T_gal / 2)
                    r50_gal_err = r50_gal * T_gal_err / (2 * T_gal)
                else:
                    r50_gal = r50_gal_err = np.nan
                output_dict[name]['r50'].append(r50_gal)
                output_dict[name]['r50_err'].append(r50_gal_err)
                output_dict[name]['r50psf'].append(results[idx]["r50_PSFo"])
                output_dict[name]['r50psf_err'].append(
                    results[idx]["r50_err_PSFo"]
                )
                output_dict[name]["g1"].append(g[0])
                output_dict[name]["g2"].append(g[1])
                output_dict[name]["g1_err"].append(np.sqrt(g_cov[0, 0]))
                output_dict[name]["g2_err"].append(np.sqrt(g_cov[1, 1]))
                output_dict[name]["flux"].append(flux)
                output_dict[name]["flux_err"].append(flux_err)
                output_dict[name]["mag"].append(mag)
                output_dict[name]["mag_err"].append(mag_err)

                if "s2n" in fit:
                    output_dict[name]["s2n"].append(fit["s2n"])
                elif "s2n_r" in fit:
                    output_dict[name]["s2n"].append(fit["s2n_r"])
                elif fit["flags"] != 0:
                    output_dict[name]["s2n"].append(np.nan)
                else:
                    raise KeyError("No SNR key (s2n, s2n_r) found in results")

                output_dict[name]["flags"].append(fit["flags"])
                output_dict[name]["mcal_flags"].append(
                    results[idx].get("mcal_flags", 0)
                )

        return output_dict

    def get_output_path(self, directory):
        """Get Output Path.

        Return path of output ngmix catalogue file.

        Parameters
        ----------
        directoy: str
            directory name

        Returns
        -------
        str
            output path

        """
        return f"{directory}/ngmix{self._file_number_string}.fits"


    def save_results(self, output_dict):
        """Save Results.

        Save the results into a FITS file.

        Parameters
        ----------
        output_dict: dict
            Dictionary containing the results

        """
        n_hdu = len(output_dict.keys())
        if n_hdu != 5:
            raise IndexError(
                f"FITS output file data has {n_hdu} HDUs,"
                + " expected are 5"
            )

        output_name = self.get_output_path(self._output_dir)
        if not os.path.exists(output_name):
            f_out = file_io.FITSCatalogue(
                output_name, open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
            )
            for key in output_dict.keys():
                f_out.save_as_fits(output_dict[key], ext_name=key.upper())
            return

        with fits.open(output_name, mode='update') as hdul:

            # Iterate through HDUs (assuming they are Binary Table HDUs with data)
            for idx, hdu in enumerate(hdul):
                if isinstance(hdu, fits.BinTableHDU):  # Check for table data

                    # HDU extension name
                    ext_name = hdu.name.lower()
                    if ext_name not in output_dict:
                        raise ValueError(
                            f"HDU extension {ext_name} from existing FITS"
                            + " file not found in data"
                        )

                    # Existing data
                    existing_data = hdu.data
                    existing_dtype = hdu.data.dtype

                    # New data
                    new_data = output_dict[ext_name]

                    # Verify that all column names in existing_data exist
                    # in new_data_dict
                    if not all(
                        colname in new_data
                        for colname in existing_dtype.names
                    ):
                        raise ValueError(
                            "Mismatch between existing columns"
                            + f" ({existing_dtype.names}) and new data"
                            + f" columns ({list(new_data)})."
                        )

                    # New data to be appended
                    structured_data = np.zeros(
                        len(next(iter(new_data.values()))),
                        dtype=existing_data.dtype,
                    )
                    for colname in existing_data.dtype.names:
                        structured_data[colname] = new_data[colname]

                    # Combine existing and new data
                    updated_data = np.append(existing_data, structured_data)

                    # Update the data in the HDU
                    hdu.data = updated_data

            # Save changes to the FITS file
            hdul.flush()

    @classmethod
    def check_key(self, expccd_name_tmp, vign_cat, vignet_path):
        if expccd_name_tmp not in vign_cat:
            raise KeyError(
                f"Key '{expccd_name_tmp}' (exposure CCD ID from PSF postage stamp list)"
                + " not found in postage stamp database"
                + f" file '{vignet_path}'"
            )

    def process(self):
        """Process.

        Funcion to processs NGMIX.
        organizes object cutouts from detection catalog in image, 
        weight, and flag files
        per object: 
            gathers wcs and psf info from exposures
            background subtracts (make this an option)
            scales by relative zeropoints
            runs metacal convolutions and ngmix fitting
        Returns
        -------
        dict
            Dictionary containing the NGMIX metacal results

        """
        tile_cat = Tile_cat(self._tile_cat_path)
        vignet_cat = self._vignet_cat

        final_res = []
        prior = self.get_prior()

        count = 0
        n_empty_cat = 0
        n_no_epoch = 0
        n_ngmix_fail = 0
        n_fitted = 0
        id_first = -1
        id_last = -1
        count_batch = 0
        saved_batch_cumul = 0

        for i_tile, obj_id in enumerate(tile_cat.obj_id):
            if self._id_obj_min > 0 and obj_id < self._id_obj_min:
                continue
            if self._id_obj_max > 0 and obj_id > self._id_obj_max:
                continue
            if id_first == -1:
                id_first = obj_id
            id_last = obj_id
            count += 1

            # Skip objects with no multi-epoch PSF or vignet data
            if (vignet_cat.psf_vign_cat[str(obj_id)] == 'empty'
                    or vignet_cat.gal_vign_cat[str(obj_id)] == 'empty'):
                n_empty_cat += 1
                continue

            stamp = prepare_postage_stamps(vignet_cat, obj_id, i_tile, tile_cat)

            if len(stamp.gals) == 0:
                n_no_epoch += 1
                continue

            try:
                flux_guess = (
                    tile_cat.flux[i_tile]
                    if tile_cat.flux is not None
                    else 1.0
                )
                res, psf_res = do_ngmix_metacal(
                    stamp,
                    prior,
                    flux_guess,
                    self._rng,
                    centroid_source=self._centroid_source,
                    blend_handling=self._blend_handling,
                    object_number=obj_id,
                )
            except Exception as ee:
                self._w_log.info(
                    f'ngmix failed for object ID={obj_id}.\nMessage: {ee}'
                )
                n_ngmix_fail += 1
                continue

            res['obj_id'] = obj_id
            # epochs that survived the PSF fit and entered the model,
            # not the number of epochs submitted (v1 contract)
            res['n_epoch_model'] = psf_res['n_epoch']
            # Count of metacal fit types (0-5) with nonzero fit flags.
            # (In ngmix v1 the same-named column counted moments-initial-guess
            # failures from get_guess, which no longer exists — hence the
            # rename to mcal_types_fail / NGMIX_MCAL_TYPES_FAIL.)
            res['mcal_types_fail'] = sum(
                1 for k in METACAL_TYPES
                if res.get(k, {}).get('flags', 0) != 0
            )
            res['mcal_flags'] = get_mcal_flags(res)
            # PSF half-light radius r50 = sqrt(2 ln 2) * sigma with
            # sigma = sqrt(T / 2); error from d sigma / d T = 1 / (4 sigma)
            sigma_psfo = np.sqrt(max(psf_res['T_psf'], 0) / 2)
            res['g_PSFo'] = psf_res['g_psf']
            res['g_err_PSFo'] = psf_res['g_psf_err']
            res['T_PSFo'] = psf_res['T_psf']
            res['T_err_PSFo'] = psf_res['T_psf_err']
            res['r50_PSFo'] = SIGMA_TO_R50 * sigma_psfo
            res['r50_err_PSFo'] = (
                SIGMA_TO_R50 * psf_res['T_psf_err'] / (4 * sigma_psfo)
                if sigma_psfo > 0 else np.nan
            )
            final_res.append(res)
            n_fitted += 1
            count_batch += 1

            if self._save_batch > 0 and count_batch == self._save_batch:
                res_dict = self.compile_results(final_res)
                self.save_results(res_dict)
                saved_batch_cumul += count_batch
                self._w_log.info(
                    f"Batch-saved {count_batch} ({len(res_dict)} valid) objects,"
                    + f" cumul={saved_batch_cumul}"
                )
                final_res = []
                count_batch = 0

            if count % 1000 == 0:
                self._w_log.info(
                    f"Progress: {count} iterated, {n_empty_cat} empty catalog,"
                    + f" {n_no_epoch} no valid epoch, {n_ngmix_fail} fit failed,"
                    + f" {n_fitted} fitted"
                )

        self._w_log.info(
            f"ngmix loop finished: {count} iterated"
            + f" (id {id_first}/{id_last}),"
            + f" {n_empty_cat} empty catalog,"
            + f" {n_no_epoch} no valid epoch,"
            + f" {n_ngmix_fail} fit failed,"
            + f" {n_fitted} fitted"
        )

        vignet_cat.close()

        # Put all results together
        res_dict = self.compile_results(final_res)

        # Save results
        self.save_results(res_dict)

def prepare_postage_stamps(vignet, obj_id, i_tile, tile_cat):
    # define per-object lists of individual exposures to go into ngmix
    stamp = Postage_stamp()
    #identify exposure and ccd number from psf catalog
    psf_expccd_names = list(vignet.psf_vign_cat[str(obj_id)].keys())
    for expccd_name in psf_expccd_names:
        exp_name, ccd_n = re.split('-', expccd_name)

        gal_vign = (
            vignet.gal_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        )

        if np.all(gal_vign == 0):
            continue
        
        if stamp.bkg_sub:
            bkg_vign = (
                vignet.bkg_vign_cat[str(obj_id)][expccd_name]['VIGNET']
            )
            gal_vign_sub_bkg = background_subtract(
                gal_vign,
                bkg_vign
            )
        else:
            gal_vign_sub_bkg = gal_vign

        tile_vign = (
            np.copy(tile_cat.vign[i_tile])
            if tile_cat.vign is not None
            else None
        )
        if stamp.megacam_flip and tile_vign is not None:
            tile_vign = Ngmix.MegaCamFlip(tile_vign, int(ccd_n))

        flag_vign = (
            vignet.flag_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        )
        if tile_vign is not None:
            flag_vign[np.where(tile_vign == -1e30)] = 2**10
        v_flag_tmp = flag_vign.ravel()
        # remove objects that are more than 1/3 masked
        if len(np.where(v_flag_tmp != 0)[0]) / v_flag_tmp.size > 1 / 3.0:
            continue

        weight_vign = vignet.weight_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        bkg_rms_vign = (
            vignet.bkg_rms_vign_cat[str(obj_id)][expccd_name]['VIGNET']
            if vignet.bkg_rms_vign_cat is not None
            else None
        )

        epoch_wcs = vignet.f_wcs_file[exp_name][int(ccd_n)]['WCS']
        jacob = get_galsim_jacobian(
            epoch_wcs,
            tile_cat.ra[i_tile],
            tile_cat.dec[i_tile]
        )

        header = fits.Header.fromstring(
            vignet.f_wcs_file[exp_name][int(ccd_n)]['header']
        )

        # rescale by relative zero-points
        (
            gal_vign_scaled,
            weight_vign_scaled,
            bkg_rms_vign_scaled,
        ) = rescale_epoch_fluxes(
            gal_vign_sub_bkg,
            weight_vign,
            header,
            bkg_rms_vign,
        )

        # gather postage stamps in all of the epochs
        stamp.gals.append(gal_vign_scaled)
        stamp.psfs.append(
            vignet.psf_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        )
        stamp.weights.append(weight_vign_scaled)
        stamp.flags.append(flag_vign)
        stamp.bkg_rms.append(bkg_rms_vign_scaled)
        stamp.jacobs.append(jacob)
        # For the "wcs" centroid source (see make_ngmix_observation).
        stamp.wcs.append(epoch_wcs)
        stamp.ra.append(tile_cat.ra[i_tile])
        stamp.dec.append(tile_cat.dec[i_tile])

    return stamp

def background_subtract(gal,bkg):
    """background subtraction.
        
    Parameters
    ----------
    gal : numpy.ndarray
        galaxy image
    bkg : numpy.ndarray
        background
        
    Returns
    -------
    numpy.ndarray
        background subtracted galaxy
    """

    # background subtraction
    gal_vign_sub_bkg = gal - bkg

    return gal_vign_sub_bkg

def rescale_epoch_fluxes(gal, weight, header, bkg_rms=None):
    """rescale epochs by relative zeropoints to be on the same flux scale
        
    Parameters
    ----------
    gal : numpy.ndarray
        background subtracted galaxy image
    weight : numpy.ndarray
        weight image
    header : 
        image header
    bkg_rms : numpy.ndarray, optional
        Background RMS image
        
    Returns
    -------
    numpy.ndarray
        rescaled galaxy image
    numpy.ndarray
        rescaled weight image
    numpy.ndarray or None
        rescaled background RMS image
    """
    Fscale = header['FSCALE']

    gal_scaled = gal * Fscale
    weight_scaled = weight * 1 / Fscale ** 2
    bkg_rms_scaled = bkg_rms * Fscale if bkg_rms is not None else None

    return gal_scaled, weight_scaled, bkg_rms_scaled

def get_galsim_jacobian(wcs, ra, dec):
    """Get local wcs.
    This produces a galsim jacobian at a point.  We call it local_wcs because we convert to a ngmix object to create the jacobian later.
    TO DO: can we do this within ngmix?

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        WCS object for which we want the Jacobian
    ra : float
        RA position of the center of the vignet (in degrees)
    dec : float
        Dec position of the center of the vignet (in degress)

    Returns
    -------
    galsim.wcs.BaseWCS.jacobian
        Jacobian of the WCS at the required position

    """
    g_wcs = galsim.fitswcs.AstropyWCS(wcs=wcs)
    world_pos = galsim.CelestialCoord(
        ra=ra * galsim.angle.degrees,
        dec=dec * galsim.angle.degrees,
    )
    galsim_jacob = g_wcs.jacobian(world_pos=world_pos)

    return galsim_jacob


def get_noise(gal, weight, guess, pixel_scale, thresh=1.2):
    """Get Noise.
    TO DO: modify guess, pixel scale
    Compute the sigma of the noise from an object postage stamp.
    Use a guess on the object size, ellipticity and flux to create a window
    function.

    Parameters
    ----------
    gal : numpy.ndarray
        Galaxy image
    weight : numpy.ndarray
        Weight image
    guess : list
        Gaussian parameters fot the window function
        ``[x0, y0, g1, g2, T, flux]``
    pixel_scale : float
        Pixel scale of the galaxy image
    thresh : float, optional
        Threshold to cut the window function,
        cut = ``thresh`` * sigma_noise;  the default is ``1.2``

    Returns
    -------
    float
        Sigma of the noise on the galaxy image

    """
    img_shape = gal.shape
    m_weight = weight != 0

    sig_tmp = sigma_mad(gal[m_weight])

    gauss_win = galsim.Gaussian(sigma=np.sqrt(guess[4] / 2), flux=guess[5])
    gauss_win = gauss_win.shear(g1=guess[2], g2=guess[3])
    gauss_win = gauss_win.drawImage(
        nx=img_shape[0], ny=img_shape[1], scale=pixel_scale
    ).array

    m_weight = weight[gauss_win < thresh * sig_tmp] != 0

    sig_noise = sigma_mad(gal[gauss_win < thresh * sig_tmp][m_weight])

    return sig_noise

def uberseg_weight(weight, seg, object_number):
    """Zero a stamp's weight on neighbour-side pixels — the UberSeg treatment.

    UberSeg is Erin Sheldon's MEDS/ngmix neighbour mask (``esheldon/meds``,
    https://github.com/esheldon/meds — ``MEDS.get_uberseg`` /
    ``meds._uberseg.uberseg_tree``, in MEDS itself "adapted from Niall MacCrann
    and Joe Zuntz", and used in the DES shear pipeline). Each stamp pixel is
    assigned to the object whose segmentation footprint it lies nearest to — a
    nearest-segment Voronoi partition — and only pixels assigned to the central
    object keep their weight.

    This function **reimplements** the partition rather than depending on
    ``meds``: it is a five-line ``scipy.spatial.cKDTree`` nearest-neighbour
    query, and pulling in the whole MEDS-file library — plus ``esutil`` and a
    C-extension build (``fitsio`` is already in the ShapePipe stack; ``esutil``
    is not) — for one standard geometric operation is disproportionate. The
    reimplementation was validated bit-for-bit against Sheldon's own reference
    ``get_uberseg`` over a battery of synthetic and random multi-object seg
    maps: the two masks agree on every pixel except exact Voronoi-boundary ties
    (a measure-zero tie-break convention — ``argmin`` first-min vs cKDTree
    order — scientifically inert). See the ``uberseg`` fiber's ``validation/``
    harness for the equivalence proof and figures. The cKDTree stands in for
    the C k-d tree; results are identical up to that boundary convention.

    Because the partition is by distance to the nearest footprint, the pixels
    surviving around a compact central object form a single connected,
    roughly circular core; the "circularisation" is emergent geometry, not a
    separate aperture. Unlike the noise-fill treatment, masked pixels are
    handed to ngmix as a hard mask (weight = 0), never replaced by a noise
    realisation.

    Parameters
    ----------
    weight : numpy.ndarray
        Per-pixel weight (inverse variance) map for the stamp.
    seg : numpy.ndarray
        Segmentation map on the same grid as ``weight``: 0 for sky, the
        SExtractor object number for each detected object's footprint.
    object_number : int
        Segmentation label of the central object (its SExtractor ``NUMBER``).

    Returns
    -------
    numpy.ndarray
        Copy of ``weight`` with neighbour-side pixels zeroed.
    """
    weight = np.copy(weight)

    obj_pix = np.argwhere(seg != 0)
    # No detected footprint, or only sky plus the central object: there is
    # no neighbour to mask, so the weight is unchanged (cf. MEDS' early
    # ``len(np.unique(seg)) == 2`` return).
    if obj_pix.shape[0] == 0:
        return weight
    labels = seg[seg != 0]
    if np.all(labels == object_number):
        return weight

    # Nearest segmentation pixel for every stamp pixel; zero the weight
    # wherever that nearest footprint belongs to a neighbour rather than to
    # the central object.
    grid = np.indices(seg.shape).reshape(2, -1).T
    _, nearest = cKDTree(obj_pix).query(grid)
    nearest_label = labels[nearest].reshape(seg.shape)
    weight[nearest_label != object_number] = 0.0

    return weight

def prepare_ngmix_weights(
    gal, weight, flag, rng, bkg_rms=None,
    blend_handling="noisefill", seg=None, object_number=None,
):
    """bookkeeping for ngmix weights. runs on a single galaxy and epoch
        pixel scale and galaxy guess
        TO DO: decide if we want galaxy guess stuff

    Parameters
    ----------
    gal : numpy.ndarray
    weight : numpy.ndarray
    flag : numpy.ndarray
    rng : numpy.random.RandomState
        Random state for the noise realisations (seeded per tile for
        reproducibility).
    bkg_rms : numpy.ndarray, optional
        Per-pixel background RMS map. If supplied, unmasked pixels use
        ``1 / bkg_rms**2`` as the ngmix inverse variance.
    blend_handling : {"noisefill", "uberseg"}, optional
        How to treat pixels shared with a neighbour. ``"noisefill"`` (default)
        replaces flagged pixels with a noise realisation and keeps their
        inverse-variance weight — the historical behaviour. ``"uberseg"``
        instead hard-masks (weight = 0) every pixel closer to a neighbour's
        segmentation footprint than to the central object's, leaving the
        image untouched (see :func:`uberseg_weight`).
    seg : numpy.ndarray, optional
        Segmentation map on the stamp grid (object NUMBERs). Required for
        ``blend_handling="uberseg"``; ignored otherwise.
    object_number : int, optional
        Central object's segmentation label. Required for
        ``blend_handling="uberseg"``; ignored otherwise.

    Returns
    -------
    numpy.ndarray
        Galaxy image. For ``"noisefill"`` masked pixels are replaced by noise;
        for ``"uberseg"`` the image is returned untouched.
    numpy.ndarray
        Variance map for NGMIX.
    numpy.ndarray
        Noise image.
    """
    if blend_handling not in BLEND_HANDLINGS:
        raise ValueError(
            f"Unknown blend_handling '{blend_handling}'; expected one of"
            + f" {BLEND_HANDLINGS}"
        )

    mask = np.copy(weight) != 0
    mask[flag != 0] = False

    if bkg_rms is None:
        sig_noise = sigma_mad(gal)
        # Guard the degenerate constant stamp (sigma_mad == 0): 0 * inf
        # would otherwise put NaN in a fully-masked weight map.
        weight_map = (
            mask.astype(float) / sig_noise ** 2
            if sig_noise > 0
            else np.zeros_like(gal, dtype=float)
        )
    else:
        valid_rms = np.isfinite(bkg_rms) & (bkg_rms > 0)
        mask &= valid_rms
        weight_map = np.zeros_like(gal, dtype=float)
        weight_map[mask] = 1.0 / bkg_rms[mask] ** 2
        # Per-pixel noise sigma for the realisations below: metacal's
        # fixnoise bookkeeping (1/w + 1/w_noise) assumes the noise image
        # is a faithful realisation of the per-pixel variance the weights
        # claim; a scalar sigma there mis-reports errors and erodes the
        # inverse-variance advantage whenever the RMS map actually varies.
        sig_noise = (
            np.where(valid_rms, bkg_rms, np.median(bkg_rms[mask]))
            if mask.any()
            else sigma_mad(gal)
        )

    noise_img = rng.standard_normal(gal.shape) * sig_noise
    noise_img_gal = rng.standard_normal(gal.shape) * sig_noise

    gal_masked = np.copy(gal)
    if blend_handling == "uberseg":
        # Hard-mask neighbour-side pixels (weight -> 0) from the segmentation
        # geometry; the image is left untouched (the masked pixels carry no
        # weight, so ngmix ignores them in the likelihood). Bad/flagged
        # pixels already sit at weight 0 from the mask above.
        if seg is None or object_number is None:
            raise ValueError(
                "blend_handling='uberseg' requires a segmentation map and the"
                + " central object_number; none reached prepare_ngmix_weights."
                + " The seg-map source is gated (see CosmoStat/shapepipe issue"
                + " on UberSeg segmentation plumbing)."
            )
        weight_map = uberseg_weight(weight_map, seg, object_number)
    elif (~mask).any():
        # noisefill (default): replace masked pixels with a noise realisation.
        gal_masked[~mask] = noise_img_gal[~mask]

    return gal_masked, weight_map, noise_img

def make_ngmix_observation(
    gal, weight, flag, psf, wcs, rng,
    bkg_rms=None, centroid_source="hsm", wcs_full=None, ra=None, dec=None,
    blend_handling="noisefill", seg=None, object_number=None,
):
    """Build an ngmix Observation for a single galaxy epoch.

    The galaxy Jacobian origin sets where the centroid prior is centered, so
    it must sit on the object. Two ways to place it, selected by
    ``centroid_source``:

    * ``"hsm"`` (default) — re-center on the HSM adaptive-moment centroid
      measured from the stamp. Robust for **galaxies**: it follows the actual
      light and so the centroid prior (centered at the Jacobian origin) does
      not bias an object that is offset from the stamp center.
    * ``"wcs"`` — place the origin at the object's catalog sky position,
      projected through the WCS to a sub-pixel pixel offset from the stamp
      center, with no shape measurement. Better for **stars**: their HSM
      moments are noisy, so trusting the astrometry is more stable than
      re-measuring the centroid.

    Parameters
    ----------
    gal : numpy.ndarray
    weight : numpy.ndarray
    flag : numpy.ndarray
    psf : numpy.ndarray
    wcs : galsim.BaseWCS
        Local WCS Jacobian at the object position.
    rng : numpy.random.RandomState
        Random state for the noise realisations (seeded per tile for
        reproducibility).
    bkg_rms : numpy.ndarray, optional
        Per-pixel background RMS map.
    centroid_source : {"hsm", "wcs"}, optional
        How to place the galaxy Jacobian origin; the default is ``"hsm"``.
    wcs_full : astropy.wcs.WCS, optional
        Full exposure WCS for the object's CCD. Required for
        ``centroid_source="wcs"`` (ignored for ``"hsm"``).
    ra, dec : float, optional
        Object sky position in degrees. Required for
        ``centroid_source="wcs"`` (ignored for ``"hsm"``).
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment passed through to :func:`prepare_ngmix_weights`;
        the default ``"noisefill"`` is the historical behaviour.
    seg : numpy.ndarray, optional
        Segmentation map on the stamp grid. Required for
        ``blend_handling="uberseg"`` (ignored otherwise).
    object_number : int, optional
        Central object's segmentation label. Required for
        ``blend_handling="uberseg"`` (ignored otherwise).

    Returns
    -------
    ngmix.observation.Observation
    """
    psf_jacob = ngmix.Jacobian(
        row=(psf.shape[0] - 1) / 2,
        col=(psf.shape[1] - 1) / 2,
        wcs=wcs,
    )
    psf_obs = Observation(psf, jacobian=psf_jacob)

    gal_masked, weight_map, noise_img = prepare_ngmix_weights(
        gal, weight, flag, rng, bkg_rms=bkg_rms,
        blend_handling=blend_handling, seg=seg, object_number=object_number,
    )

    if centroid_source == "hsm":
        # Re-center Jacobian on HSM centroid (pixel offset from stamp center).
        # Fixes: centroid prior biases fit when galaxy is offset from stamp
        # center. Robust for galaxies; noisy for stars (use "wcs" there).
        try:
            _hsm = galsim.hsm.FindAdaptiveMom(
                galsim.Image(gal, scale=1.0), strict=False
            )
            if _hsm.error_message != "":
                raise galsim.hsm.GalSimHSMError(_hsm.error_message)
            _cen = _hsm.moments_centroid - galsim.Image(gal, scale=1.0).center
            cen_row, cen_col = _cen.y, _cen.x
        except Exception:
            cen_row, cen_col = 0.0, 0.0
    elif centroid_source == "wcs":
        # Place the origin at the catalog sky position projected through the
        # WCS — no shape measurement. Stars have noisy HSM moments, so trust
        # the astrometry instead of re-measuring the centroid.
        g_wcs = galsim.fitswcs.AstropyWCS(wcs=wcs_full)
        world_pos = galsim.CelestialCoord(
            ra * galsim.degrees, dec * galsim.degrees
        )
        pos = g_wcs.toImage(world_pos)
        cen_col = pos.x - np.round(pos.x).astype(int)
        cen_row = pos.y - np.round(pos.y).astype(int)
    else:
        raise ValueError(
            f"Unknown centroid_source '{centroid_source}'; expected"
            + " 'hsm' or 'wcs'"
        )

    gal_jacob = ngmix.Jacobian(
        row=(gal.shape[0] - 1) / 2 + cen_row,
        col=(gal.shape[1] - 1) / 2 + cen_col,
        wcs=wcs,
    )

    return Observation(
        gal_masked,
        weight=weight_map,
        jacobian=gal_jacob,
        psf=psf_obs,
        noise=noise_img,
    )

def average_multiepoch_psf(obsdict):
    """ averages psf information over multiple epochs
    we may need to do this for original psf as well
    Parameters
    ----------
    obsdict : dict
        Observation dict returned by MetacalBootstrapper.go().

    Returns
    -------
    dict
        Keys: 'g_psf', 'g_psf_err', 'T_psf', 'T_psf_err' (weighted
        averages over the epochs whose PSF fit succeeded) and 'n_epoch'
        (the number of those surviving epochs).
    """
    psf_dict = {}
    nepoch = len(obsdict['noshear'])
    n_epoch_used = 0
    wsum = 0
    g_psf_sum = np.array([0., 0.])
    g_psf_err_sum = np.array([0., 0.])
    T_psf_sum = 0
    T_psf_err_sum = 0
    for n_e in np.arange(nepoch):
        result = obsdict['noshear'][n_e].psf.meta['result']
        # ignore_failed_psf=True drops failed-PSF epochs from the galaxy
        # fit but keeps them in obsdict; their result carries only
        # flags/pars (no T/g), so skip them here too.
        if result['flags'] != 0:
            continue
        ne_wsum = obsdict['noshear'][n_e].weight.sum()

        n_epoch_used += 1
        wsum += ne_wsum
        g_psf_sum += result['g'] * ne_wsum
        g_psf_err_sum += result['g_err'] * ne_wsum
        T_psf_sum += result['T'] * ne_wsum
        T_psf_err_sum += result['T_err'] * ne_wsum

    if wsum == 0:
        raise ZeroDivisionError('Sum of weights = 0, division by zero')

    psf_dict['g_psf'] = g_psf_sum / wsum
    psf_dict['g_psf_err'] = g_psf_err_sum / wsum
    psf_dict['T_psf'] = T_psf_sum / wsum
    psf_dict['T_psf_err'] = T_psf_err_sum / wsum
    psf_dict['n_epoch'] = n_epoch_used

    return psf_dict


def make_runners(prior, flux_guess, rng):
    """Build the module's galaxy and PSF runners.

    Single source of truth for the fitter configuration (Gaussian galaxy
    and PSF models, guessers, retry counts), shared by the metacal
    bootstrap below and by validation tests that fit module-built
    observations directly.

    Parameters
    ----------
    prior : ngmix.joint_prior.PriorSimpleSep
        Priors for the fitting parameters.
    flux_guess : float
        Initial flux guess.
    rng : numpy.random.RandomState
        Random state for the guessers.

    Returns
    -------
    tuple
        (runner, psf_runner) : ngmix.runners.Runner, ngmix.runners.PSFRunner
    """
    fitter = ngmix.fitting.Fitter(model='gauss', prior=prior)
    guesser = ngmix.guessers.TPSFFluxAndPriorGuesser(rng=rng, T=0.25, prior=prior)

    psf_fitter = ngmix.fitting.Fitter(model='gauss', prior=prior)
    psf_guesser = ngmix.guessers.TFluxGuesser(rng=rng, T=0.25, prior=prior, flux=flux_guess)

    return (
        ngmix.runners.Runner(fitter=fitter, guesser=guesser, ntry=5),
        ngmix.runners.PSFRunner(fitter=psf_fitter, guesser=psf_guesser, ntry=2),
    )


def do_ngmix_metacal(
    stamp, prior, flux_guess, rng, centroid_source="hsm",
    blend_handling="noisefill", object_number=None,
):
    """Do Ngmix Metacal.

    Performs metacalibration on a single multi-epoch object and returns the
    joint shape measurement with NGMIX.

    Parameters
    ----------
    stamp : Postage_stamp
        Postage stamps for all epochs of one galaxy.
    prior : ngmix.joint_prior.PriorSimpleSep
        Priors for the fitting parameters.
    flux_guess : float
        Initial flux guess.
    rng : numpy.random.RandomState
        Random state for guesses and priors.
    centroid_source : {"hsm", "wcs"}, optional
        How to place the galaxy Jacobian origin; passed through to
        :func:`make_ngmix_observation`. The default is ``"hsm"`` (HSM
        adaptive-moment centroid); ``"wcs"`` uses the catalog sky position
        projected through the WCS — see that function for the star-vs-galaxy
        rationale.
    blend_handling : {"noisefill", "uberseg"}, optional
        Neighbour treatment passed through to
        :func:`make_ngmix_observation`; the default ``"noisefill"`` is the
        historical behaviour. ``"uberseg"`` consumes ``stamp.segs`` and
        ``object_number``.
    object_number : int, optional
        Central object's segmentation label (its SExtractor ``NUMBER``).
        Required for ``blend_handling="uberseg"`` (ignored otherwise).

    Returns
    -------
    tuple
        (resdict, psf_res) where resdict is the MetacalBootstrapper result
        dict and psf_res is the averaged PSF dict from average_multiepoch_psf.
    """
    n_epoch = len(stamp.gals)
    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    gal_obs_list = ObsList()
    for n_e in range(n_epoch):
        bkg_rms = stamp.bkg_rms[n_e] if len(stamp.bkg_rms) > n_e else None
        gal_obs = make_ngmix_observation(
            stamp.gals[n_e],
            stamp.weights[n_e],
            stamp.flags[n_e],
            stamp.psfs[n_e],
            stamp.jacobs[n_e],
            rng,
            bkg_rms=bkg_rms,
            centroid_source=centroid_source,
            wcs_full=stamp.wcs[n_e] if n_e < len(stamp.wcs) else None,
            ra=stamp.ra[n_e] if n_e < len(stamp.ra) else None,
            dec=stamp.dec[n_e] if n_e < len(stamp.dec) else None,
            blend_handling=blend_handling,
            seg=stamp.segs[n_e] if n_e < len(stamp.segs) else None,
            object_number=object_number,
        )
        gal_obs_list.append(gal_obs)

    runner, psf_runner = make_runners(prior, flux_guess, rng)

    metacal_pars = {
        'types': ['noshear', '1p', '1m', '2p', '2m'],
        'step': 0.01,
        'psf': 'fitgauss',
        'fixnoise': True,
        'use_noise_image': True,
    }

    boot = ngmix.metacal.MetacalBootstrapper(
        runner=runner,
        psf_runner=psf_runner,
        ignore_failed_psf=True,
        rng=rng,
        **metacal_pars,
    )
    resdict, obsdict = boot.go(gal_obs_list)
    psf_res = average_multiepoch_psf(obsdict)
    return resdict, psf_res
