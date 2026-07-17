"""NGMIX.

This module contains a class for ngmix shape measurement.

:Authors: Lucie Baumont, Axel Guinot

"""

import os
import re
from typing import NamedTuple

import ngmix
import galsim
import numpy as np
from astropy.io import fits
from cs_util import size as cs_size
from modopt.math.stats import sigma_mad
from ngmix.observation import Observation, ObsList
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io

METACAL_TYPES = ('noshear', '1p', '1m', '2p', '2m')

# Bit set in mcal_flags when a metacal type carries no fit result at all
# (absent entry or missing 'flags' key), or when a nominally clean fit
# returns non-finite shear. Absence of evidence of success is failure,
# not success: previously an absent flag defaulted to 0 and a silently
# failing fitter produced an all-zero flag column.
FLAG_NO_RESULT = 2**30

# Noise budget for the PSF observation's flat weight map (psf_wt =
# 1/PSF_NOISE**2). Mirrors the esheldon/aguinot pattern (sigma ~ 1e-5/1e-6);
# 1e-5 is the value Axel Guinot's #749 reproduction used. The fit is driven
# by the *relative* weighting of the likelihood vs the prior, so the exact
# value is non-critical once it is finite (validated on the digital twin: the
# recovered PSF shape/size are flat across 1e-4..1e-6). See
# make_ngmix_observation.
PSF_NOISE = 1e-5


class MetacalResult(NamedTuple):
    """Return of :func:`do_ngmix_metacal`: the metacal fit plus both PSFs.

    A NamedTuple so the two PSF families are named, not positional — guarding
    against a reconv/orig transposition at call sites. Still a plain tuple, so
    ``resdict, psf_res, psf_orig_res = do_ngmix_metacal(...)`` keeps working.

    Attributes
    ----------
    resdict : dict
        MetacalBootstrapper result dict (one entry per metacal type).
    reconv : dict
        Averaged metacal *reconvolution*-kernel PSF (round, enlarged):
        :func:`average_multiepoch_psf`.
    orig : dict
        Averaged *original* image-PSF (psfex/mccd, its true shape and size):
        :func:`average_original_psf`.
    """

    resdict: dict
    reconv: dict
    orig: dict


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
        [res.get(name, {}).get('flags', FLAG_NO_RESULT)
         for name in METACAL_TYPES]
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


def position_seed(ra, dec, ccd):
    """Deterministic RNG seed from an object's sky position (ngmix#796).

    For image-simulation m-bias with the Pujol estimator, the same scene is
    simulated under different applied shears ("image branches") and the shear
    response is read from the branch difference of the SAME objects. Metacal's
    ``fixnoise`` adds a counter-noise realisation drawn from an RNG; if that RNG
    is seeded per tile, an object gets *different* added noise in different
    branches (detection order differs), and the noise fails to cancel in the
    difference, inflating ``sigma_m``. Seeding the per-object RNG from sky
    position instead makes the same object draw the same added noise (and the
    same fit guesses) in every branch, so both cancel and the m-bias error
    shrinks. Off in production; a knob for the sim path only.

    Box math (kept exactly as Fabian's issue #796)::

        box_x = floor(ra  * 3600 / 3) + (ccd + 1)
        box_y = floor(dec * 3600 / 3) + (ccd + 2)

    The 3-arcsec boxes make the seed robust to detection-order changes and to
    small centroid jitter within a box: the same physical object lands in the
    same box across branches and so gets the same seed. The ``+ccd`` offsets
    disambiguate exposure overlap.

    **Deviation from #796.** Fabian combines the boxes as ``box_x + box_y``,
    which collides along anti-diagonals — e.g. boxes ``(10, 20)`` and
    ``(11, 19)`` both give seed 30, so two distinct objects share a noise
    stream. We instead combine them with a Cantor pairing (a bijection on
    non-negative integers), after a zig-zag fold that maps signed box indices
    (dec can be negative) onto non-negative ones. The result is reduced mod
    ``2**32`` to land in the valid ``numpy.random.RandomState`` seed range
    ``[0, 2**32)``. The box math is untouched; only the collision-prone sum is
    replaced.

    Parameters
    ----------
    ra, dec : float
        Object sky position in degrees (first-epoch coordinates; all epochs of
        one object share them).
    ccd : int
        CCD number of the first epoch.

    Returns
    -------
    int
        Seed in ``[0, 2**32)`` for ``numpy.random.RandomState``.
    """
    box_x = int(np.floor((ra * 3600) / 3) + (ccd + 1))
    box_y = int(np.floor((dec * 3600) / 3) + (ccd + 2))
    # Zig-zag fold signed -> non-negative (0,-1,1,-2,... -> 0,1,2,3,...) so the
    # Cantor pairing, which is a bijection on the non-negative integers, stays
    # injective for southern (dec<0) boxes.
    zx = 2 * box_x if box_x >= 0 else -2 * box_x - 1
    zy = 2 * box_y if box_y >= 0 else -2 * box_y - 1
    cantor = (zx + zy) * (zx + zy + 1) // 2 + zy
    return cantor % (2 ** 32)


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
        self.jacobs = []
        # Per-epoch full WCS and the object's sky position, used only by the
        # "wcs" centroid source (skipped for the default "hsm" path).
        self.wcs = []
        self.ra = []
        self.dec = []
        # CCD number of the first epoch, used only to build the per-object
        # position seed (see :func:`position_seed`, ngmix#796). ``None`` when
        # position seeding is off, so the field costs nothing on the hot path.
        self.ccd = None
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
        self.bkg_vign_cat = SqliteDict(bkg_vignet_path) if bkg_vignet_path is not None else None
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
        if self.bkg_vign_cat is not None:
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
    seed_from_position : bool, optional
        If ``True``, replace the tile-level RNG with a per-object RNG seeded
        from the object's sky position (:func:`position_seed`) inside the
        object loop, so metacal's ``fixnoise`` counter-noise and the fit
        guesses cancel across Pujol image-simulation branches (ngmix#796). The
        default ``False`` leaves the production path byte-identical. See
        :func:`position_seed` for the physics and the seed construction.

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
        bkg_sub=True,
        centroid_source="hsm",
        seed_from_position=False,
        metacal_psf="fitgauss",
    ):

        # Base count = catalogue + vignets, excluding the f_wcs headers (passed
        # separately). One fewer when the background vignet is absent
        # (``bkg_sub=False``, image sims). An extra trailing slot may carry the
        # optional background-rms vignet (#779), so both counts are valid.
        n_base = 6 if bkg_sub else 5
        if len(input_file_list) not in {n_base, n_base + 1}:
            raise IndexError(
                f"Input file list has length {len(input_file_list)},"
                + f" required is {n_base} or {n_base + 1}"
            )

        self._tile_cat_path = input_file_list[0]
        if bkg_sub:
            bkg_path, psf_path, weight_path, flag_path = (
                input_file_list[2], input_file_list[3],
                input_file_list[4], input_file_list[5],
            )
        else:
            bkg_path, psf_path, weight_path, flag_path = (
                None, input_file_list[2],
                input_file_list[3], input_file_list[4],
            )
        bkg_rms_vignet_path = (
            input_file_list[n_base]
            if len(input_file_list) == n_base + 1
            else None
        )
        self._vignet_cat = Vignet(
            input_file_list[1],
            bkg_path,
            psf_path,
            weight_path,
            flag_path,
            f_wcs_path,
            bkg_rms_vignet_path,
        )

        self._output_dir = output_dir
        self._file_number_string = file_number_string

        self._zero_point = zero_point
        self._pixel_scale = pixel_scale

        self._f_wcs_path = f_wcs_path

        self._save_batch = save_batch
        self._id_obj_min = id_obj_min
        self._id_obj_max = id_obj_max
        self._bkg_sub = bkg_sub
        self._centroid_source = centroid_source
        self._seed_from_position = seed_from_position
        self._metacal_psf = metacal_psf

        self._w_log = w_log

        # Initiatlise random generator
        seed = int(''.join(re.findall(r'\d+', self._file_number_string)))
        self._rng = np.random.RandomState(seed)
        self._w_log.info(f'Random generator initialisation seed = {seed}')
        if self._seed_from_position:
            self._w_log.info(
                'SEED_FROM_POSITION on: per-object RNG seeded from sky position'
                ' for Pujol noise cancellation (image sims, ngmix#796)'
            )

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
            Compiled results ready to be written to a file.

            Two PSF column families — each carrying ellipticity *and* size,
            for *different* PSFs (shapepipe#749). See the function docstrings
            for what each family IS:

            * ``*_psf_orig`` (``g1``/``g2`` + ``*_err``, ``T``) — the original
              image PSF, fit by :func:`average_original_psf`.
            * ``*_psf_reconv`` — the metacal reconvolution kernel, fit by
              :func:`average_multiepoch_psf`.

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
            # galaxy
            'g1',
            'g1_err',
            'g2',
            'g2_err',
            'T',
            'T_err',
            'flux',
            'flux_err',
            's2n',
            'mag',
            'mag_err',
            'flags',
            'mcal_flags',
            # original image PSF (psfex/mccd), fit by average_original_psf
            'g1_psf_orig',
            'g2_psf_orig',
            'g1_err_psf_orig',
            'g2_err_psf_orig',
            'T_psf_orig',
            'T_err_psf_orig',
            # metacal reconvolution kernel, fit by average_multiepoch_psf
            'g1_psf_reconv',
            'g2_psf_reconv',
            'g1_err_psf_reconv',
            'g2_err_psf_reconv',
            'T_psf_reconv',
            'T_err_psf_reconv',
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
                # 1-5 retry count, so the column is named accordingly. Fits
                # that fail before nfev is ever set (e.g. the NaN-filled
                # branch above) fall back to the same -1 sentinel, keeping
                # the column int64 (FITS 'K') across every batch: nfev_fit
                # is int64 with -1 meaning failed/absent.
                output_dict[name]["nfev_fit"].append(
                    fit.get("nfev", -1)
                )
                # The two PSF families are object-level (one value per
                # object, not per shear type) and self-named: every key
                # below is copied straight through from compile-loop input to
                # output, so the column name *is* the value's provenance.
                #   *_psf_orig   = original image PSF (average_original_psf)
                #   *_psf_reconv = reconvolution kernel (average_multiepoch_psf)
                for psf_key in (
                    'g1_psf_orig', 'g2_psf_orig',
                    'g1_err_psf_orig', 'g2_err_psf_orig',
                    'T_psf_orig', 'T_err_psf_orig',
                    'g1_psf_reconv', 'g2_psf_reconv',
                    'g1_err_psf_reconv', 'g2_err_psf_reconv',
                    'T_psf_reconv', 'T_err_psf_reconv',
                ):
                    output_dict[name][psf_key].append(results[idx][psf_key])

                output_dict[name]["g1"].append(g[0])
                output_dict[name]["g2"].append(g[1])
                output_dict[name]["g1_err"].append(np.sqrt(g_cov[0, 0]))
                output_dict[name]["g2_err"].append(np.sqrt(g_cov[1, 1]))
                output_dict[name]["T"].append(T_gal)
                output_dict[name]["T_err"].append(T_gal_err)
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

    def log_mean_ellipticity(self):
        """Log mean ellipticity from NOSHEAR HDU to the run log.

        Reports <e1>, <e2> with standard errors for all objects and for
        objects passing the default metacal cuts (flags==0, mcal_flags==0,
        10 < SNR < 500, T/Tpsf > 0.5).
        """
        output_path = self.get_output_path(self._output_dir)
        try:
            with fits.open(output_path) as hdul:
                d = hdul['NOSHEAR'].data
                g1 = d['g1'].astype(float)
                g2 = d['g2'].astype(float)
                flags = d['flags']
                mcal_flags = d['mcal_flags']
                s2n = d['s2n'].astype(float)
                T = d['T'].astype(float)
                Tpsf = d['Tpsf'].astype(float)
        except Exception as e:
            self._w_log.warning(f"Could not compute mean ellipticity: {e}")
            return

        n_total = len(g1)
        if n_total == 0:
            self._w_log.info("Mean ellipticity: no objects in output catalogue")
            return

        def _log_stats(g1_sel, g2_sel, label):
            n = len(g1_sel)
            if n == 0:
                self._w_log.info(f"  {label}: 0 objects")
                return
            mean_g1 = g1_sel.mean()
            mean_g2 = g2_sel.mean()
            err_g1 = g1_sel.std() / np.sqrt(n)
            err_g2 = g2_sel.std() / np.sqrt(n)
            self._w_log.info(
                f"  {label} (N={n}):"
                f"  <e1> = {mean_g1:+.4e} +/- {err_g1:.4e},"
                f"  <e2> = {mean_g2:+.4e} +/- {err_g2:.4e}"
            )

        self._w_log.info(f"Mean ellipticity (NOSHEAR, N_total={n_total}):")
        _log_stats(g1, g2, "no cuts")

        with np.errstate(invalid='ignore'):
            mask = (
                (flags == 0)
                & (mcal_flags == 0)
                & (s2n >= 10.0)
                & (s2n <= 500.0)
                & (T / Tpsf >= 0.5)
            )
        _log_stats(g1[mask], g2[mask], "SNR in [10, 500], T/Tpsf > 0.5")

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
        n_flagged = 0
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

            stamp = prepare_postage_stamps(vignet_cat, obj_id, i_tile, tile_cat, self._bkg_sub)

            if len(stamp.gals) == 0:
                n_no_epoch += 1
                continue

            # Position-seeded per-object RNG for Pujol noise cancellation in
            # image sims (ngmix#796): the same object gets the same fixnoise
            # counter-noise and fit guesses in every shear branch, so both
            # cancel in the branch difference. The prior is rebuilt from the
            # same per-object RNG because the guesser draws its initial guess
            # via prior.sample() (ngmix guessers.py), which consumes the RNG the
            # prior was CONSTRUCTED with — so a per-object rng alone would leave
            # the guess drawing from the shared tile stream and break
            # cancellation. Off in production, where the single tile-level
            # self._rng and the tile-level prior carry the whole loop.
            if self._seed_from_position:
                obj_rng = np.random.RandomState(
                    position_seed(stamp.ra[0], stamp.dec[0], stamp.ccd)
                )
                obj_prior = get_prior(self._pixel_scale, obj_rng)
            else:
                obj_rng = self._rng
                obj_prior = prior

            try:
                flux_guess = (
                    tile_cat.flux[i_tile]
                    if tile_cat.flux is not None
                    else 1.0
                )
                res, psf_res, psf_orig_res = do_ngmix_metacal(
                    stamp,
                    obj_prior,
                    flux_guess,
                    obj_rng,
                    centroid_source=self._centroid_source,
                    metacal_psf=self._metacal_psf,
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
                if res.get(k, {}).get('flags', FLAG_NO_RESULT) != 0
            )
            res['mcal_flags'] = get_mcal_flags(res)
            # A "successful" fit whose noshear shear is non-finite is a
            # silent failure (sentinel output with flags == 0): flag it.
            if res['mcal_flags'] == 0 and not np.all(np.isfinite(
                np.asarray(
                    res.get('noshear', {}).get('g', (np.nan, np.nan)),
                    dtype=float,
                )
            )):
                res['mcal_flags'] = FLAG_NO_RESULT
            if res['mcal_flags'] != 0:
                n_flagged += 1
            # Two distinct PSF families (shapepipe#749), each carrying its own
            # ellipticity AND size: the metacal reconvolution kernel (psf_res)
            # and the original image PSF (psf_orig_res). Tag both into res from
            # ONE template so the families stay symmetric — same quantities,
            # parallel ``*_psf_{family}`` names — and cannot drift apart by a
            # hand-edit to one and not the other.
            for family, psf in (("reconv", psf_res), ("orig", psf_orig_res)):
                res[f"g1_psf_{family}"] = psf["g_psf"][0]
                res[f"g2_psf_{family}"] = psf["g_psf"][1]
                res[f"g1_err_psf_{family}"] = psf["g_psf_err"][0]
                res[f"g2_err_psf_{family}"] = psf["g_psf_err"][1]
                res[f"T_psf_{family}"] = psf["T_psf"]
                res[f"T_err_psf_{family}"] = psf["T_psf_err"]
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

        if count > 0 and n_fitted == 0:
            raise RuntimeError(
                f'ngmix: all {count} objects failed the metacal fit'
                ' (0 fitted) -- likely an upstream library/PSF problem,'
                ' not a data property; aborting instead of writing an'
                ' empty catalogue.'
            )
        if n_fitted > 0 and n_flagged == n_fitted:
            self._w_log.error(
                f'ngmix: 100% of {n_fitted} fitted objects carry nonzero'
                ' mcal_flags -- the metacal fit failed wholesale; outputs'
                ' are unusable.'
            )

        vignet_cat.close()

        # Put all results together
        res_dict = self.compile_results(final_res)

        # Save results
        self.save_results(res_dict)

        # Log mean ellipticity statistics
        self.log_mean_ellipticity()

def prepare_postage_stamps(vignet, obj_id, i_tile, tile_cat, bkg_sub=True):
    # define per-object lists of individual exposures to go into ngmix
    stamp = Postage_stamp(bkg_sub=bkg_sub)
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

        # Skip epochs where sigma_mad=0 (CCD edge: mostly-zero stamp).
        # prepare_ngmix_weights divides by sig_noise; zero sigma causes NaN/inf
        # weights that corrupt GalSim's C-level FFT allocations.
        if not sigma_mad(gal_vign_sub_bkg) > 0:
            continue

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
        # CCD of the first surviving epoch — Fabian's coord_list[0] convention
        # for the position seed (ngmix#796). All epochs of one object share the
        # ra/dec above, so first-epoch CCD pins one deterministic seed stream.
        if stamp.ccd is None:
            stamp.ccd = int(ccd_n)

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

    gauss_win = galsim.Gaussian(sigma=cs_size.T_to_sigma(guess[4]), flux=guess[5])
    gauss_win = gauss_win.shear(g1=guess[2], g2=guess[3])
    gauss_win = gauss_win.drawImage(
        nx=img_shape[0], ny=img_shape[1], scale=pixel_scale
    ).array

    m_weight = weight[gauss_win < thresh * sig_tmp] != 0

    sig_noise = sigma_mad(gal[gauss_win < thresh * sig_tmp][m_weight])

    return sig_noise

def prepare_ngmix_weights(gal, weight, flag, rng, bkg_rms=None):
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

    Returns
    -------
    numpy.ndarray
        Galaxy image with masked pixels replaced by noise.
    numpy.ndarray
        Variance map for NGMIX.
    numpy.ndarray
        Noise image.
    """
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

    # Guard: sig_noise=0 means galaxy is at the CCD edge (mostly-zero stamp).
    # Division by zero would make weight_map NaN/inf, crashing GalSim C code.
    # np.all keeps the guard valid for the per-pixel bkg_rms path (any zero
    # pixel would produce the same NaN weight the scalar case guards against).
    if not np.all(sig_noise > 0):
        return np.zeros_like(gal), np.zeros_like(weight_map), np.zeros_like(gal)

    noise_img = rng.standard_normal(gal.shape) * sig_noise
    noise_img_gal = rng.standard_normal(gal.shape) * sig_noise

    gal_masked = np.copy(gal)
    if (~mask).any():
        gal_masked[~mask] = noise_img_gal[~mask]

    return gal_masked, weight_map, noise_img

def make_ngmix_observation(
    gal, weight, flag, psf, wcs, rng,
    bkg_rms=None, centroid_source="hsm", wcs_full=None, ra=None, dec=None,
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

    Returns
    -------
    ngmix.observation.Observation
    """
    psf_jacob = ngmix.Jacobian(
        row=(psf.shape[0] - 1) / 2,
        col=(psf.shape[1] - 1) / 2,
        wcs=wcs,
    )
    # A model fit always needs a weight: without one ngmix defaults to unit
    # weights on a unit-flux PSF stamp, so the galaxy GPriorBA(0.4) prior
    # handed to the diagnostic PSF fitter swamps the (tiny) likelihood and the
    # recovered PSF shape collapses toward the prior (#749). PSFEx/MCCD models
    # already carry a little noise, so a finite flat weight matching the
    # esheldon/aguinot noise budget (psf_noise ~ 1e-5) suffices to restore the
    # likelihood — no explicit stamp-noise injection needed (#774). This is the
    # PSF analogue of the galaxy weight built just below; force float so an
    # integer-typed stamp cannot truncate the weight (matching that build).
    # Metacal's own PSF fit is prior-free (AdmomFitter) and so is insensitive
    # to this flat weight's scale, leaving the calibration untouched.
    psf_wt = np.full_like(psf, 1.0 / PSF_NOISE ** 2, dtype=float)
    psf_obs = Observation(psf, weight=psf_wt, jacobian=psf_jacob)

    gal_masked, weight_map, noise_img = prepare_ngmix_weights(
        gal, weight, flag, rng, bkg_rms=bkg_rms
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

def _average_psf_fits(results_and_weights):
    """Weight-average a set of per-epoch ngmix PSF-fit results.

    Shared core for both PSF families this module exports: the metacal
    reconvolution kernel (:func:`average_multiepoch_psf`) and the original
    image PSF (:func:`average_original_psf`). Epochs whose PSF fit failed
    (``flags != 0``, carrying only flags/pars and no T/g) are dropped.

    Parameters
    ----------
    results_and_weights : iterable of (dict, float)
        Per-epoch ``(result, weight)`` pairs, where ``result`` is an ngmix
        Fitter result with keys ``flags``, ``g``, ``g_err``, ``T``,
        ``T_err`` and ``weight`` is the epoch's averaging weight.

    Returns
    -------
    dict
        Keys ``g_psf``, ``g_psf_err``, ``T_psf``, ``T_psf_err`` (weighted
        averages over the surviving epochs) and ``n_epoch`` (their count).
    """
    n_epoch_used = 0
    wsum = 0
    g_psf_sum = np.array([0., 0.])
    g_psf_err_sum = np.array([0., 0.])
    T_psf_sum = 0
    T_psf_err_sum = 0
    for result, weight in results_and_weights:
        if result['flags'] != 0:
            continue
        n_epoch_used += 1
        wsum += weight
        g_psf_sum += result['g'] * weight
        g_psf_err_sum += result['g_err'] * weight
        T_psf_sum += result['T'] * weight
        T_psf_err_sum += result['T_err'] * weight

    if wsum == 0:
        raise ZeroDivisionError('Sum of weights = 0, division by zero')

    return {
        'g_psf': g_psf_sum / wsum,
        'g_psf_err': g_psf_err_sum / wsum,
        'T_psf': T_psf_sum / wsum,
        'T_psf_err': T_psf_err_sum / wsum,
        'n_epoch': n_epoch_used,
    }


def average_multiepoch_psf(obsdict):
    """Average the metacal *reconvolution* PSF over epochs.

    The PSF carried by each metacal observation (``obs.psf``) is the
    Gaussian reconvolution kernel that metacal fit and convolved back in —
    round by construction and slightly enlarged relative to the original
    PSF. This is the kernel defining the sheared galaxy images, exported to
    the reconvolution-kernel columns
    (``NGMIX_G1/G2_PSF_RECONV``, ``NGMIX_T_PSF_RECONV``). The independent fit
    of the *original* image PSF is :func:`average_original_psf`.

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
    # ignore_failed_psf=True drops failed-PSF epochs from the galaxy fit but
    # keeps them in obsdict; _average_psf_fits skips them on flags != 0.
    return _average_psf_fits(
        (obs.psf.meta['result'], obs.weight.sum())
        for obs in obsdict['noshear']
    )


def average_original_psf(gal_obs_list, psf_runner):
    """Fit and average the *original* image PSF over epochs.

    The original PSF is the psfex/mccd model stamp handed to ngmix
    (``gal_obs.psf``), fit here with the same ``psf_runner`` (hence the same
    fit prior and guesser) used inside metacal, but on the PSF *before*
    metacal's reconvolution. Exported to the original-PSF columns
    (``NGMIX_G1/G2_PSF_ORIG``, ``NGMIX_T_PSF_ORIG``). Distinct from the
    reconvolution-kernel fit (:func:`average_multiepoch_psf`): the original
    PSF retains its true ellipticity and size, whereas the reconvolution
    kernel is round and enlarged by construction. This is the PSF whose true
    shape and size enter object-wise PSF-leakage diagnostics.

    Weighted by the raw galaxy inverse variance (``gal_obs.weight.sum()`` per
    epoch); :func:`average_multiepoch_psf` uses the same scheme but the
    fixnoise-combined metacal-image weight instead.

    The fit runs on a *copy* of each PSF observation so ``gal_obs.psf`` —
    the object metacal later deep-copies and consumes via
    ``boot.go(gal_obs_list)`` — is never mutated. ``PSFRunner.go`` sets both
    ``.meta['result']`` and, on success, the ``.gmix`` attribute of the
    observation it fits; were that ``gal_obs.psf`` itself, the stray gmix
    would survive metacal's deep copy and be reused as the
    ``MetacalFitGaussPSF`` fallback when its own admom+ML PSF fits both fail,
    silently rescuing objects the base branch dropped (``BootPSFFailure``) and
    changing the galaxy/shear result set. Fitting a copy closes this
    PSF-aliasing channel; combined with :func:`do_ngmix_metacal` seeding this
    pre-fit from a *snapshot* of the metacal RNG (so the pre-fit does not
    advance it), the add-column refactor stays bit-identical on the galaxy
    results.

    Parameters
    ----------
    gal_obs_list : ngmix.observation.ObsList
        Per-epoch galaxy observations; each ``gal_obs.psf`` is the original
        (pre-metacal) PSF observation to fit, with no further ``.psf`` of
        its own so the runner fits the stamp itself. Left pristine.
    psf_runner : ngmix.runners.PSFRunner
        The module's PSF runner (built by :func:`make_runners` from the
        shared ``prior``).

    Returns
    -------
    dict
        Same keys as :func:`average_multiepoch_psf`.
    """
    def fit(gal_obs):
        # Fit a COPY so gal_obs.psf stays pristine for metacal — see docstring.
        # Failed fits keep flags != 0 and are dropped by _average_psf_fits.
        psf_obs = gal_obs.psf.copy()
        psf_runner.go(psf_obs)
        return psf_obs.meta['result'], gal_obs.weight.sum()

    return _average_psf_fits(fit(gal_obs) for gal_obs in gal_obs_list)


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
    stamp, prior, flux_guess, rng, centroid_source="hsm", metacal_psf="fitgauss"
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
    metacal_psf : {"fitgauss", "gauss", "dilate", "azgauss"}, optional
        The metacal reconvolution-kernel scheme passed to ngmix as
        ``metacal_pars['psf']``. The default ``"fitgauss"`` fits a Gaussian to
        the PSF and rounds it (``MetacalFitGaussPSF``); ``"gauss"`` reconvolves
        with a fixed round Gaussian sized from the PSF (``MetacalGaussPSF``);
        ``"dilate"`` dilates the original PSF; ``"azgauss"`` (ngmix >= 2.4.1) is
        a noise-robust variant of ``"gauss"``. This kernel is a SEPARATE object
        from the PSF-model fit (the ``psf_runner`` above) — it only sets the
        round PSF that metacal reconvolves with after shearing, so it moves the
        metacal *response* (and therefore the recovered shear) but never reaches
        the deconvolution, which is by the PSF image.

    Returns
    -------
    MetacalResult
        Named 3-tuple ``(resdict, reconv, orig)``: the MetacalBootstrapper
        result dict, the averaged metacal *reconvolution*-kernel PSF dict
        (:func:`average_multiepoch_psf`), and the averaged *original* image-PSF
        dict (:func:`average_original_psf`). The two PSF dicts share keys but
        describe different PSFs; the named fields guard against transposing
        them. Unpacks positionally as ``resdict, psf_res, psf_orig_res``.
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
        )
        gal_obs_list.append(gal_obs)

    runner, psf_runner = make_runners(prior, flux_guess, rng)

    # Fit the ORIGINAL (psfex/mccd) PSF before metacal reconvolves it. Use a
    # psf_runner seeded from a snapshot of rng's state so this prefit does NOT
    # advance the rng that boot.go consumes below — keeping the galaxy/shear
    # results bit-identical to the no-prefit branch. average_original_psf fits a
    # COPY of each gal_obs.psf, so gal_obs_list also reaches boot.go pristine.
    prefit_rng = np.random.RandomState()
    prefit_rng.set_state(rng.get_state())
    _, prefit_psf_runner = make_runners(prior, flux_guess, prefit_rng)
    psf_orig_res = average_original_psf(gal_obs_list, prefit_psf_runner)

    metacal_pars = {
        'types': ['noshear', '1p', '1m', '2p', '2m'],
        'step': 0.01,
        'psf': metacal_psf,
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
    del obsdict

    # Each FitModel in resdict retains the full metacal observations
    # (``.obs``, a plain attribute set by FitModel._set_obs) plus the pixel
    # arrays derived from them (``._pixels_list``) — several MB per object.
    # Results accumulate until the SAVE_BATCH flush, and compile_results
    # reads only scalar dict keys, so release the arrays now to keep
    # resident memory flat.
    for _fm in resdict.values():
        if hasattr(_fm, 'obs'):
            _fm.obs = None
        if hasattr(_fm, '_pixels_list'):
            _fm._pixels_list = None

    return MetacalResult(resdict=resdict, reconv=psf_res, orig=psf_orig_res)
