"""NGMIX.

This module contains a class for ngmix shape measurement.

:Author: Axel Guinot

"""

import re

import galsim
import ngmix
import numpy as np
from astropy.io import fits
from modopt.math.stats import sigma_mad
from ngmix.fitting import LMSimple
from ngmix.observation import MultiBandObsList, Observation, ObsList
from numpy.random import uniform as urand
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io


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
    id_obj_min : int, optional
        First galaxy ID to process, not used if the value is set to ``-1``;
        the default is ``-1``
    id_obj_max : int, optional
        Last galaxy ID to process, not used if the value is set to ``-1``;
        the default is ``-1``

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
        id_obj_min=-1,
        id_obj_max=-1,
    ):

        if len(input_file_list) != 6:
            raise IndexError(
                f"Input file list has length {len(input_file_list)},"
                + " required is 6"
            )

        self._tile_cat_path = input_file_list[0]
        self._gal_vignet_path = input_file_list[1]
        self._bkg_vignet_path = input_file_list[2]
        self._psf_vignet_path = input_file_list[3]
        self._weight_vignet_path = input_file_list[4]
        self._flag_vignet_path = input_file_list[5]

        self._output_dir = output_dir
        self._file_number_string = file_number_string

        self._zero_point = zero_point
        self._pixel_scale = pixel_scale

        self._f_wcs_path = f_wcs_path
        self._id_obj_min = id_obj_min
        self._id_obj_max = id_obj_max

        self._w_log = w_log

        # Initiatlise random generator
        seed = int("".join(re.findall(r"\d+", self._file_number_string)))
        np.random.seed(seed)
        self._w_log.info(f"Random generator initialisation seed = {seed}")

    @classmethod
    def MegaCamFlip(self, vign, ccd_nb):
        """Flip for MegaCam.

        MegaCam has CCDs that are upside down. This function flips the
        postage stamps in these CCDs.

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

    def get_prior(self):
        """Get Prior.

        Return prior for the different parameters.

        Returns
        -------
        ngmix.priors
            Priors for the different parameters

        """
        # Prior on ellipticity. Details do not matter, as long
        # as it regularizes the fit. From Bernstein & Armstrong 2014
        g_sigma = 0.4
        g_prior = ngmix.priors.GPriorBA(g_sigma)

        # 2-d Gaussian prior on the center row and column center
        # (relative to the center of the jacobian, which
        # would be zero) and the sigma of the Gaussians.
        # Units same as jacobian, probably arcsec
        row, col = 0.0, 0.0
        row_sigma, col_sigma = self._pixel_scale, self._pixel_scale
        cen_prior = ngmix.priors.CenPrior(row, col, row_sigma, col_sigma)

        # Size prior. Instead of flat, two-sided error function (TwoSidedErf)
        # could be used
        Tminval = -10.0  # arcsec squared
        Tmaxval = 1.0e6
        T_prior = ngmix.priors.FlatPrior(Tminval, Tmaxval)

        # Flux prior. Bounds need to make sense for
        # images in question
        Fminval = -1.0e4
        Fmaxval = 1.0e9
        F_prior = ngmix.priors.FlatPrior(Fminval, Fmaxval)

        # Joint prior, combine all individual priors
        prior = ngmix.joint_prior.PriorSimpleSep(
            cen_prior, g_prior, T_prior, F_prior
        )

        return prior

    def compile_results(self, results):
        """Compile Results.

        Prepare the results of NGMIX before saving.

        Parameters
        ----------
        results : dict
            Results of NGMIX metacal

        Returns
        -------
        dict
            Compiled results ready to be written to a file

        Raises
        ------
        KeyError
            If SNR key not found

        """
        names = ["1m", "1p", "2m", "2p", "noshear"]
        names2 = [
            "id",
            "n_epoch_model",
            "moments_fail",
            "ntry_fit",
            "g1_psfo_ngmix",
            "g2_psfo_ngmix",
            "T_psfo_ngmix",
            "g1_err_psfo_ngmix",
            "g2_err_psfo_ngmix",
            "T_err_psfo_ngmix",
            "g1",
            "g1_err",
            "g2",
            "g2_err",
            "T",
            "T_err",
            "Tpsf",
            "g1_psf",
            "g2_psf",
            "flux",
            "flux_err",
            "s2n",
            "mag",
            "mag_err",
            "flags",
            "mcal_flags",
        ]
        output_dict = {k: {kk: [] for kk in names2} for k in names}
        for idx in range(len(results)):
            for name in names:

                mag = (
                    -2.5 * np.log10(results[idx][name]["flux"])
                    + self._zero_point
                )
                mag_err = np.abs(
                    -2.5
                    * results[idx][name]["flux_err"]
                    / (results[idx][name]["flux"] * np.log(10))
                )

                output_dict[name]["id"].append(results[idx]["obj_id"])
                output_dict[name]["n_epoch_model"].append(
                    results[idx]["n_epoch_model"]
                )
                output_dict[name]["moments_fail"].append(
                    results[idx]["moments_fail"]
                )
                output_dict[name]["ntry_fit"].append(results[idx][name]["ntry"])
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
                output_dict[name]["T_psfo_ngmix"].append(results[idx]["T_PSFo"])
                output_dict[name]["T_err_psfo_ngmix"].append(
                    results[idx]["T_err_PSFo"]
                )
                output_dict[name]["g1"].append(results[idx][name]["g"][0])
                output_dict[name]["g1_err"].append(
                    results[idx][name]["pars_err"][2]
                )
                output_dict[name]["g2"].append(results[idx][name]["g"][1])
                output_dict[name]["g2_err"].append(
                    results[idx][name]["pars_err"][3]
                )
                output_dict[name]["T"].append(results[idx][name]["T"])
                output_dict[name]["T_err"].append(results[idx][name]["T_err"])
                output_dict[name]["Tpsf"].append(results[idx][name]["Tpsf"])
                output_dict[name]["g1_psf"].append(
                    results[idx][name]["gpsf"][0]
                )
                output_dict[name]["g2_psf"].append(
                    results[idx][name]["gpsf"][1]
                )
                output_dict[name]["flux"].append(results[idx][name]["flux"])
                output_dict[name]["flux_err"].append(
                    results[idx][name]["flux_err"]
                )
                output_dict[name]["mag"].append(mag)
                output_dict[name]["mag_err"].append(mag_err)

                if "s2n" in results[idx][name]:
                    output_dict[name]["s2n"].append(results[idx][name]["s2n"])
                elif "s2n_r" in results[idx][name]:
                    output_dict[name]["s2n"].append(results[idx][name]["s2n_r"])
                else:
                    raise KeyError("No SNR key (s2n, s2n_r) found in results")

                output_dict[name]["flags"].append(results[idx][name]["flags"])
                output_dict[name]["mcal_flags"].append(
                    results[idx]["mcal_flags"]
                )

        return output_dict

    def save_results(self, output_dict):
        """Save Results.

        Save the results into a FITS file.

        Parameters
        ----------
        output_dict
            Dictionary containing the results

        """
        output_name = f"{self._output_dir}/ngmix{self._file_number_string}.fits"

        f = file_io.FITSCatalogue(
            output_name, open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite
        )

        n_hdu = len(output_dict.keys())
        if n_hdu != 5:
            raise IndexError(
                f"FITS output file data has {n_hdu} HDUs,"
                + " expected are 5"
            )
        for key in output_dict.keys():
            f.save_as_fits(output_dict[key], ext_name=key.upper())


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

        Returns
        -------
        dict
            Dictionary containing the NGMIX metacal results

        """
        tile_cat = file_io.FITSCatalogue(
            self._tile_cat_path,
            SEx_catalogue=True,
        )
        tile_cat.open()
        obj_id = np.copy(tile_cat.get_data()["NUMBER"])
        tile_vign = np.copy(tile_cat.get_data()["VIGNET"])
        tile_ra = np.copy(tile_cat.get_data()["XWIN_WORLD"])
        tile_dec = np.copy(tile_cat.get_data()["YWIN_WORLD"])
        tile_cat.close()

        f_wcs_file = SqliteDict(self._f_wcs_path)
        gal_vign_cat = SqliteDict(self._gal_vignet_path)
        bkg_vign_cat = SqliteDict(self._bkg_vignet_path)
        psf_vign_cat = SqliteDict(self._psf_vignet_path)
        weight_vign_cat = SqliteDict(self._weight_vignet_path)
        flag_vign_cat = SqliteDict(self._flag_vignet_path)

        final_res = []
        prior = self.get_prior()

        count = 0
        id_first = -1
        id_last = -1

        self._w_log.info(f"Processing objects # {self._id_obj_min} ... {self._id_obj_max}")
        for i_tile, id_tmp in enumerate(obj_id):

            if self._id_obj_min > 0 and id_tmp < self._id_obj_min:
                continue
            if self._id_obj_max > 0 and id_tmp > self._id_obj_max:
                continue

            if id_first == -1:
                id_first = id_tmp
            id_last = id_tmp
            str_id_tmp = str(id_tmp)

            gal_vign = []
            psf_vign = []
            sigma_psf = []
            weight_vign = []
            flag_vign = []
            jacob_list = []
            if psf_vign_cat[str_id_tmp] == "empty":
                self._w_log.info(f"Skipping object {id_tmp}: empty PSF vignet")
                continue

            self.check_key(str_id_tmp, gal_vign_cat, self._gal_vignet_path)
            self.check_key(str_id_tmp, bkg_vign_cat, self._bkg_vignet_path)
            self.check_key(str_id_tmp, flag_vign_cat, self._flag_vignet_path)
            self.check_key(str_id_tmp, weight_vign_cat, self._weight_vignet_path)
            if gal_vign_cat[str_id_tmp] == "empty":
                self._w_log.info(
                    f"Skipping object {id_tmp}: empty galaxy vignet"
                )
                continue

            psf_expccd_name = list(psf_vign_cat[str(id_tmp)].keys())
            for expccd_name_tmp in psf_expccd_name:

                self.check_key(expccd_name_tmp, gal_vign_cat[str_id_tmp], self._gal_vignet_path)
                self.check_key(expccd_name_tmp, bkg_vign_cat[str_id_tmp], self._bkg_vignet_path)
                self.check_key(expccd_name_tmp, flag_vign_cat[str_id_tmp], self._flag_vignet_path)
                self.check_key(expccd_name_tmp, weight_vign_cat[str_id_tmp], self._weight_vignet_path)

                exp_name, ccd_n = re.split("-", expccd_name_tmp)

                gal_vign_tmp = gal_vign_cat[str_id_tmp][expccd_name_tmp][
                    "VIGNET"
                ]
                if len(np.where(gal_vign_tmp.ravel() == 0)[0]) != 0:
                    self._w_log.info(
                        f"Skipping exp {expccd_name_tmp} for object {id_tmp}: zero-length galaxy vignet"
                    )
                    continue

                bkg_vign_tmp = bkg_vign_cat[str_id_tmp][expccd_name_tmp][
                    "VIGNET"
                ]
                gal_vign_sub_bkg = gal_vign_tmp - bkg_vign_tmp

                tile_vign_tmp = Ngmix.MegaCamFlip(
                    np.copy(tile_vign[i_tile]), int(ccd_n)
                )

                flag_vign_tmp = flag_vign_cat[str_id_tmp][expccd_name_tmp][
                    "VIGNET"
                ]
                flag_vign_tmp[np.where(tile_vign_tmp == -1e30)] = 2**10
                v_flag_tmp = flag_vign_tmp.ravel()
                if len(np.where(v_flag_tmp != 0)[0]) / (51 * 51) > 1 / 3.0:
                    self._w_log.info(
                        f"Skipping exp {expccd_name_tmp} for object {id_tmp}: mask > 1/3"
                    )
                    continue

                weight_vign_tmp = weight_vign_cat[str_id_tmp][expccd_name_tmp][
                    "VIGNET"
                ]

                jacob_tmp = get_jacob(
                    f_wcs_file[exp_name][int(ccd_n)]["WCS"],
                    tile_ra[i_tile],
                    tile_dec[i_tile],
                )

                header_tmp = fits.Header.fromstring(
                    f_wcs_file[exp_name][int(ccd_n)]["header"]
                )
                Fscale = header_tmp["FSCALE"]

                gal_vign_scaled = gal_vign_sub_bkg * Fscale
                weight_vign_scaled = weight_vign_tmp * 1 / Fscale**2

                gal_vign.append(gal_vign_scaled)
                psf_vign.append(
                    psf_vign_cat[str_id_tmp][expccd_name_tmp]["VIGNET"]
                )
                sigma_psf.append(
                    psf_vign_cat[str_id_tmp][expccd_name_tmp]["SHAPES"][
                        "SIGMA_PSF_HSM"
                    ]
                )
                weight_vign.append(weight_vign_scaled)
                flag_vign.append(flag_vign_tmp)
                jacob_list.append(jacob_tmp)

            if len(gal_vign) == 0:
                self._w_log.info(
                    f"Skipping object {id_tmp}: no exposure vignets added"
                )
                continue
            try:
                res = do_ngmix_metacal(
                    gal_vign,
                    psf_vign,
                    sigma_psf,
                    weight_vign,
                    flag_vign,
                    jacob_list,
                    prior,
                    self._pixel_scale,
                )
            except Exception as ee:
                self._w_log.info(
                    f"ngmix failed for object ID={id_tmp}.\nMessage: {ee}"
                )
                continue

            count = count + 1

            res["obj_id"] = id_tmp
            res["n_epoch_model"] = len(gal_vign)
            final_res.append(res)

        self._w_log.info(
            f"ngmix loop over objects finished, measured {count} "
            + f"objects, id first/last={id_first}/{id_last}"
        )

        f_wcs_file.close()
        gal_vign_cat.close()
        bkg_vign_cat.close()
        flag_vign_cat.close()
        weight_vign_cat.close()
        psf_vign_cat.close()

        # Put all results together
        res_dict = self.compile_results(final_res)

        # Save results
        self.save_results(res_dict)


def get_guess(
    img,
    pixel_scale,
    guess_flux_unit="img",
    guess_size_type="T",
    guess_size_unit="sky",
    guess_centroid=True,
    guess_centroid_unit="sky",
):
    r"""Get Guess.

    Get the guess vector for the NGMIX shape measurement
    ``[center_x, center_y, g1, g2, size_T, flux]``.
    No guesses are given for the ellipticity ``(0, 0)``.

    Parameters
    ----------
    img : numpy.ndarray
        Array containing the image
    pixel_scale : float
        Approximation of the pixel scale
    guess_flux_unit : str
        If ``img`` returns the flux in pixel units, otherwise if ``sky``
        returns the flux in :math:`{\rm arcsec}^{-2}`
    guess_size_type : str
        If ``T`` returns the size in quadrupole moments definition
        :math:`2\sigma^2`, otherwise if ``sigma`` returns the moments
        :math:`\sigma`
    guess_size_unit : str
        If ``img`` returns the size in pixel units, otherwise if ``sky``
        returns the size in arcsec
    guess_centroid : bool
        If ``True``, will return a guess on the object centroid, otherwise if
        ``False``, will return the image centre
    guess_centroid_unit : str
        If ``img`` returns the centroid in pixel unit, otherwise if ``sky``
        returns the centroid in arcsec

    Returns
    -------
    numpy.ndarray
        Return the guess array ``[center_x, center_y, g1, g2, size_T, flux]``

    Raises
    ------
    GalSimHSMError
        For an error in the computation of adaptive moments
    ValueError
        For invalid unit guess types

    """
    galsim_img = galsim.Image(img, scale=pixel_scale)

    hsm_shape = galsim.hsm.FindAdaptiveMom(galsim_img, strict=False)

    error_msg = hsm_shape.error_message

    if error_msg != "":
        raise galsim.hsm.GalSimHSMError(
            f"Error in adaptive moments :\n{error_msg}"
        )

    if guess_flux_unit == "img":
        guess_flux = hsm_shape.moments_amp
    elif guess_flux_unit == "sky":
        guess_flux = hsm_shape.moments_amp / pixel_scale**2
    else:
        raise ValueError(
            f"invalid guess_flux_unit '{guess_flux_unit}',"
            + " must be one of 'img', 'sky'"
        )

    if guess_size_unit == "img":
        size_unit = 1.0
    elif guess_size_unit == "sky":
        size_unit = pixel_scale
    else:
        raise ValueError(
            "invalid guess_size_unit '{guess_size_unit}',"
            + "must be one of 'img', 'sky'"
        )

    if guess_size_type == "sigma":
        guess_size = hsm_shape.moments_sigma * size_unit
    elif guess_size_type == "T":
        guess_size = 2 * (hsm_shape.moments_sigma * size_unit) ** 2

    if guess_centroid_unit == "img":
        centroid_unit = 1
    elif guess_centroid_unit == "sky":
        centroid_unit = pixel_scale
    else:
        raise ValueError(
            f"invalid guess_centroid_unit '{guess_centroid_unit}',"
            + "  must be one of 'img', 'sky'"
        )

    if guess_centroid:
        guess_centroid = (
            hsm_shape.moments_centroid - galsim_img.center
        ) * centroid_unit
    else:
        guess_centroid = galsim_img.center * centroid_unit

    guess = np.array(
        [guess_centroid.x, guess_centroid.y, 0.0, 0.0, guess_size, guess_flux]
    )

    return guess


def make_galsimfit(obs, model, guess0, prior=None, ntry=5):
    """Make GalSim Fit.

    Fit image using simple GalSim model.

    Parameters
    ----------
    obs : ngmix.observation.Observation
        Image to fit
    model : str
        Model for fit
    guess0 : numpy.ndarray
        Parameters of first model guess
    prior : ngmix.prior, optional
        Prior for fit paraemeters
    ntry : int, optional
        Number of tries for fit, the default is ``5``

    Returns
    -------
    dict
        Results

    Raises
    ------
    ngmix.BootGalFailure
        Failure to bootstrap galaxy

    """
    limit = 0.1

    guess = np.copy(guess0)
    fres = {}
    for it in range(ntry):
        guess[0:5] += urand(low=-limit, high=limit)
        guess[5:] *= 1 + urand(low=-limit, high=limit)
        fres["flags"] = 1
        try:
            fitter = ngmix.galsimfit.GalsimSimple(
                obs,
                model,
                prior=prior,
            )
            fitter.go(guess)
            fres = fitter.get_result()
        except Exception:
            continue

        if fres["flags"] == 0:
            break

    if fres["flags"] != 0:
        raise ngmix.gexceptions.BootGalFailure(
            "Failed to fit galaxy image with galsimfit"
        )

    fres["ntry"] = it + 1

    return fres


def get_jacob(wcs, ra, dec):
    """Get Jacobian.

    Return the Jacobian of the WCS at the required position.

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
    r"""Get Noise.

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
        cut = ``thresh`` * :math:`\sigma_{\rm noise}`;  the default is ``1.2``

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


def do_ngmix_metacal(
    gals, psfs, psfs_sigma, weights, flags, jacob_list, prior, pixel_scale
):
    """Do Ngmix Metacal.

    Perform the metacalibration on a multi-epoch object and return the joint
    shape measurement with NGMIX.

    Parameters
    ----------
    gals : list
        List of the galaxy vignets
    psfs : list
        List of the PSF vignets
    psfs_sigma : list
        List of the sigma PSFs
    weights : list
        List of the weight vignets
    flags : list
        List of the flag vignets
    jacob_list : list
        List of the Jacobians
    prior : ngmix.priors
        Priors for the fitting parameters
    pixel_scale : float
        pixel scale in arcsec

    Returns
    -------
    dict
        Dictionary containing the results of NGMIX metacal

    """
    n_epoch = len(gals)

    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    # Make observation
    gal_obs_list = ObsList()
    T_guess_psf = []
    psf_res_gT = {
        "g_PSFo": np.array([0.0, 0.0]),
        "g_err_PSFo": np.array([0.0, 0.0]),
        "T_PSFo": 0.0,
        "T_err_PSFo": 0.0,
    }
    gal_guess = []
    gal_guess_flag = True
    wsum = 0
    for n_e in range(n_epoch):

        psf_jacob = ngmix.Jacobian(
            row=(psfs[0].shape[0] - 1) / 2,
            col=(psfs[0].shape[1] - 1) / 2,
            wcs=jacob_list[n_e],
        )

        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)

        psf_T = psfs_sigma[n_e] * 1.17741 * pixel_scale

        weight_map = np.copy(weights[n_e])
        weight_map[np.where(flags[n_e] != 0)] = 0.0
        weight_map[weight_map != 0] = 1

        psf_guess = np.array([0.0, 0.0, 0.0, 0.0, psf_T, 1.0])
        try:
            psf_res = make_galsimfit(psf_obs, "gauss", psf_guess)
        except Exception:
            continue

        # Gal guess
        try:
            gal_guess_tmp = get_guess(
                gals[n_e], pixel_scale, guess_size_type="sigma"
            )
        except Exception:
            gal_guess_flag = False
            gal_guess_tmp = np.array([0.0, 0.0, 0.0, 0.0, 1, 100])

        # Recenter jacobian if necessary
        gal_jacob = ngmix.Jacobian(
            row=(gals[0].shape[0] - 1) / 2 + gal_guess_tmp[0],
            col=(gals[0].shape[1] - 1) / 2 + gal_guess_tmp[1],
            wcs=jacob_list[n_e],
        )

        # Noise handling
        if gal_guess_flag:
            sig_noise = get_noise(
                gals[n_e],
                weight_map,
                gal_guess_tmp,
                pixel_scale,
            )
        else:
            sig_noise = sigma_mad(gals[n_e])

        noise_img = np.random.randn(*gals[n_e].shape) * sig_noise
        noise_img_gal = np.random.randn(*gals[n_e].shape) * sig_noise

        gal_masked = np.copy(gals[n_e])
        if len(np.where(weight_map == 0)[0]) != 0:
            gal_masked[weight_map == 0] = noise_img_gal[weight_map == 0]

        weight_map *= 1 / sig_noise**2

        # Original PSF fit
        w_tmp = np.sum(weight_map)
        psf_res_gT["g_PSFo"] += psf_res["g"] * w_tmp
        psf_res_gT["g_err_PSFo"] += (
            np.array([psf_res["pars_err"][2], psf_res["pars_err"][3]]) * w_tmp
        )
        psf_res_gT["T_PSFo"] += psf_res["T"] * w_tmp
        psf_res_gT["T_err_PSFo"] += psf_res["T_err"] * w_tmp
        wsum += w_tmp

        gal_obs = Observation(
            gal_masked,
            weight=weight_map,
            jacobian=gal_jacob,
            psf=psf_obs,
            noise=noise_img,
        )

        if gal_guess_flag:
            gal_guess_tmp[:2] = 0
            gal_guess.append(gal_guess_tmp)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)
        gal_guess_flag = True

    if wsum == 0:
        raise ZeroDivisionError("Sum of weights = 0, division by zero")

    # Normalize PSF fit output
    for key in psf_res_gT.keys():
        psf_res_gT[key] /= wsum

    # Gal guess handling
    fail_get_guess = False
    if len(gal_guess) == 0:
        fail_get_guess = True
        gal_pars = [0.0, 0.0, 0.0, 0.0, 1, 100]
    else:
        gal_pars = np.mean(gal_guess, 0)

    psf_model = "gauss"
    gal_model = "gauss"

    # metacal specific parameters
    metacal_pars = {
        "types": ["noshear", "1p", "1m", "2p", "2m"],
        "step": 0.01,
        "psf": "gauss",
        "fixnoise": True,
        "cheatnoise": False,
        "symmetrize_psf": False,
        "use_noise_image": True,
    }

    Tguess = np.mean(T_guess_psf)

    # retry the fit twice
    ntry = 2

    obs_dict_mcal = ngmix.metacal.get_all_metacal(gal_obs_list, **metacal_pars)
    res = {"mcal_flags": 0}

    ntry = 5

    for key in sorted(obs_dict_mcal):

        fres = make_galsimfit(
            obs_dict_mcal[key], gal_model, gal_pars, prior=prior
        )

        res["mcal_flags"] |= fres["flags"]
        tres = {}

        for name in fres.keys():
            tres[name] = fres[name]
        tres["flags"] = fres["flags"]

        wsum = 0
        Tpsf_sum = 0
        gpsf_sum = np.zeros(2)
        npsf = 0
        for obs in obs_dict_mcal[key]:

            if hasattr(obs, "psf_nopix"):
                try:
                    psf_res = make_galsimfit(
                        obs.psf_nopix,
                        psf_model,
                        np.array([0.0, 0.0, 0.0, 0.0, Tguess, 1.0]),
                        ntry=ntry,
                    )
                except Exception:
                    continue
                g1, g2 = psf_res["g"]
                T = psf_res["T"]
            else:
                try:
                    psf_res = make_galsimfit(
                        obs.psf,
                        psf_model,
                        np.array([0.0, 0.0, 0.0, 0.0, Tguess, 1.0]),
                    )
                except Exception:
                    continue
                g1, g2 = psf_res["g"]
                T = psf_res["T"]

            # TODO we sometimes use other weights
            twsum = obs.weight.sum()

            wsum += twsum
            gpsf_sum[0] += g1 * twsum
            gpsf_sum[1] += g2 * twsum
            Tpsf_sum += T * twsum
            npsf += 1

        tres["gpsf"] = gpsf_sum / wsum
        tres["Tpsf"] = Tpsf_sum / wsum

        res[key] = tres

    # result dictionary, keyed by the types in metacal_pars above
    metacal_res = res

    metacal_res.update(psf_res_gT)
    metacal_res["moments_fail"] = fail_get_guess

    return metacal_res
