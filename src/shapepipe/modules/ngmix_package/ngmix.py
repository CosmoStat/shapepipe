"""NGMIX.

This module contains a class for ngmix shape measurement.

:Authors: Lucie Baumont, Axel Guinot

"""

import os
import re
import ngmix
import galsim
import numpy as np
from shutil import copyfile
from astropy.io import fits
from modopt.math.stats import sigma_mad
from ngmix.observation import Observation, ObsList
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io


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


# I still don't know how to handle this
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
        self.size = np.copy(data['FWHM_WORLD']) if 'FWHM_WORLD' in cols else None
        self.e = np.copy(data['ELLIPTICITY']) if 'ELLIPTICITY' in cols else None
        self.theta = np.copy(data['THETA_WIN_WORLD']) if 'THETA_WIN_WORLD' in cols else None

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
        self.jacobs = []
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
    """
    def __init__(
        self,
        gal_vignet_path,
        bkg_vignet_path,
        psf_vignet_path,
        weight_vignet_path,
        flag_vignet_path,
        f_wcs_path

    ):
        self.f_wcs_file = SqliteDict(f_wcs_path)
        self.gal_vign_cat = SqliteDict(gal_vignet_path)
        self.bkg_vign_cat = SqliteDict(bkg_vignet_path)
        self.psf_vign_cat = SqliteDict(psf_vignet_path)
        self.weight_vign_cat = SqliteDict(weight_vignet_path)
        self.flag_vign_cat = SqliteDict(flag_vignet_path)

    def close(self):
        self.f_wcs_file.close()
        self.gal_vign_cat.close()
        self.bkg_vign_cat.close()
        self.flag_vign_cat.close()
        self.weight_vign_cat.close()
        self.psf_vign_cat.close()

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
    ):

        if len(input_file_list) != 6:
            raise IndexError(
                f"Input file list has length {len(input_file_list)},"
                + " required is 6"
            )

        self._tile_cat_path = input_file_list[0]
        self._vignet_cat = Vignet(
            input_file_list[1],
            input_file_list[2],
            input_file_list[3],
            input_file_list[4],
            input_file_list[5],
            f_wcs_path
        )
        #self._gal_vignet_path = input_file_list[1]
        #self._bkg_vignet_path = input_file_list[2]
        #self._psf_vignet_path = input_file_list[3]
        #self._weight_vignet_path = input_file_list[4]
        #self._flag_vignet_path = input_file_list[5]

      

        self._output_dir = output_dir
        self._file_number_string = file_number_string

        self._zero_point = zero_point
        self._pixel_scale = pixel_scale

        self._f_wcs_path = f_wcs_path

        self._save_batch = save_batch
        self._id_obj_min = id_obj_min
        self._id_obj_max = id_obj_max

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
        names = ["1m", "1p", "2m", "2p", "noshear"]
        names2 = [
            'id',
            'n_epoch_model',
            'moments_fail',
            'ntry_fit',
            'g1_psfo_ngmix',
            'g2_psfo_ngmix',
            'T_psfo_ngmix',
            'T_err_psfo_ngmix',
            'r50_psfo_ngmix',
            'g1_err_psfo_ngmix',
            'g2_err_psfo_ngmix',
            'r50_err_psfo_ngmix',
            'g1',
            'g1_err',
            'g2',
            'g2_err',
            'T',
            'T_err',
            'Tpsf',
            'r50',
            'r50_err',
            'r50psf',
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
                output_dict[name]["ntry_fit"].append(results[idx][name]["nfev"])
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
                output_dict[name]['r50_psfo_ngmix'].append(
                    results[idx]['r50_PSFo']
                )
                output_dict[name]['r50_err_psfo_ngmix'].append(
                    results[idx]['r50_err_PSFo']
                )
                output_dict[name]["T"].append(results[idx][name]["T"])
                output_dict[name]["T_err"].append(results[idx][name]["T_err"])
                output_dict[name]["Tpsf"].append(results[idx]["T_PSFo"])
                output_dict[name]["g1_psf"].append(results[idx]["g_PSFo"][0])
                output_dict[name]["g2_psf"].append(results[idx]["g_PSFo"][1])
                output_dict[name]['r50'].append(results[idx][name]['pars'][4])
                output_dict[name]['r50_err'].append(results[idx][name]['pars_err'][4])
                output_dict[name]['r50psf'].append(results[idx]["r50_PSFo"])
                output_dict[name]["g1"].append(results[idx][name]["g"][0])
                output_dict[name]["g2"].append(results[idx][name]["g"][1])
                output_dict[name]["g1_err"].append(
                    np.sqrt(results[idx][name]["g_cov"][0, 0])
                )
                output_dict[name]["g2_err"].append(
                    np.sqrt(results[idx][name]["g_cov"][1, 1])
                )
                output_dict[name]["flux"].append(results[idx][name]["flux"])
                output_dict[name]["flux_err"].append(results[idx][name]["flux_err"])
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
                            "HDU extension {ext_name} from existing FITS file"
                            + f" not found in data"
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
                        print("output_dict", [col for col in new_data])
                        print("existing_da", existing_data.dtype.names)
                        raise ValueError(
                            "Mismatch between existing columns and new data columns."
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
                )
            except Exception as ee:
                self._w_log.info(
                    f'ngmix failed for object ID={obj_id}.\nMessage: {ee}'
                )
                n_ngmix_fail += 1
                continue

            res['obj_id'] = obj_id
            res['n_epoch_model'] = len(stamp.gals)
            res['moments_fail'] = sum(
                1 for k in ['noshear', '1p', '1m', '2p', '2m']
                if res.get(k, {}).get('flags', 0) != 0
            )
            r50_psfo = np.sqrt(max(psf_res['T_psf'], 0) / 2)
            res['g_PSFo'] = psf_res['g_psf']
            res['g_err_PSFo'] = psf_res['g_psf_err']
            res['T_PSFo'] = psf_res['T_psf']
            res['T_err_PSFo'] = psf_res['T_psf_err']
            res['r50_PSFo'] = r50_psfo
            res['r50_err_PSFo'] = (
                psf_res['T_psf_err'] / (2 * r50_psfo) if r50_psfo > 0 else np.nan
            )
            res['mcal_flags'] = 0
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

        weight_vign = (
            vignet.weight_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        )

        jacob = get_galsim_jacobian(
            vignet.f_wcs_file[exp_name][int(ccd_n)]['WCS'],
            tile_cat.ra[i_tile],
            tile_cat.dec[i_tile]
        )

        header = fits.Header.fromstring(
            vignet.f_wcs_file[exp_name][int(ccd_n)]['header']
        )

        # rescale by relative zero-points
        gal_vign_scaled, weight_vign_scaled = rescale_epoch_fluxes(
            gal_vign_sub_bkg,
            weight_vign,
            header
            )

        # gather postage stamps in all of the epochs
        stamp.gals.append(gal_vign_scaled)
        stamp.psfs.append(
            vignet.psf_vign_cat[str(obj_id)][expccd_name]['VIGNET']
        )
        stamp.weights.append(weight_vign_scaled)
        stamp.flags.append(flag_vign)
        stamp.jacobs.append(jacob)
                
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

def rescale_epoch_fluxes(gal,weight,header):
    """rescale epochs by relative zeropoints to be on the same flux scale
        
    Parameters
    ----------
    gal : numpy.ndarray
        background subtracted galaxy image
    weight : numpy.ndarray
        weight image
    header : 
        image header
        
    Returns
    -------
    numpy.ndarray
        rescaled galaxy image
    numpy.ndarray
        rescaled weight image
    """
    Fscale = header['FSCALE']

    gal_scaled = gal * Fscale
    weight_scaled = weight * 1 / Fscale ** 2

    return gal_scaled, weight_scaled

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
        ``[x0, y0, g1, g2, r50, flux]``
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

def prepare_ngmix_weights(gal, weight, flag, rng):
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

    Returns
    -------
    numpy.ndarray
        Galaxy image with masked pixels replaced by noise.
    numpy.ndarray
        Variance map for NGMIX.
    numpy.ndarray
        Noise image.
    """
    weight_map = np.copy(weight)
    weight_map[flag != 0] = 0.0

    sig_noise = sigma_mad(gal)

    noise_img = rng.standard_normal(gal.shape) * sig_noise
    noise_img_gal = rng.standard_normal(gal.shape) * sig_noise

    gal_masked = np.copy(gal)
    if (weight_map == 0).any():
        gal_masked[weight_map == 0] = noise_img_gal[weight_map == 0]

    weight_map *= 1 / sig_noise ** 2

    return gal_masked, weight_map, noise_img

def make_ngmix_observation(gal, weight, flag, psf, wcs, rng):
    """Build an ngmix Observation for a single galaxy epoch.

    The galaxy Jacobian is re-centered on the HSM centroid so that the
    centroid prior (centered at the Jacobian origin) does not bias the fit.

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
        gal, weight, flag, rng
    )

    # Re-center Jacobian on HSM centroid (pixel offset from stamp center).
    # Fixes: centroid prior biases fit when galaxy is offset from stamp center.
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
        Keys: 'g_psf', 'g_psf_err', 'T_psf', 'T_psf_err' (weighted averages).
    """
    # create dictionary
    names = ['T_psf', 'T_psf_err', 'g_psf', 'g_psf_err']
    psf_dict = {k: [] for k in names}
    nepoch = len(obsdict['noshear'])
    wsum = 0
    g_psf_sum = np.array([0., 0.])
    g_psf_err_sum = np.array([0., 0.])
    T_psf_sum = 0
    T_psf_err_sum = 0
    for n_e in np.arange(nepoch):
        T_psf=obsdict['noshear'][n_e].psf.meta['result']['T']
        T_psf_err=obsdict['noshear'][n_e].psf.meta['result']['T_err']
        g_psf=obsdict['noshear'][n_e].psf.meta['result']['g']
        g_psf_err=obsdict['noshear'][n_e].psf.meta['result']['g_err']
        ne_wsum = obsdict['noshear'][n_e].weight.sum()

        wsum += ne_wsum
        g_psf_sum += g_psf * ne_wsum
        g_psf_err_sum += g_psf_err * ne_wsum
        T_psf_sum += T_psf * ne_wsum
        T_psf_err_sum += T_psf_err * ne_wsum

    if wsum == 0:
        raise ZeroDivisionError('Sum of weights = 0, division by zero')

    psf_dict['g_psf'] = g_psf_sum / wsum
    psf_dict['g_psf_err'] = g_psf_err_sum / wsum
    psf_dict['T_psf'] = T_psf_sum / wsum
    psf_dict['T_psf_err'] = T_psf_err_sum / wsum

    return psf_dict


def do_ngmix_metacal(stamp, prior, flux_guess, rng):
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

    Returns
    -------
    tuple
        (resdict, psf_res) where resdict is the MetacalBootstrapper result
        dict and psf_res is the averaged PSF dict from average_multiepoch_psf.
    """
    n_epoch = len(stamp.gals)
    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    psf_model = 'gauss'
    gal_model = 'gauss'

    gal_obs_list = ObsList()
    for n_e in range(n_epoch):
        gal_obs = make_ngmix_observation(
            stamp.gals[n_e],
            stamp.weights[n_e],
            stamp.flags[n_e],
            stamp.psfs[n_e],
            stamp.jacobs[n_e],
            rng,
        )
        gal_obs_list.append(gal_obs)

    fitter = ngmix.fitting.Fitter(model=gal_model, prior=prior)
    guesser = ngmix.guessers.TPSFFluxAndPriorGuesser(rng=rng, T=0.25, prior=prior)

    psf_fitter = ngmix.fitting.Fitter(model=psf_model, prior=prior)
    psf_guesser = ngmix.guessers.TFluxGuesser(rng=rng, T=0.25, prior=prior, flux=flux_guess)

    psf_runner = ngmix.runners.PSFRunner(fitter=psf_fitter, guesser=psf_guesser, ntry=2)
    runner = ngmix.runners.Runner(fitter=fitter, guesser=guesser, ntry=5)

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
