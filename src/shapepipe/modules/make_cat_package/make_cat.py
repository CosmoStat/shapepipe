"""MAKE CATALOGUE.

This module contains a class to create a shear catalogue.

:Authors: Axel Guinot, Martin Kilbinger

"""

import os
import re

import numpy as np
from astropy import coordinates as coords
from astropy import units as u
from astropy.wcs import WCS
from cs_util import size as cs_size
from sqlitedict import SqliteDict

from shapepipe.pipeline import file_io


def get_output_name(output_dir, file_number_string):
    """Get Output Name.

    Return output file name.

    Parameters
    ----------
    output_dir : str
        directory name
    file_number_string : str
        ShapePipe pipeline number string

    Returns
    -------
    str
        output path name

    """
    return f"{output_dir}/final_cat{file_number_string}.fits"


def prepare_final_cat_file(output_path, file_number_string):
    """Prepare Final Catalogue File.

    Create a ``FITSCatalogue`` object for the current file.

    Parameters
    ----------
    output_path : str
        Output file path
    file_number_string : str
        String with current file numbering

    Returns
    -------
    file_io.FITSCatalogue
        Output FITS file

    """

    output_name = get_output_name(output_path, file_number_string)

    return file_io.FITSCatalogue(
        output_name,
        open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
    )


def remove_field_name(arr, name):
    """Remove Field Name.

    Remove a column of a structured array from the given name.

    Parameters
    ----------
    arr : numpy.ndarray
        A numpy strucured array
    name : str
        Name of the field to remove

    Returns
    -------
    numpy.ndarray
        The structured array with the field removed

    """
    names = list(arr.dtype.names)
    if name in names:
        names.remove(name)
    arr2 = arr[names]
    return arr2


def save_sextractor_data(final_cat_file, sexcat_path, remove_vignet=True):
    """Save SExtractor Data.

    Save the SExtractor catalogue into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalogue
        Final catalogue
    sexcat_path : str
        Path to SExtractor catalogue to save
    remove_vignet : bool
        If ``True`` will not save the ``VIGNET`` field into the final catalogue

    Returns
    -------
    int
        Number of objects saved

    """
    sexcat_file = file_io.FITSCatalogue(sexcat_path, SEx_catalogue=True)
    sexcat_file.open()
    data = np.copy(sexcat_file.get_data())
    if remove_vignet:
        data = remove_field_name(data, "VIGNET")

    final_cat_file.save_as_fits(data, ext_name="RESULTS")

    cat_size = len(data)

    tile_id = float(
        ".".join(
            re.split("-", os.path.splitext(os.path.split(sexcat_path)[1])[0])[
                1:
            ]
        )
    )
    tile_id_array = np.ones(cat_size) * tile_id

    final_cat_file.open()
    final_cat_file.add_col("TILE_ID", tile_id_array)

    sexcat_file.close()

    return cat_size


def save_sm_data(
    final_cat_file,
    sexcat_sm_path,
    do_classif=True,
    star_thresh=0.003,
    gal_thresh=0.01,
    n_obj=-1,
):
    r"""Save Spread-Model Data.

    Save the spread-model data into the final catalogue.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalogue
        Final catalogue
    sexcat_sm_path : str
        Path to spread-model catalogue to save. If ``None``, spread_model is
        set to 99
    do_classif : bool
        If ``True`` objects will be classified into stars, galaxies, and other,
        using the classifier
        :math:`{\rm class} = {\rm sm} + 2 * {\rm sm}_{\rm err}`
    star_thresh : float
        Threshold for star selection; object is classified as star if
        :math:`|{\rm class}| <` ``star_thresh``
    gal_thresh : float
        Threshold for galaxy selection; object is classified as galaxy if
        :math:`{\rm class} >` ``gal_thresh``
    nobj : int, optional
        Number of objects, only used if sexcat_sm_path is ``-1``

    Returns
    -------
    int
        Number of objects saved
    """
    final_cat_file.open()

    if sexcat_sm_path is not None:
        sexcat_sm_file = file_io.FITSCatalogue(
            sexcat_sm_path,
            SEx_catalogue=True,
        )
        sexcat_sm_file.open()

        sm = np.copy(sexcat_sm_file.get_data()["SPREAD_MODEL"])
        sm_err = np.copy(sexcat_sm_file.get_data()["SPREADERR_MODEL"])

        sexcat_sm_file.close()

    else:
        sm = np.ones(n_obj) * 99
        sm_err = np.ones(n_obj) * 99

    final_cat_file.add_col("SPREAD_MODEL", sm)
    final_cat_file.add_col("SPREADERR_MODEL", sm_err)

    if do_classif:
        obj_flag = np.ones_like(sm, dtype="int16") * 2
        classif = sm + 2.0 * sm_err
        obj_flag[np.where(np.abs(classif) < star_thresh)] = 0
        obj_flag[np.where(classif > gal_thresh)] = 1

        final_cat_file.add_col("SPREAD_CLASS", obj_flag)

    final_cat_file.close()

    return n_obj


class SaveCatalogue:
    """Save Catalogue.

    Class to save catalogue.

    Parameters
    ----------
    final_cat_file : str
        Final catalogue file name
    cat_size_target : int
        target catalogue size
    w_log : logging.Logger
        Logging instance

    """

    def __init__(self, final_cat_file, cat_size_target, w_log):

        self._final_cat_file = final_cat_file
        self._cat_size_target = cat_size_target
        self._w_log = w_log

    def process(
        self,
        mode="",
        cat_path=None,
        moments=False,
    ):
        """Process Catalogue.

        Parameters
        ----------
        mode : str
            Run mode, options are ``ngmix``, ``galsim`` or ``psf``
        cat_path : str
            Path to input catalogue
        moments : bool
            Option to run ``ngmix`` mode with moments

        Returns
        --------
        str
            error message if failure; `None` if success

        """
        self._output_dict = {}

        self._final_cat_file.open()
        self._obj_id = np.copy(self._final_cat_file.get_data()["NUMBER"])

        err_msg = None
        if mode == "ngmix":
            err_msg = self._save_ngmix_data(cat_path, moments)
        elif mode == "galsim":
            self._save_galsim_shapes(cat_path)
        elif mode == "psf":
            self._save_psf_data(cat_path)
        else:
            err_msg = (
                f"Invalid process mode ({mode}) for "
                + '``make_cat.Savecatalogue``. Options are "ngmix", '
                + '"galsim" or "psf".'
            )

        if err_msg is None:

            for key in self._output_dict.keys():
                self._final_cat_file.add_col(key, self._output_dict[key])

        self._final_cat_file.close()

        return err_msg

    def _update_dict(self, key_string, value):
        """Update Dictionary.

        Update dictionary with value for all keys matching key string.

        Parameters
        ----------
        key_string : str
            Key string
        value : numpy.ndarray
            Value to be assigned to the keys

        """
        self._output_dict = {
            **self._output_dict,
            **{
                f"{key_string}{key_end}": np.copy(value)
                for key_end in self._key_ends
            },
        }

    def _add2dict(self, key, value, index=None):
        """Add to Dictionary.

        Add key, value pair to output dictionary.

        Parameters
        ----------
        key : str
            Dictionary key
        value : any
            Dictionary value
        index : int
            Dictionary element index

        """
        if not isinstance(index, type(None)):
            self._output_dict[key][index] = value
        else:
            self._output_dict[key] = value

    def _save_ngmix_data(self, ngmix_cat_path, moments=False):
        """Save NGMIX Data.

        Save the NGMIX catalogue into the final one.

        Column grammar: ``NGMIX[m]_<COMPONENT>[_ERR][_<OBJECT>]_<SHEAR>``,
        plus three OBJECT/SHEAR-less per-object metadata columns
        (``NGMIX[m]_MCAL_FLAGS``, ``NGMIX_N_EPOCH``,
        ``NGMIX_MCAL_TYPES_FAIL``). The galaxy is the implicit default object
        and carries NO ``OBJECT`` token (``NGMIX_G1_NOSHEAR``, dropping the
        ``GAL`` segment carried by the pre-#761 names). The explicit PSF
        objects are ``PSF_ORIG``
        (the original image PSF, fit by
        :func:`ngmix.average_original_psf`) and ``PSF_RECONV`` (the metacal
        reconvolution kernel, fit by :func:`ngmix.average_multiepoch_psf`);
        see those functions for what each PSF family IS. ``PSF_ORIG`` and
        ``PSF_RECONV`` are independent fits of different PSFs, no longer the
        single aliased value of the pre-fix code (shapepipe#749).

        Parameters
        ----------
        ngmix_cat_path : str
            Path to NGMIX catalogue
        moments : bool, optional
            If True, write the parallel ``NGMIXm_*`` (moments-branch) columns.

        """
        self._key_ends = ["1M", "1P", "2M", "2P", "NOSHEAR"]

        ngmix_cat_file = file_io.FITSCatalogue(ngmix_cat_path)
        ngmix_cat_file.open()

        ngmix_n_epoch = ngmix_cat_file.get_data()["n_epoch_model"]
        # Low number of ngmix objects could be due to
        # (1) shape measurement failures (e.g. missing PSF): ok, continue
        # (2) previous processing errors, e.g. premature run of
        # merge_sep_cats_runner: raise error
        if len(ngmix_n_epoch) / self._cat_size_target < 0.1:
            err_msg = (
                f"Merged shape catalogue {ngmix_cat_path} has very different"
                + f" size ({len(ngmix_n_epoch)}) compared to target size"
                + f" {self._cat_size_target})"
            )
            self._w_log.info(err_msg)
            #ngmix_cat_file.close()
            #return err_msg

        ngmix_mcal_types_fail = ngmix_cat_file.get_data()["mcal_types_fail"]
        # Needed in both moments and non-moments modes (used unconditionally
        # below), so read them outside the branch.
        ngmix_mcal_flags = ngmix_cat_file.get_data()["mcal_flags"]
        ngmix_id = ngmix_cat_file.get_data()["id"]

        n_obj = len(self._obj_id)
        self._w_log.info(f"writing ngmix info for {n_obj} objects")

        if moments:
            m = "m"
        else:
            m = ""

            self._add2dict("NGMIX_N_EPOCH", np.zeros(n_obj))
            self._add2dict("NGMIX_MCAL_TYPES_FAIL", np.zeros(n_obj))

        prefix = f"NGMIX{m}"

        # Galaxy (no object token), original image PSF (PSF_ORIG) and metacal
        # reconvolution kernel (PSF_RECONV); see ngmix.average_original_psf /
        # average_multiepoch_psf for what each PSF family is. G1/G2 are scalar
        # reduced-shear components, not a 2-vector. Sentinels:
        # sizes/fluxes/mags/flags 0, *_ERR fluxes/mags -1, ellipticities -10,
        # *_ERR sizes 1e30.
        for key_str in (
            f"{prefix}_T_",
            f"{prefix}_SNR_",
            f"{prefix}_FLUX_",
            f"{prefix}_MAG_",
            f"{prefix}_FLAGS_",
            f"{prefix}_T_PSF_ORIG_",
            f"{prefix}_T_PSF_RECONV_",
        ):
            self._update_dict(key_str, np.zeros(n_obj))
        for key_str in (
            f"{prefix}_FLUX_ERR_",
            f"{prefix}_MAG_ERR_",
        ):
            self._update_dict(key_str, np.ones(n_obj) * -1)
        for key_str in (
            f"{prefix}_G1_",
            f"{prefix}_G2_",
            f"{prefix}_G1_ERR_",
            f"{prefix}_G2_ERR_",
            f"{prefix}_G1_PSF_ORIG_",
            f"{prefix}_G2_PSF_ORIG_",
            f"{prefix}_G1_ERR_PSF_ORIG_",
            f"{prefix}_G2_ERR_PSF_ORIG_",
            f"{prefix}_G1_PSF_RECONV_",
            f"{prefix}_G2_PSF_RECONV_",
            f"{prefix}_G1_ERR_PSF_RECONV_",
            f"{prefix}_G2_ERR_PSF_RECONV_",
        ):
            self._update_dict(key_str, np.ones(n_obj) * -10.0)
        for key_str in (
            f"{prefix}_T_ERR_",
            f"{prefix}_T_ERR_PSF_ORIG_",
            f"{prefix}_T_ERR_PSF_RECONV_",
        ):
            self._update_dict(key_str, np.ones(n_obj) * 1e30)
        self._add2dict(f"{prefix}_MCAL_FLAGS", np.zeros(n_obj))

        for idx, id_tmp in enumerate(self._obj_id):
            ind = np.where(id_tmp == ngmix_id)[0]
            if len(ind) > 0:

                for key in self._key_ends:

                    ncf_data = ngmix_cat_file.get_data(key)

                    # Galaxy shape, size and photometry.
                    self._add2dict(
                        f"{prefix}_G1_{key}", ncf_data["g1"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_G2_{key}", ncf_data["g2"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_G1_ERR_{key}",
                        ncf_data["g1_err"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_G2_ERR_{key}",
                        ncf_data["g2_err"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_T_{key}", ncf_data["T"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_T_ERR_{key}",
                        ncf_data["T_err"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_SNR_{key}", ncf_data["s2n"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_FLUX_{key}",
                        ncf_data["flux"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_FLUX_ERR_{key}",
                        ncf_data["flux_err"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_MAG_{key}", ncf_data["mag"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_MAG_ERR_{key}",
                        ncf_data["mag_err"][ind[0]], idx
                    )
                    self._add2dict(
                        f"{prefix}_FLAGS_{key}",
                        ncf_data["flags"][ind[0]], idx
                    )

                    # Original image PSF (average_original_psf) and metacal
                    # reconvolution kernel (average_multiepoch_psf). Both PSF
                    # families share ONE write template, so the FITS column
                    # name (``{obj}``) and the res-key it reads (``{family}``)
                    # are generated from the same pair and cannot drift apart.
                    for family, obj in (
                        ("orig", "PSF_ORIG"),
                        ("reconv", "PSF_RECONV"),
                    ):
                        self._add2dict(
                            f"{prefix}_G1_{obj}_{key}",
                            ncf_data[f"g1_psf_{family}"][ind[0]], idx
                        )
                        self._add2dict(
                            f"{prefix}_G2_{obj}_{key}",
                            ncf_data[f"g2_psf_{family}"][ind[0]], idx
                        )
                        self._add2dict(
                            f"{prefix}_G1_ERR_{obj}_{key}",
                            ncf_data[f"g1_err_psf_{family}"][ind[0]], idx
                        )
                        self._add2dict(
                            f"{prefix}_G2_ERR_{obj}_{key}",
                            ncf_data[f"g2_err_psf_{family}"][ind[0]], idx
                        )
                        self._add2dict(
                            f"{prefix}_T_{obj}_{key}",
                            ncf_data[f"T_psf_{family}"][ind[0]], idx
                        )
                        self._add2dict(
                            f"{prefix}_T_ERR_{obj}_{key}",
                            ncf_data[f"T_err_psf_{family}"][ind[0]], idx
                        )

                self._add2dict(
                    f"NGMIX{m}_MCAL_FLAGS",
                    ngmix_mcal_flags[ind[0]],
                    idx,
                )

                if not moments:
                    self._add2dict(
                        f"NGMIX{m}_N_EPOCH",
                        ngmix_n_epoch[ind[0]],
                        idx,
                    )
                    self._add2dict(
                        f"NGMIX{m}_MCAL_TYPES_FAIL",
                        ngmix_mcal_types_fail[ind[0]],
                        idx,
                    )

        ngmix_cat_file.close()

        return None

    def _save_galsim_shapes(self, galsim_cat_path):
        """Save GalSim Shapes.

        Save the GalSim catalogue into the final one.

        Parameters
        ----------
        galsim_cat_path : str
            Path to GalSim catalogue to save

        """
        galsim_cat_file = file_io.FITSCatalogue(galsim_cat_path)
        galsim_cat_file.open()

        self._key_ends = galsim_cat_file.get_ext_name()[1:]

        galsim_id = galsim_cat_file.get_data()["id"]

        for key_str in (
            "GALSIM_T_",
            "GALSIM_T_PSF_",
            "GALSIM_FLUX_",
            "GALSIM_MAG_",
        ):
            self._update_dict(key_str, np.zeros(len(self._obj_id)))
        for key_str in ("GALSIM_FLUX_ERR_", "GALSIM_MAG_ERR_", "GALSIM_RES_"):
            self._update_dict(key_str, np.ones(len(self._obj_id)) * -1)
        for key_str in (
            "GALSIM_G1_",
            "GALSIM_G2_",
            "GALSIM_G1_ERR_",
            "GALSIM_G2_ERR_",
            "GALSIM_G1_UNCORR_",
            "GALSIM_G2_UNCORR_",
            "GALSIM_G1_PSF_",
            "GALSIM_G2_PSF_",
        ):
            self._update_dict(key_str, np.ones(len(self._obj_id)) * -10.0)
        self._update_dict(
            "GALSIM_FLAGS_",
            np.ones(len(self._obj_id), dtype="int16"),
        )

        for idx, id_tmp in enumerate(self._obj_id):
            ind = np.where(id_tmp == galsim_id)[0]
            if len(ind) > 0:

                for key in self._key_ends:

                    gcf_data = galsim_cat_file.get_data(key)

                    if key == "ORIGINAL_PSF":

                        # PSF columns sourced from the galaxy uncorr fields
                        # for this special extension (asymmetry preserved).
                        self._add2dict(
                            f"GALSIM_G1_PSF_{key}",
                            gcf_data["gal_uncorr_g1"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G2_PSF_{key}",
                            gcf_data["gal_uncorr_g2"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_T_PSF_{key}",
                            cs_size.sigma_to_T(gcf_data["gal_sigma"][ind[0]]),
                            idx
                        )

                    else:

                        self._add2dict(
                            f"GALSIM_G1_{key}", gcf_data["gal_g1"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G2_{key}", gcf_data["gal_g2"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G1_ERR_{key}",
                            gcf_data["gal_g1_err"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G2_ERR_{key}",
                            gcf_data["gal_g2_err"][ind[0]], idx
                        )

                        self._add2dict(
                            f"GALSIM_G1_UNCORR_{key}",
                            gcf_data["gal_uncorr_g1"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G2_UNCORR_{key}",
                            gcf_data["gal_uncorr_g2"][ind[0]], idx
                        )

                        self._add2dict(
                            f"GALSIM_T_{key}",
                            cs_size.sigma_to_T(gcf_data["gal_sigma"][ind[0]]),
                            idx
                        )

                        self._add2dict(
                            f"GALSIM_G1_PSF_{key}",
                            gcf_data["psf_g1"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_G2_PSF_{key}",
                            gcf_data["psf_g2"][ind[0]], idx
                        )
                        self._add2dict(
                            f"GALSIM_T_PSF_{key}",
                            cs_size.sigma_to_T(gcf_data["psf_sigma"][ind[0]]),
                            idx
                        )

                        flux = gcf_data["gal_flux"][ind[0]]
                        flux_err = gcf_data["gal_flux_err"][ind[0]]
                        self._add2dict(f"GALSIM_FLUX_{key}", flux, idx)
                        self._add2dict(f"GALSIM_FLUX_ERR_{key}", flux_err, idx)

                        mag = gcf_data["gal_mag"][ind[0]]
                        mag_err = gcf_data["gal_mag_err"][ind[0]]
                        self._add2dict(f"GALSIM_MAG_{key}", mag, idx)
                        self._add2dict(f"GALSIM_MAG_ERR_{key}", mag_err, idx)

                        flags = gcf_data["gal_flag"][ind[0]]
                        self._add2dict(f"GALSIM_FLAGS_{key}", flags, idx)

                        res = gcf_data["gal_resolution"][ind[0]]
                        self._add2dict(f"GALSIM_RES_{key}", res, idx)

        galsim_cat_file.close()

    def _save_psf_data(self, galaxy_psf_path):
        """Save PSF data.

        Save the PSF catalogue into the final one.

        Parameters
        ----------
        galaxy_psf_path : str
            Path to the PSF catalogue to save

        """
        galaxy_psf_cat = SqliteDict(galaxy_psf_path)

        max_epoch = np.max(self._final_cat_file.get_data()["N_EPOCH"]) + 1

        self._output_dict = {
            f"HSM_G1_PSF_{idx + 1}": np.ones(len(self._obj_id)) * -10.0
            for idx in range(max_epoch)
        }
        self._output_dict = {
            **self._output_dict,
            **{
                f"HSM_G2_PSF_{idx + 1}": np.ones(len(self._obj_id)) * -10.0
                for idx in range(max_epoch)
            },
        }
        self._output_dict = {
            **self._output_dict,
            **{
                f"HSM_T_PSF_{idx + 1}": np.zeros(len(self._obj_id))
                for idx in range(max_epoch)
            },
        }
        self._output_dict = {
            **self._output_dict,
            **{
                f"HSM_FLAG_PSF_{idx + 1}": np.ones(
                    len(self._obj_id), dtype="int16"
                )
                for idx in range(max_epoch)
            },
        }

        for idx, id_tmp in enumerate(self._obj_id):

            if galaxy_psf_cat[str(id_tmp)] == "empty":
                continue

            for epoch, key in enumerate(galaxy_psf_cat[str(id_tmp)].keys()):

                gpc_data = galaxy_psf_cat[str(id_tmp)][key]

                if gpc_data["SHAPES"]["HSM_FLAG_PSF"] != 0:
                    continue

                self._add2dict(
                    f"HSM_G1_PSF_{epoch + 1}",
                    gpc_data["SHAPES"]["HSM_G1_PSF"], idx
                )
                self._add2dict(
                    f"HSM_G2_PSF_{epoch + 1}",
                    gpc_data["SHAPES"]["HSM_G2_PSF"], idx
                )

                # HSM_T_PSF already holds T (sigma_to_T applied at the
                # producer's _interpolate_me); read straight through.
                self._add2dict(
                    f"HSM_T_PSF_{epoch + 1}",
                    gpc_data["SHAPES"]["HSM_T_PSF"], idx
                )

                self._add2dict(
                    f"HSM_FLAG_PSF_{epoch + 1}",
                    gpc_data["SHAPES"]["HSM_FLAG_PSF"], idx
                )

        galaxy_psf_cat.close()
