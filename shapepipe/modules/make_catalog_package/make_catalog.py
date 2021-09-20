# -*- coding: utf-8 -*-

"""MAKE CATALOG

This module contains a class to create a shear catalog.

:Author: Axel Guinot

"""

from astropy import units as u
from astropy import coordinates as coords
from astropy.wcs import WCS
from shapepipe.pipeline import file_io
from sqlitedict import SqliteDict

import os
import re
import numpy as np


def prepare_final_cat_file(output_path, file_number_string):
    """Prepare Final Catalog File.

    Create a FITSCatalog object for the current file.

    Parameters
    ----------
    output_path : str
        Output file path
    file_number_string : str
        String with current file numbering

    Returns
    -------
    shapepipe.pipeline.file_io.FITSCatalog
        Output FITS file

    """
    output_name = (
        f'{output_path}/final_cat{file_number_string}.fits'
    )

    return file_io.FITSCatalog(
        output_name,
        open_mode=file_io.BaseCatalog.OpenMode.ReadWrite,
    )


def remove_field_name(arr, name):
    """Remove Field Name.

    Remove a column of a structured array from the given name.

    Parameters
    ----------
    arr : numpy.ndarray
        A numpy strucured array.
    name : str
        Name of the field to remove.

    Returns
    -------
    numpy.ndarray
        The structured with the field removed.

    """
    names = list(arr.dtype.names)
    if name in names:
        names.remove(name)
    arr2 = arr[names]
    return arr2


def save_sextractor_data(final_cat_file, sexcat_path, remove_vignet=True):
    """Save SExtractor Data.

    Save the SExtractor catalog into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    sexcat_path : str
        Path to SExtractor catalog to save.
    remove_vignet : bool
        If True will not save the 'VIGNET' field into the final catalog.

    """
    sexcat_file = file_io.FITSCatalog(sexcat_path, SEx_catalog=True)
    sexcat_file.open()
    data = np.copy(sexcat_file.get_data())
    if remove_vignet:
        data = remove_field_name(data, 'VIGNET')

    final_cat_file.save_as_fits(data, ext_name='RESULTS')

    cat_size = len(data)

    tile_id = float('.'.join(
        re.split('-', os.path.splitext(os.path.split(sexcat_path)[1])[0])[1:]
    ))
    tile_id_array = np.ones(cat_size) * tile_id

    final_cat_file.open()
    final_cat_file.add_col('TILE_ID', tile_id_array)

    sexcat_file.close()


def save_sm_data(
    final_cat_file,
    sexcat_sm_path,
    do_classif=True,
    star_thresh=0.003,
    gal_thresh=0.01,
):
    """Save Spread-Model Data.

    Save the spread-model data into the final catalog.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    sexcat_sm_path : str
        Path to spread-model catalog to save.
    do_classif : bool
        If True will make a star/galaxy classification. Based on :
        class = sm + 5/3 * sm_err
    star_thresh : float
        Threshold for star selection. |class| < star_thresh
    gal_thresh : float
        Threshold for galaxy selection. class > gal_thresh

    """

    final_cat_file.open()

    sexcat_sm_file = file_io.FITSCatalog(sexcat_sm_path, SEx_catalog=True)
    sexcat_sm_file.open()

    sm = np.copy(sexcat_sm_file.get_data()['SPREAD_MODEL'])
    sm_err = np.copy(sexcat_sm_file.get_data()['SPREADERR_MODEL'])

    sexcat_sm_file.close()

    final_cat_file.add_col('SPREAD_MODEL', sm)
    final_cat_file.add_col('SPREADERR_MODEL', sm_err)

    if do_classif:
        obj_flag = np.ones_like(sm, dtype='int16') * 2
        classif = sm + 2. * sm_err
        obj_flag[np.where(np.abs(classif) < star_thresh)] = 0
        obj_flag[np.where(classif > gal_thresh)] = 1

        final_cat_file.add_col('SPREAD_CLASS', obj_flag)

    final_cat_file.close()


def get_ra_dec(xxx, yyy):
    """Get RA Dec.

    Transform Stephen Gwyn naming into RA and Dec position

    Parameters
    ----------
    xxx : int
        First 3 numbers in the tile name.
    yyy : int
        Last 3 numbers in the tile name.

    Returns
    -------
    ra, dec : tuple
        Ra and dec positions for the center of tile.

    """

    return xxx / 2 / np.cos((yyy / 2 - 90) * np.pi / 180), yyy / 2 - 90


def get_tile_wcs(xxx, yyy):
    """Get Tile WCS.

    Create an astropy.wcs.WCS object from the name of the tile.

    Parameters
    ----------
    xxx : int
        First 3 numbers in the tile name.
    yyy : int
        Last 3 numbers in the tile name.

    Returns
    -------
    w : astropy.wcs.WCS
        WCS for the tile.

    """
    ra, dec = get_ra_dec(xxx, yyy)

    w = WCS(naxis=2)
    w.wcs.crval = np.array([ra, dec])
    w.wcs.crpix = np.array([5000, 5000])
    w.wcs.cd = np.array(
        [[0.187 / 3600, 0], [0, 0.187 / 3600]]
    )
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cunit = ['deg', 'deg']
    w._naxis = [10000, 10000]

    return w


def remove_common_elements(
    final_cat_file,
    tiles_id_file,
    pos_param=['XWIN_WORLD', 'YWIN_WORLD']
):
    """Remove Common Elements.

    This function create a mask for the overlapping region between 2 stack
    images.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    tile_id_file : str
        Path to the file containing all the tile IDs.
    pos_param : list
        List with the column name for ra and dec positions.

    """
    key_id = 'TILE_ID'

    final_cat_file.open()
    ra_tile = final_cat_file.get_data()[pos_param[0]]
    dec_tile = final_cat_file.get_data()[pos_param[1]]
    tile_id = final_cat_file.get_data()[key_id][0] * 1000

    all_id = np.loadtxt(tiles_id_file, unpack=True) * 1000

    all_tile_ra, all_tile_dec = get_ra_dec(
        (all_id / 1000).astype(int),
        all_id - (all_id / 1000).astype(int) * 1000,
    )
    all_tile_coord = coords.SkyCoord(
        ra=all_tile_ra * u.deg,
        dec=all_tile_dec * u.deg,
    )

    ra_m_tile, dec_m_tile = get_ra_dec(
        int(tile_id / 1000),
        tile_id - int(tile_id / 1000) * 1000,
    )

    mask_tile = np.ones_like(ra_tile, dtype=bool)

    tile_main_coord = coords.SkyCoord(
        ra=ra_m_tile * u.deg,
        dec=dec_m_tile * u.deg,
    )
    catalog_coord = coords.SkyCoord(ra=ra_tile * u.deg, dec=dec_tile * u.deg)

    sep_ref = tile_main_coord.separation(catalog_coord)

    close_tiles = np.argsort(tile_main_coord.separation(all_tile_coord))[1:9]

    for id_tile_check in all_id[close_tiles]:

        ra_m_check, dec_m_check = get_ra_dec(
            int(id_tile_check / 1000),
            id_tile_check - int(id_tile_check / 1000) * 1000,
        )

        tile_check_coord = coords.SkyCoord(
            ra=ra_m_check * u.deg,
            dec=dec_m_check * u.deg,
        )

        sep_check = tile_check_coord.separation(catalog_coord)

        w_tile_check = get_tile_wcs(
            int(id_tile_check / 1000),
            id_tile_check - int(id_tile_check / 1000) * 1000,
        )

        mask_tile &= (
            np.less(sep_ref, sep_check)
            | np.invert(w_tile_check.footprint_contains(catalog_coord))
        )

    final_cat_file.add_col('FLAG_TILING', mask_tile)

    final_cat_file.close()


class SaveCatalog:
    """Save Catalog.

    Class to save catalog.

    Parameters
    ----------
    final_cat_file : str
        Final catalogue file name

    """

    def __init__(self, final_cat_file):

        self.final_cat_file = final_cat_file

    def process(
        self,
        mode='',
        cat_path=None,
        moments=False,
    ):
        """Process Catalog.

        Parameters
        ----------
        mode : str
            Run mode, options are ('ngmix', 'galsim', 'psf')
        cat_path : str
            Path to input catalog
        moments : bool
            Opotion to run ngmix mode with moment

        """
        self._output_dict = {}

        self.final_cat_file.open()
        self._obj_id = np.copy(self.final_cat_file.get_data()['NUMBER'])

        if mode == 'ngmix':
            self._save_ngmix_data(cat_path, moments)
        elif mode == 'galsim':
            self._save_galsim_shapes(cat_path)
        elif mode == 'psf':
            self._save_psf_data(cat_path)
        else:
            raise ValueError(
                f'Invalid process mode ({mode}) for '
                + '``make_catalog.SaveCatalog``. Options are "ngmix", '
                + '"galsim" or "psf".'
            )

        for key in self._output_dict.keys():
            self.final_cat_file.add_col(key, self._output_dict[key])

        self.final_cat_file.close()

    def _update_dict(self, key_string, value):
        """Update Dictionary.

        Update dictionary with value for all keys matching key string.

        Parameters
        ----------
        key_string : str
            Key string
        value : np.ndarray
            Value to be assigned to the keys

        """
        self._output_dict = {
            **self._output_dict,
            **{
                f'{key_string}{key_end}': np.copy(value)
                for key_end in self._key_ends
            }
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
        """Save ngmix Data.

        Save the ngmix catalog into the final one.

        Parameters
        ----------
        ngmix_cat_path : str
            Path to ngmix catalog

        """
        self._key_ends = ['1M', '1P', '2M', '2P', 'NOSHEAR']

        ngmix_cat_file = file_io.FITSCatalog(ngmix_cat_path)
        ngmix_cat_file.open()

        ngmix_n_epoch = ngmix_cat_file.get_data()['n_epoch_model']
        ngmix_mom_fail = ngmix_cat_file.get_data()['moments_fail']

        if moments:
            m = 'm'
        else:
            m = ''

            ngmix_mcal_flags = ngmix_cat_file.get_data()['mcal_flags']
            ngmix_id = ngmix_cat_file.get_data()['id']

            self._add2dict('NGMIX_N_EPOCH', np.zeros(len(self._obj_id)))
            self._add2dict('NGMIX_MOM_FAIL', np.zeros(len(self._obj_id)))

        prefix = f'NGMIX{m}'

        for key_str in (
            f'{prefix}_T_',
            f'{prefix}_Tpsf_',
            f'{prefix}_SNR_',
            f'{prefix}_FLUX_',
            f'{prefix}_MAG_',
            f'{prefix}_FLAGS_',
            f'{prefix}_T_PSFo_',
        ):
            self._update_dict(key_str, np.zeros(len(self._obj_id)))
        for key_str in (f'NGMIX{m}_FLUX_ERR_', f'NGMIX{m}_MAG_ERR_'):
            self._update_dict(key_str, np.ones(len(self._obj_id)) * -1)
        for key_str in (
            f'NGMIX{m}_ELL_',
            f'NGMIX{m}_ELL_ERR_',
            f'NGMIX{m}_ELL_PSFo_'
        ):
            self._update_dict(key_str, np.ones((len(self._obj_id), 2)) * -10.0)
        self._update_dict(
            f'NGMIX{m}_T_ERR_',
            np.ones(len(self._obj_id)) * 1e30,
        )
        self._add2dict(f'NGMIX{m}_MCAL_FLAGS', np.zeros(len(self._obj_id)))

        for idx, _ in enumerate(self._obj_id):
            for key in self._key_ends:
                x = self._output_dict[f'NGMIX{m}_ELL_{key}'][idx]
                if np.all(x != np.array([-10., -10.])):
                    print(x)

        for idx, id_tmp in enumerate(self._obj_id):
            ind = np.where(id_tmp == ngmix_id)[0]
            if len(ind) > 0:

                for key in self._key_ends:

                    ncf_data = ngmix_cat_file.get_data(key)

                    g = (ncf_data['g1'][ind[0]], ncf_data['g2'][ind[0]])
                    g_err = (
                        ncf_data['g1_err'][ind[0]],
                        ncf_data['g2_err'][ind[0]]
                    )

                    self._add2dict(f'NGMIX{m}_ELL_{key}', g, idx)
                    self._add2dict(f'NGMIX{m}_ELL_ERR_{key}', g_err, idx)

                    t = ncf_data['T'][ind[0]]
                    t_err = ncf_data['T_err'][ind[0]]
                    tpsf = ncf_data['Tpsf'][ind[0]]
                    self._add2dict(f'NGMIX{m}_T_{key}', t, idx)
                    self._add2dict(f'NGMIX{m}_T_ERR_{key}', t_err, idx)
                    self._add2dict(f'NGMIX{m}_Tpsf_{key}', tpsf, idx)

                    s2n = ncf_data['s2n'][ind[0]]
                    self._add2dict(f'NGMIX{m}_SNR_{key}', s2n, idx)

                    flux = ncf_data['flux'][ind[0]]
                    flux_err = ncf_data['flux_err'][ind[0]]
                    self._add2dict(f'NGMIX{m}_FLUX_{key}', flux, idx)
                    self._add2dict(f'NGMIX{m}_FLUX_ERR_{key}', flux_err, idx)

                    mag = ncf_data['mag'][ind[0]]
                    mag_err = ncf_data['mag_err'][ind[0]]
                    self._add2dict(f'NGMIX{m}_MAG_{key}', mag, idx)
                    self._add2dict(f'NGMIX{m}_MAG_ERR_{key}', mag_err, idx)

                    flags = ncf_data['flags'][ind[0]]
                    self._add2dict(f'NGMIX{m}_FLAGS_{key}', flags, idx)

                    g_psf = (
                        ncf_data['g1_psfo_ngmix'][ind[0]],
                        ncf_data['g2_psfo_ngmix'][ind[0]]
                    )
                    self._add2dict(f'NGMIX{m}_ELL_PSFo_{key}', g_psf, idx)

                    t_psfo = ncf_data['T_psfo_ngmix'][ind[0]]
                    self._add2dict(f'NGMIX{m}_T_PSFo_{key}', t_psfo, idx)

                self._add2dict(
                    f'NGMIX{m}_MCAL_FLAGS',
                    ngmix_mcal_flags[ind[0]],
                    idx,
                )

                if not moments:
                    self._add2dict(
                        f'NGMIX{m}_N_EPOCH',
                        ngmix_n_epoch[ind[0]],
                        idx,
                    )
                    self._add2dict(
                        f'NGMIX{m}_MOM_FAIL',
                        ngmix_mom_fail[ind[0]],
                        idx,
                    )

        ngmix_cat_file.close()

    def _save_galsim_shapes(self, galsim_cat_path):
        """Save GalSim Shapes.

        Save the GalSim catalog into the final one.

        Parameters
        ----------
        ngmix_cat_path : str
            Path to ngmix catalog to save.

        """
        galsim_cat_file = file_io.FITSCatalog(galsim_cat_path)
        galsim_cat_file.open()

        self._key_ends = galsim_cat_file.get_ext_name()[1:]

        galsim_id = galsim_cat_file.get_data()['id']

        for key_str in (
            'GALSIM_GAL_SIGMA_',
            'GALSIM_PSF_SIGMA_',
            'GALSIM_FLUX_',
            'GALSIM_MAG_',
        ):
            self._update_dict(key_str, np.zeros(len(self._obj_id)))
        for key_str in ('GALSIM_FLUX_ERR_', 'GALSIM_MAG_ERR_', 'GALSIM_RES_'):
            self._update_dict(key_str, np.ones(len(self._obj_id)) * -1)
        for key_str in (
            'GALSIM_GAL_ELL_',
            'GALSIM_GAL_ELL_ERR_',
            'GALSIM_GAL_ELL_UNCORR_',
            'GALSIM_PSF_ELL_'
        ):
            self._update_dict(key_str, np.ones((len(self._obj_id), 2)) * -10.0)
        self._update_dict(
            'GALSIM_FLAGS_',
            np.ones(len(self._obj_id), dtype='int16'),
        )

        for idx, id_tmp in enumerate(self._obj_id):
            ind = np.where(id_tmp == galsim_id)[0]
            if len(ind) > 0:

                for key in self._key_ends:

                    gcf_data = galsim_cat_file.get_data(key)

                    if key == 'ORIGINAL_PSF':

                        uncorr_g = (
                            gcf_data['gal_uncorr_g1'][ind[0]],
                            gcf_data['gal_uncorr_g2'][ind[0]]
                        )
                        psf_sig = gcf_data['gal_sigma'][ind[0]]
                        self._add2dict(f'GALSIM_PSF_ELL_{key}', uncorr_g, idx)
                        self._add2dict(f'GALSIM_PSF_SIGMA_{key}', psf_sig, idx)

                    else:

                        g = (
                            gcf_data['gal_g1'][ind[0]],
                            gcf_data['gal_g2'][ind[0]]
                        )
                        g_err = (
                            gcf_data['gal_g1_err'][ind[0]],
                            gcf_data['gal_g2_err'][ind[0]]
                        )
                        self._add2dict(f'GALSIM_GAL_ELL_{key}', g, idx)
                        self._add2dict(f'GALSIM_GAL_ELL_ERR_{key}', g_err, idx)

                        uncorr_g = (
                            gcf_data['gal_uncorr_g1'][ind[0]],
                            gcf_data['gal_uncorr_g2'][ind[0]],
                        )
                        self._add2dict(
                            f'GALSIM_GAL_ELL_UNCORR_{key}',
                            uncorr_g,
                            idx,
                        )

                        sigma = gcf_data['gal_sigma'][ind[0]]
                        self._add2dict(f'GALSIM_GAL_SIGMA_{key}', sigma, idx)

                        psf_g = (
                            gcf_data['psf_g1'][ind[0]],
                            gcf_data['psf_g2'][ind[0]]
                        )
                        psf_sigma = gcf_data['psf_sigma'][ind[0]]
                        self._add2dict(f'GALSIM_PSF_ELL_{key}', psf_g, idx)
                        self._add2dict(
                            f'GALSIM_PSF_SIGMA_{key}',
                            psf_sigma,
                            idx,
                        )

                        flux = gcf_data['gal_flux'][ind[0]]
                        flux_err = gcf_data['gal_flux_err'][ind[0]]
                        self._add2dict(f'GALSIM_FLUX_{key}', flux, idx)
                        self._add2dict(f'GALSIM_FLUX_ERR_{key}', flux_err, idx)

                        mag = gcf_data['gal_mag'][ind[0]]
                        mag_err = gcf_data['gal_mag_err'][ind[0]]
                        self._add2dict(f'GALSIM_MAG_{key}', mag, idx)
                        self._add2dict(f'GALSIM_MAG_ERR_{key}', mag_err, idx)

                        flags = gcf_data['gal_flag'][ind[0]]
                        self._add2dict(f'GALSIM_FLAGS_{key}', flags, idx)

                        res = gcf_data['gal_resolution'][ind[0]]
                        self._add2dict(f'GALSIM_RES_{key}', res, idx)

        galsim_cat_file.close()

    def _save_psf_data(self, galaxy_psf_path):
        """Save PSF data.

        Save the PSF catalog into the final one.

        Parameters
        ----------
        galaxy_psf_path : str
            Path to the PSF catalog to save.

        """
        galaxy_psf_cat = SqliteDict(galaxy_psf_path)

        max_epoch = np.max(self.final_cat_file.get_data()['N_EPOCH']) + 1

        self._output_dict = {
            f'PSF_ELL_{idx + 1}': np.ones((len(self._obj_id), 2)) * -10.
            for idx in range(max_epoch)
        }
        self._output_dict = {
            **self._output_dict,
            **{
                f'PSF_FWHM_{idx + 1}': np.zeros(len(self._obj_id))
                for idx in range(max_epoch)
            }
        }
        self._output_dict = {
            **self._output_dict,
            **{
                f'PSF_FLAG_{idx + 1}': np.ones(
                    len(self._obj_id), dtype='int16'
                )
                for idx in range(max_epoch)
            }
        }

        for idx, id_tmp in enumerate(self._obj_id):

            if galaxy_psf_cat[str(id_tmp)] == 'empty':
                continue

            for epoch, key in enumerate(galaxy_psf_cat[str(id_tmp)].keys()):

                gpc_data = galaxy_psf_cat[str(id_tmp)][key]

                if gpc_data['SHAPES']['FLAG_PSF_HSM'] != 0:
                    continue

                e_psf = (
                    gpc_data['SHAPES']['E1_PSF_HSM'],
                    gpc_data['SHAPES']['E2_PSF_HSM']
                )
                self._add2dict(f'PSF_ELL_{epoch + 1}', e_psf, idx)

                psf_fwhm = gpc_data['SHAPES']['SIGMA_PSF_HSM'] * 2.355
                self._add2dict(f'PSF_FWHM_{epoch + 1}', psf_fwhm, idx)

                flag_psf = gpc_data['SHAPES']['FLAG_PSF_HSM']
                self._add2dict(f'PSF_FLAG_{epoch + 1}', flag_psf, idx)

        galaxy_psf_cat.close()
