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
        w.wcs.cd = np.array([[0.187/3600, 0],
                            [0, 0.187/3600]])
        w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        w.wcs.cunit = ['deg', 'deg']
        w._naxis = [10000, 10000]

        return w

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


def save_ngmix_data(final_cat_file, ngmix_cat_path):
    """Save ngmix Data

    Save the ngmix catalog into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    ngmix_cat_path : str
        Path to ngmix catalog to save.

    """

    final_cat_file.open()
    obj_id = np.copy(final_cat_file.get_data()['NUMBER'])

    ngmix_cat_file =file_io.FITSCatalog(ngmix_cat_path)
    ngmix_cat_file.open()
    ngmix_n_epoch = ngmix_cat_file.get_data()['n_epoch_model']
    ngmix_mcal_flags = ngmix_cat_file.get_data()['mcal_flags']

    ngmix_mom_fail = ngmix_cat_file.get_data()['moments_fail']

    ngmix_id = ngmix_cat_file.get_data()['id']

    keys = ['1M', '1P', '2M', '2P', 'NOSHEAR']
    output_dict = {
        f'NGMIX_ELL_{key}': np.ones((len(obj_id), 2)) * -10.0
        for key in keys
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIX_ELL_ERR_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_T_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_T_ERR_{key}': np.ones(len(obj_id)) * 1e30 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_Tpsf_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_SNR_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_FLUX_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_FLUX_ERR_{key}': np.ones(len(obj_id)) * -1 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_MAG_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_MAG_ERR_{key}': np.ones(len(obj_id)) * -1 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIX_FLAGS_{key}': np.ones(len(obj_id), dtype='int16')
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIX_ELL_PSFo_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'NGMIX_T_PSFo_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict['NGMIX_N_EPOCH'] = np.zeros(len(obj_id))
    output_dict['NGMIX_MCAL_FLAGS'] = np.zeros(len(obj_id))
    output_dict['NGMIX_MOM_FAIL'] = np.zeros(len(obj_id))

    for idx, id_tmp in enumerate(obj_id):
        ind = np.where(id_tmp == ngmix_id)[0]
        if len(ind) > 0:

            for key in keys:
                output_dict[f'NGMIX_ELL_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1'][ind[0]]
                )
                output_dict[f'NGMIX_ELL_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2'][ind[0]]
                )
                output_dict[f'NGMIX_ELL_ERR_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1_err'][ind[0]]
                )
                output_dict[f'NGMIX_ELL_ERR_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2_err'][ind[0]]
                )
                output_dict[f'NGMIX_T_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T'][ind[0]]
                )
                output_dict[f'NGMIX_T_ERR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T_err'][ind[0]]
                )
                output_dict[f'NGMIX_Tpsf_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['Tpsf'][ind[0]]
                )
                output_dict[f'NGMIX_SNR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['s2n'][ind[0]]
                )
                output_dict[f'NGMIX_FLUX_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flux'][ind[0]]
                )
                output_dict[f'NGMIX_FLUX_ERR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flux_err'][ind[0]]
                )
                output_dict[f'NGMIX_MAG_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['mag'][ind[0]]
                )
                output_dict[f'NGMIX_MAG_ERR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['mag_err'][ind[0]]
                )
                output_dict[f'NGMIX_FLAGS_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flags'][ind[0]]
                )
                output_dict[f'NGMIX_ELL_PSFo_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1_psfo_ngmix'][ind[0]]
                )
                output_dict[f'NGMIX_ELL_PSFo_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2_psfo_ngmix'][ind[0]]
                )
                output_dict[f'NGMIX_T_PSFo_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T_psfo_ngmix'][ind[0]]
                )

            output_dict[f'NGMIX_N_EPOCH'][idx] = ngmix_n_epoch[ind[0]]
            output_dict[f'NGMIX_MCAL_FLAGS'][idx] = ngmix_mcal_flags[ind[0]]
            output_dict[f'NGMIX_MOM_FAIL'][idx] = ngmix_mom_fail[ind[0]]

    for key in output_dict.keys():
        final_cat_file.add_col(key, output_dict[key])

    final_cat_file.close()
    ngmix_cat_file.close()


def save_ngmix_mom_shapes(final_cat_file, ngmix_cat_path):
    """Save ngmix MOM Shapes.

    Save the ngmix catalog into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    ngmix_cat_path : str
        Path to ngmix catalog to save.

    """

    final_cat_file.open()
    obj_id = np.copy(final_cat_file.get_data()['NUMBER'])

    ngmix_cat_file =file_io.FITSCatalog(ngmix_cat_path)
    ngmix_cat_file.open()
    ngmix_mcal_flags = ngmix_cat_file.get_data()['mcal_flags']

    ngmix_id = ngmix_cat_file.get_data()['id']

    keys = ['1M', '1P', '2M', '2P', 'NOSHEAR']
    output_dict = {
        f'NGMIXM_ELL_{key}': np.ones((len(obj_id), 2)) * -10. for key in keys
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIXM_ELL_ERR_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_T_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_T_ERR_{key}': np.ones(len(obj_id)) * 1e30 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_Tpsf_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_SNR_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_FLUX_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_FLUX_ERR_{key}': -1 * np.ones(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIXM_FLAGS_{key}': np.ones(len(obj_id), dtype='int16')
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{
            f'NGMIXM_ELL_PSFo_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'NGMIXM_T_PSFo_{key}': np.zeros(len(obj_id)) for key in keys}
    }

    output_dict['NGMIXM_MCAL_FLAGS'] = np.zeros(len(obj_id))

    for idx, id_tmp in enumerate(obj_id):
        ind = np.where(id_tmp == ngmix_id)[0]
        if len(ind) > 0:
            for key in keys:
                output_dict[f'NGMIXM_ELL_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1'][ind[0]]
                )
                output_dict[f'NGMIXM_ELL_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2'][ind[0]]
                )
                output_dict[f'NGMIXM_ELL_ERR_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1_err'][ind[0]]
                )
                output_dict[f'NGMIXM_ELL_ERR_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2_err'][ind[0]]
                )
                output_dict[f'NGMIXM_T_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T'][ind[0]]
                )
                output_dict[f'NGMIXM_T_ERR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T_err'][ind[0]]
                )
                output_dict[f'NGMIXM_Tpsf_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['Tpsf'][ind[0]]
                )
                output_dict[f'NGMIXM_SNR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['s2n'][ind[0]]
                )
                output_dict[f'NGMIXM_FLUX_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flux'][ind[0]]
                )
                output_dict[f'NGMIXM_FLUX_ERR_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flux_err'][ind[0]]
                )
                output_dict[f'NGMIXM_FLAGS_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['flags'][ind[0]]
                )
                output_dict[f'NGMIXM_ELL_PSFo_{key}'][idx][0] = (
                    ngmix_cat_file.get_data(key)['g1_psfo_ngmix'][ind[0]]
                )
                output_dict[f'NGMIXM_ELL_PSFo_{key}'][idx][1] = (
                    ngmix_cat_file.get_data(key)['g2_psfo_ngmix'][ind[0]]
                )
                output_dict[f'NGMIXM_T_PSFo_{key}'][idx] = (
                    ngmix_cat_file.get_data(key)['T_psfo_ngmix'][ind[0]]
                )

            output_dict['NGMIXM_MCAL_FLAGS'][idx] = ngmix_mcal_flags[ind[0]]

    for key in output_dict.keys():
        final_cat_file.add_col(key, output_dict[key])

    final_cat_file.close()
    ngmix_cat_file.close()


def save_galsim_shapes(final_cat_file, galsim_cat_path):
    """Save GalSim Shapes.

    Save the GalSim catalog into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    ngmix_cat_path : str
        Path to ngmix catalog to save.

    """

    final_cat_file.open()
    obj_id = np.copy(final_cat_file.get_data()['NUMBER'])

    galsim_cat_file =file_io.FITSCatalog(galsim_cat_path)
    galsim_cat_file.open()

    galsim_id = galsim_cat_file.get_data()['id']

    keys = galsim_cat_file.get_ext_name()[1:]
    output_dict = {
        f'GALSIM_GAL_ELL_{key}': np.ones((len(obj_id), 2)) * -10.
        for key in keys
    }
    output_dict = {
        **output_dict,
        **{
            f'GALSIM_GAL_ELL_ERR_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{
            f'GALSIM_GAL_ELL_UNCORR_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_GAL_SIGMA_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{
            f'GALSIM_PSF_ELL_{key}': np.ones((len(obj_id), 2)) * -10.
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_PSF_SIGMA_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_FLUX_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_FLUX_ERR_{key}': np.ones(len(obj_id)) * -1 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_MAG_{key}': np.zeros(len(obj_id)) for key in keys}
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_MAG_ERR_{key}': np.ones(len(obj_id)) * -1 for key in keys}
    }
    output_dict = {
        **output_dict,
        **{
            f'GALSIM_FLAGS_{key}': np.ones(len(obj_id), dtype='int16')
            for key in keys
        }
    }
    output_dict = {
        **output_dict,
        **{f'GALSIM_RES_{key}': np.ones(len(obj_id)) * -1. for key in keys}
    }
    for idx, id_tmp in enumerate(obj_id):
        ind = np.where(id_tmp == galsim_id)[0]
        if len(ind) > 0:
            for key in keys:
                if key == 'ORIGINAL_PSF':
                    output_dict[f'GALSIM_PSF_ELL_{key}'][idx][0] = (
                        galsim_cat_file.get_data(key)['gal_uncorr_g1'][ind[0]]
                    )
                    output_dict[f'GALSIM_PSF_ELL_{key}'][idx][1] = (
                        galsim_cat_file.get_data(key)['gal_uncorr_g2'][ind[0]]
                    )
                    output_dict[f'GALSIM_PSF_SIGMA_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_sigma'][ind[0]]
                    )
                else:
                    output_dict[f'GALSIM_GAL_ELL_{key}'][idx][0] = (
                        galsim_cat_file.get_data(key)['gal_g1'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_ELL_{key}'][idx][1] = (
                        galsim_cat_file.get_data(key)['gal_g2'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_ELL_ERR_{key}'][idx][0] = (
                        galsim_cat_file.get_data(key)['gal_g1_err'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_ELL_ERR_{key}'][idx][1] = (
                        galsim_cat_file.get_data(key)['gal_g2_err'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_ELL_UNCORR_{key}'][idx][0] = (
                        galsim_cat_file.get_data(key)['gal_uncorr_g1'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_ELL_UNCORR_{key}'][idx][1] = (
                        galsim_cat_file.get_data(key)['gal_uncorr_g2'][ind[0]]
                    )
                    output_dict[f'GALSIM_GAL_SIGMA_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_sigma'][ind[0]]
                    )
                    output_dict[f'GALSIM_PSF_ELL_{key}'][idx][0] = (
                        galsim_cat_file.get_data(key)['psf_g1'][ind[0]]
                    )
                    output_dict[f'GALSIM_PSF_ELL_{key}'][idx][1] = (
                        galsim_cat_file.get_data(key)['psf_g2'][ind[0]]
                    )
                    output_dict[f'GALSIM_PSF_SIGMA_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['psf_sigma'][ind[0]]
                    )
                    output_dict[f'GALSIM_FLUX_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_flux'][ind[0]]
                    )
                    output_dict[f'GALSIM_FLUX_ERR_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_flux_err'][ind[0]]
                    )
                    output_dict[f'GALSIM_MAG_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_mag'][ind[0]]
                    )
                    output_dict[f'GALSIM_MAG_ERR_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_mag_err'][ind[0]]
                    )
                    output_dict[f'GALSIM_FLAGS_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_flag'][ind[0]]
                    )
                    output_dict[f'GALSIM_RES_{key}'][idx] = (
                        galsim_cat_file.get_data(key)['gal_resolution'][ind[0]]
                    )

    for key in output_dict.keys():
        final_cat_file.add_col(key, output_dict[key])

    final_cat_file.close()
    galsim_cat_file.close()


def save_psf_data(final_cat_file, galaxy_psf_path):
    """ Save PSF data

    Save the PSF catalog into the final one.

    Parameters
    ----------
    final_cat_file : file_io.FITSCatalog
        Final catalog.
    galaxy_psf_path : str
        Path to the PSF catalog to save.

    """

    final_cat_file.open()
    obj_id = np.copy(final_cat_file.get_data()['NUMBER'])
    max_epoch = np.max(final_cat_file.get_data()['N_EPOCH']) + 1

    galaxy_psf_cat = SqliteDict(galaxy_psf_path)
    output_dict = {
        f'PSF_ELL_{idx + 1}': np.ones((len(obj_id), 2)) * -10.
        for idx in range(max_epoch)
    }
    output_dict = {
        **output_dict,
        **{
            f'PSF_FWHM_{idx + 1}': np.zeros(len(obj_id))
            for idx in range(max_epoch)
        }
    }
    output_dict = {
        **output_dict,
        **{
            f'PSF_FLAG_{idx + 1}': np.ones(len(obj_id), dtype='int16')
            for idx in range(max_epoch)
        }
    }
    for idx, id_tmp in enumerate(obj_id):
        if galaxy_psf_cat[str(id_tmp)] == 'empty':
            continue
        for epoch, key in enumerate(galaxy_psf_cat[str(id_tmp)].keys()):
            if galaxy_psf_cat[str(id_tmp)][key]['SHAPES']['FLAG_PSF_HSM'] != 0:
                continue
            output_dict[f'PSF_ELL_{epoch + 1}'][idx][0] = (
                galaxy_psf_cat[str(id_tmp)][key]['SHAPES']['E1_PSF_HSM']
            )
            output_dict[f'PSF_ELL_{epoch + 1}'][idx][1] = (
                galaxy_psf_cat[str(id_tmp)][key]['SHAPES']['E2_PSF_HSM']
            )
            output_dict[f'PSF_FWHM_{epoch + 1}'][idx] = (
                galaxy_psf_cat[str(id_tmp)][key]['SHAPES']['SIGMA_PSF_HSM']
                * 2.355
            )
            output_dict[f'PSF_FLAG_{epoch + 1}'][idx] = (
                galaxy_psf_cat[str(id_tmp)][key]['SHAPES']['FLAG_PSF_HSM']
            )

    for key in output_dict.keys():
        final_cat_file.add_col(key, output_dict[key])

    final_cat_file.close()
    galaxy_psf_cat.close()
