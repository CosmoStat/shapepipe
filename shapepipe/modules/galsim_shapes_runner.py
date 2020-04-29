# -*- coding: utf-8 -*-

"""GALSIM SHAPES RUNNER

This file contains methods to measure shapes with Galsim.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io
from sqlitedict import SqliteDict

import re

import numpy as np

import galsim


def get_jacob(wcs, ra, dec):
    """ Get jacobian

    Return the jacobian of the wcs at the required position.

    Parameters
    ----------
    wcs : astropy.wcs.WCS
        WCS object for wich we want the jacobian.
    ra : float
        Ra position of the center of the vignet (in Deg).
    dec : float
        Dec position of the center of the vignet (in Deg).

    Returns
    -------
    galsim_jacob : galsim.wcs.BaseWCS.jacobian
        Jacobian of the WCS at the required position.

    """

    g_wcs = galsim.fitswcs.AstropyWCS(wcs=wcs)
    world_pos = galsim.CelestialCoord(ra=ra*galsim.angle.degrees,
                                      dec=dec*galsim.angle.degrees)
    galsim_jacob = g_wcs.jacobian(world_pos=world_pos)

    return galsim_jacob


def stack_psfs(psfs, weights):
    """ Stack PSFs

    Perform the weighted average stacking of the PSFs.

    Parameters
    ----------
    psfs : numpy.ndarray
        Array containing the PSF for all epochs of one object.
    weights : numpy.ndarray
        Array containing the weights for all epochs of one objects.

    Returns
    -------
    psf_sum : np.ndarray
        Stacked PSF.

    """

    n_epoch = len(psfs)

    w_sum = 0
    psf_sum = np.zeros_like(psfs[0])
    for i in range(n_epoch):
        s = np.shape(weights[i])
        cx, cy = round(s[0]/2.), round(s[1]/2.)
        w = np.mean(weights[i][cx-3:cx+3, cy-3:cy+3])
        if w < 0:
            print("Weight < 0 !")
            return 'Error'
        if w == 0:
            return 'Error'
        psf_tmp = psfs[i]/np.sum(psfs[i])
        psf_sum += w * psfs[i]
        w_sum += w

    psf_sum /= w_sum

    return psf_sum


def do_galsim_shapes(gal, gal_sig, psfs, psfs_sigma, weights, flags,
                     pixel_scale):
    """ Do ngmix metacal

    Do the metacalibration on a multi-epoch object and return the join shape
    measurement with ngmix

    Parameters
    ---------
    gal : numpy.ndarray
        Galaxy vignet from the stack.
    psfs : list
        List of the PSF vignets.
    psfs_sigma : list
        List of the sigma PSFs.
    weights : list
        List of the weight vignets.
    flags : list
        List of the flag vignets.
    jacob_list : list
        List of the jacobians.
    prior : ngmix.priors
        Priors for the fitting parameters.

    Returns
    -------
    metacal_res : dict
        Dictionary containing the results of ngmix metacal.

    """

    g_gal = galsim.Image(gal, scale=pixel_scale)

    psf = stack_psfs(psfs, weights)
    if psf == 'Error':
        return 'Error'
    psf_sig = np.mean(psfs_sigma)
    g_psf = galsim.Image(psf, scale=pixel_scale)

    weight = np.sum(weights, 0)
    flag = np.sum(flags, 0)
    weight[np.where(flag != 0)] = 0
    g_weight = galsim.Image(weight)

    res_gal = galsim.hsm.EstimateShear(g_gal, g_psf, weight=g_weight,
                                       strict=False)

    return res_gal


def compile_results(results):
    """ Compile results

    Prepare the results of ngmix before saving.

    Parameters
    ----------
    results : dict
        Dictionary containing the results of ngmix metacal.

    Returns
    -------
    output_dict : dict
        Dictionary containing ready to be saved.

    """

    output_dict = {'id': [],
                   'gal_g1': [], 'gal_g2': [],
                   'gal_uncorr_g1': [], 'gal_uncorr_g2': [],
                   'gal_sigma': [],
                   'gal_resolution': [],
                   'gal_flag': [],
                   'psf_g1': [], 'psf_g2': [],
                   'psf_sigma': []}
    for i in range(len(results)):
        output_dict['id'].append(results[i]['obj_id'])
        if (results[i]['gal'].error_message == ''):
            try:
                gal_shapes = galsim.Shear(e1=results[i]['gal'].corrected_e1,
                                          e2=results[i]['gal'].corrected_e2)
                output_dict['gal_g1'].append(gal_shapes.g1)
                output_dict['gal_g2'].append(gal_shapes.g2)
                gal_err = 0
            except Exception:
                output_dict['gal_g1'].append(results[i]['gal'].corrected_e1)
                output_dict['gal_g2'].append(results[i]['gal'].corrected_e2)
                gal_err = 2
            # output_dict['gal_g1'].append(results[i]['gal'].corrected_g1)
            # output_dict['gal_g2'].append(results[i]['gal'].corrected_g2)
            # gal_err = 0
        else:
            output_dict['gal_g1'].append(-10.)
            output_dict['gal_g2'].append(-10.)
            gal_err = 1

        output_dict['gal_uncorr_g1'].append(results[i]['gal'].observed_shape.g1)
        output_dict['gal_uncorr_g2'].append(results[i]['gal'].observed_shape.g2)
        output_dict['gal_sigma'].append(results[i]['gal'].moments_sigma)
        output_dict['gal_flag'].append(gal_err)
        output_dict['gal_resolution'].append(results[i]['gal'].resolution_factor)
        output_dict['psf_g1'].append(results[i]['gal'].psf_shape.g1)
        output_dict['psf_g2'].append(results[i]['gal'].psf_shape.g2)
        output_dict['psf_sigma'].append(results[i]['gal'].psf_sigma)

    return output_dict


def save_results(output_dict, output_name):
    """ Save results

    Save the results into a fits file.

    Parameters
    ----------
    output_dict : dict
        Dictionary containing the results.
    output_name : str
        Name of the output file.

    """

    f = io.FITSCatalog(output_name,
                       open_mode=io.BaseCatalog.OpenMode.ReadWrite)

    # for key in output_dict.keys():
    f.save_as_fits(output_dict, ext_name='RESULTS')


def process(tile_cat_path, sm_cat_path, gal_vignet_path, bkg_vignet_path,
            psf_vignet_path, weight_vignet_path, flag_vignet_path,
            w_log):
    """ Process

    Process function.

    Parameters
    ----------
    tile_cat_path : str
        Path to the tile SExtractor catalog.
    gal_vignet_path : str
        Path to the galaxy vignets catalog.
    bkg_vignet_path : str
        Path to the background vignets catalog.
    psf_vignet_path : str
        Path to the PSF vignets catalog.
    weight_vignet_path : str
        Path to the weight vignets catalog.
    flag_vignet_path : str
        Path to the flag vignets catalog.
    f_wcs_path : str
        Path to the log file containing the WCS for each CCDs.

    Returns
    -------
    final_res : dict
        Dictionary containing the ngmix metacal results.

    """

    tile_cat = io.FITSCatalog(tile_cat_path, SEx_catalog=True)
    tile_cat.open()
    obj_id = np.copy(tile_cat.get_data()['NUMBER'])
    tile_vign = np.copy(tile_cat.get_data()['VIGNET'])
    tile_flag = np.copy(tile_cat.get_data()['FLAGS'])
    tile_imaflag = np.copy(tile_cat.get_data()['IMAFLAGS_ISO'])
    tile_ra = np.copy(tile_cat.get_data()['XWIN_WORLD'])
    tile_dec = np.copy(tile_cat.get_data()['YWIN_WORLD'])
    tile_n_epoch = np.copy(tile_cat.get_data()['N_EPOCH'])
    tile_fwhm = np.copy(tile_cat.get_data()['FWHM_IMAGE'])
    tile_cat.close()
    sm_cat = io.FITSCatalog(sm_cat_path, SEx_catalog=True)
    sm_cat.open()
    sm = np.copy(sm_cat.get_data()['SPREAD_MODEL'])
    sm_err = np.copy(sm_cat.get_data()['SPREADERR_MODEL'])
    sm_cat.close()
    # f_wcs_file = np.load(f_wcs_path, allow_pickle=True).item()
    bkg_vign_cat = SqliteDict(bkg_vignet_path)
    psf_vign_cat = SqliteDict(psf_vignet_path)
    weight_vign_cat = SqliteDict(weight_vignet_path)
    flag_vign_cat = SqliteDict(flag_vignet_path)

    final_res = []
    # prior = get_prior()
    output_vignet = {'PSF': [], 'WEIGHT': [], 'FLAG': [], 'GAL': [], 'id': [],
                     'gal_flag': []}
    for i_tile, id_tmp in enumerate(obj_id):
        res = {}
        # Preselection step
        # if (tile_flag[i_tile] > 1) or (tile_imaflag[i_tile] > 0):
        #     continue
        # if (sm[i_tile] + (5. / 3.) * sm_err[i_tile] < 0.01) and
        # (np.abs(sm[i_tile] + (5. / 3.) * sm_err[i_tile]) > 0.003):
        #     continue
        # if sm[i_tile] + (5. / 3.) * sm_err[i_tile] > 0.01:
        #     gal_flag = 1
        # else:
        #     gal_flag = 0

        psf_vign = []
        sigma_psf = []
        weight_vign = []
        flag_vign = []
        if (psf_vign_cat[str(id_tmp)] == 'empty'):
            continue

        skip = False
        psf_expccd_name = list(psf_vign_cat[str(id_tmp)].keys())
        for expccd_name_tmp in psf_expccd_name:

            psf_vign.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            sigma_psf.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['SHAPES']['SIGMA_PSF_HSM'])

            weight_vign.append(weight_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])

            tile_vign_tmp = np.copy(tile_vign[i_tile])
            flag_vign_tmp = flag_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            flag_vign_tmp[np.where(tile_vign_tmp == -1e30)] = 2**10
            v_flag_tmp = flag_vign_tmp.ravel()
            if len(np.where(v_flag_tmp != 0)[0])/(51*51) > 1/3.:
                skip = True
                continue
            flag_vign.append(flag_vign_tmp)

        if len(psf_vign) != tile_n_epoch[i_tile]:
            continue
        if skip:
            skip = False
            continue

        try:
            res['gal'] = do_galsim_shapes(tile_vign[i_tile],
                                          tile_fwhm[i_tile]/2.335,
                                          psf_vign,
                                          sigma_psf,
                                          weight_vign,
                                          flag_vign,
                                          0.186)
        except Exception:
            w_log.info('Galsim fail on object {}'.format(id_tmp))
            continue

        if res['gal'] == 'Error':
            w_log.info('Something went wrong with the psf on object id : '
                       '{}.'.format(id_tmp))
            continue

        res['obj_id'] = id_tmp

        final_res.append(res)

    bkg_vign_cat.close()
    flag_vign_cat.close()
    weight_vign_cat.close()
    psf_vign_cat.close()

    return final_res


@module_runner(input_module=['sextractor_runner', 'psfexinterp_runner',
                             'vignetmaker_runner'],
               version='0.0.1',
               file_pattern=['tile_sexcat', 'image', 'exp_background',
                             'galaxy_psf', 'weight', 'flag'],
               file_ext=['.fits', '.sqlite', '.sqlite', '.sqlite', '.sqlite',
                         '.sqlite'],
               depends=['numpy', 'ngmix', 'galsim'])
def galsim_shapes_runner(input_file_list, run_dirs, file_number_string,
                         config, w_log):

    output_name = (run_dirs['output'] + '/' + 'galsim' +
                   file_number_string + '.fits')

    # f_wcs_path = config.getexpanded('NGMIX_RUNNER', 'LOG_WCS')

    metacal_res = process(*input_file_list, w_log)
    res_dict = compile_results(metacal_res)
    save_results(res_dict, output_name)

    return None, None
