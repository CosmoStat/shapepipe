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
from astropy.io import fits
from astropy.wcs import WCS
import reproject
import ngmix


def mad(data, axis=None):
    """
    """

    return np.median(np.abs(data - np.median(data, axis)), axis)*1.4826


def get_gauss_2D(sigma, center=(0, 0), shape=(51, 51)):
    """
    """

    x, y = np.meshgrid(np.linspace(0, shape[0]-1, shape[0]), np.linspace(0, shape[1]-1, shape[1]))
    return np.exp(-(((x-center[0])**2. + (y-center[1])**2.))/(2. * sigma**2.)) / (sigma**2. * 2. * np.pi)


def get_wcs_from_sexcat(header_list):
    """ Get wcs from SExtractor catalog

    Read the image header from SExtractor catalog and create a wcs object.

    Parameters
    ----------
    header_list : list
        List containing all the header parameters.

    Returns
    -------
    new_wcs : astropy.wcs.WCS
        WCS object.
    """

    keys = ['CTYPE1', 'CUNIT1', 'CRVAL1', 'CRPIX1', 'CD1_1', 'CD1_2',
            'CTYPE2', 'CUNIT2', 'CRVAL2', 'CRPIX2', 'CD2_1', 'CD2_2']
    wcs_dict = {}
    for line in header_list:
        tmp = re.split('=|/', line.replace(' ', ''))
        if tmp[0] not in keys:
            continue
        wcs_dict[tmp[0]] = tmp[1].replace("'", "")

    new_wcs = WCS(naxis=2)
    new_wcs.wcs.ctype = [wcs_dict['CTYPE1'], wcs_dict['CTYPE2']]
    new_wcs.wcs.cunit = [wcs_dict['CUNIT1'], wcs_dict['CUNIT2']]
    new_wcs.wcs.crpix = [wcs_dict['CRPIX1'], wcs_dict['CRPIX2']]
    new_wcs.wcs.crval = [wcs_dict['CRVAL1'], wcs_dict['CRVAL2']]
    new_wcs.wcs.cd = [[wcs_dict['CD1_1'], wcs_dict['CD1_2']],
                      [wcs_dict['CD2_1'], wcs_dict['CD2_2']]]

    return new_wcs


def get_local_wcs(wcs, ra, dec, vign_shape):
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
    galsim_jacob = g_wcs.local(world_pos=world_pos)

    loc_wcs = WCS(naxis=2)
    loc_wcs.wcs.pc = galsim_jacob.getMatrix()/3600.
    loc_wcs.wcs.crpix = np.array(vign_shape)/2. + 0.5
    loc_wcs.wcs.crval = np.array([ra, dec])
    loc_wcs.wcs.ctype = wcs.wcs.ctype

    return loc_wcs


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


def stack_psfs(tile_loc_wcs, psfs, psfs_sigma, weights, loc_wcs):
    """ Stack PSFs

    Perform the weighted average stacking of the PSFs.

    Parameters
    ----------
    tile_loc_wcs : astropy.wcs.WCS
        Local WCS from the tile.
    psfs : numpy.ndarray
        Array containing the PSF for all epochs of one object.
    psfs_sigma : list
        List of the sigma PSFs.
    weights : numpy.ndarray
        Array containing the weights for all epochs of one objects.
    loc_wcs : list
        List of local WCS.

    Returns
    -------
    psf_sum : np.ndarray
        Stacked PSF.

    """

    n_epoch = len(psfs)

    psf_list_stack = []

    for psf, wcs in zip(psfs, loc_wcs):
        res = reproject.reproject_interp((psf, wcs), tile_loc_wcs, shape_out=psfs[0].shape)
        new_psf = res[0]
        new_psf[np.isnan(new_psf)] = 0
        psf_list_stack.append(new_psf)

    w_sum = 0
    psf_sum = np.zeros_like(psfs[0])
    psf_sigma = 0
    for i in range(n_epoch):
        s = np.shape(weights[i])
        cx, cy = int(s[0]/2.), int(s[1]/2.)
        gauss_img = get_gauss_2D(psfs_sigma[i], center=(cx, cy))
        # w = np.average(weights[i], weights=gauss_img)
        # w2 = np.sum(weights[i]*gauss_img)/
        w = np.sum(weights[i]*gauss_img)/np.sum(gauss_img[np.where(weights[i] != 0)])
        # w_old = np.mean(weights[i][cx-3:cx+3, cy-3:cy+3])
        # print("old : {} | new : {} | new2 : {}".format(w_old, w, w2))
        if w <= 0:
            raise ValueError('Error weight <= 0')
        psf_tmp = psf_list_stack[i]/np.sum(psf_list_stack[i])
        psf_sum += w * psf_tmp
        w_sum += w
        psf_sigma += psfs_sigma[i] * w

    psf_sum /= w_sum
    psf_sigma /= w_sum

    return psf_sum, psf_sigma


def check_galsim_shapes(galsim_shape, obj_id, w_log):
    """
    """

    if (galsim_shape.error_message == ''):
        if (galsim_shape.corrected_g1 != -10.) & (galsim_shape.corrected_g2 != -10):
            g1 = galsim_shape.corrected_g1
            g2 = galsim_shape.corrected_g2
            gal_err = 0
        else:
            try:
                gal_shapes = galsim.Shear(e1=galsim_shape.corrected_e1, e2=galsim_shape.corrected_e2)
                g1 = gal_shapes.g1
                g2 = gal_shapes.g2
                gal_err = 0
            except:
                g1 = galsim_shape.corrected_e1
                g2 = galsim_shape.corrected_e2
                gal_err = 2
    else:
        w_log.info('Object : {}    Error : {}'.format(obj_id, galsim_shape.error_message))
        g1 = -10.
        g2 = -10.
        gal_err = 1

    return g1, g2, gal_err


#####################
# Metacal functions #
#####################

def psf_fitter(psf_vign):
    """Psf fitter

    Function used to create a gaussian fit of the PSF.

    Parameters
    ----------
    psf_vign : numpy.array
        Array containg one vignet of psf

    """

    psf_obs = ngmix.Observation(psf_vign)
    pfitter = ngmix.fitting.LMSimple(psf_obs, 'gauss')

    shape = psf_vign.shape
    psf_pars = np.array([0., 0., 0., 0., 4., 1.])
    pfitter.go(psf_pars)

    psf_gmix_fit = pfitter.get_gmix()
    psf_obs.set_gmix(psf_gmix_fit)

    return psf_obs


def make_metacal(gal_vign, psf_vign, weight_vign, tile_jacob, option_dict):
    """Make the metacalibration

    This function call different ngmix functions to create images needed for the metacalibration.

    Parameters
    ----------
    gal_vign : numpy.array
        Array containing one vignet of galaxy
    psf_vign : numpy.array
        Array containg one vignet of psf
    option_dict : dict
        Dictionnary containg option for ngmix (keys : ['TYPES', 'FIXNOISE', 'CHEATNOISE', 'SYMMETRIZE_PSF', 'STEP'])

    """

    jacob = ngmix.Jacobian(row=(gal_vign.shape[0]-1)/2.,
                           col=(gal_vign.shape[1]-1)/2.,
                           wcs=tile_jacob)

    # psf_obs = psf_fitter(psf_vign)
    psf_obs = ngmix.Observation(psf_vign, jacobian=jacob)

    obs = ngmix.Observation(gal_vign, psf=psf_obs, weight=weight_vign, jacobian=jacob)

    obs_out = ngmix.metacal.get_all_metacal(obs,
                                            types=option_dict['TYPES'],
                                            fixnoise=option_dict['FIXNOISE'],
                                            cheatnoise=option_dict['CHEATNOISE'],
                                            step=option_dict['STEP'],
                                            psf=option_dict['PSF_KIND'])

    return obs_out


def do_galsim_shapes(gal, gal_weight, gal_sig, psfs, tile_loc_wcs, tile_jacob, loc_wcs, psfs_sigma, weights, flags, pixel_scale, do_metacal):
    """ Do ngmix metacal

    Do the metacalibration on a multi-epoch object and return the join shape
    measurement with ngmix

    Parameters
    ---------
    gal : numpy.ndarray
        Galaxy vignet from the stack.
    psfs : list
        List of the PSF vignets.
    tile_loc_wcs : astropy.wcs.WCS
        Local WCS from the tile.
    loc_wcs : list
        List of local WCS.
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

    psf, psf_sig = stack_psfs(tile_loc_wcs, psfs, psfs_sigma, weights, loc_wcs)
    if psf == 'Error':
        return 'Error'
    g_psf = galsim.Image(psf, scale=pixel_scale)

    weight = np.copy(gal_weight)

    flag = np.sum(flags, 0)
    flag[flag != 0] = 1
    flag[gal == -1e30] = 1
    g_flag = galsim.ImageI(flag, scale=pixel_scale)

    weight[np.where(flag != 0)] = 0
    weight[gal == -1e30] = 0

    inv_flag = np.ones_like(flag)
    inv_flag[flag != 0] = 0
    g_weight = galsim.Image(weight, scale=pixel_scale)

    gal[gal == -1e30] = 0

    s = np.shape(weight)
    cx, cy = int(s[0]/2.), int(s[1]/2.)
    sky_var = 1./np.average(weight, weights=get_gauss_2D(psf_sig, center=(cx, cy)))

    res_gal = {}

    if do_metacal:
        option_dict = {'TYPES': ['noshear', '1p', '1m', '2p', '2m'],
                       'FIXNOISE': True,
                       'CHEATNOISE': False,
                       'STEP': 0.01,
                       'PSF_KIND': 'gauss'}
        res = make_metacal(gal, psf, weight, tile_jacob, option_dict)

        for key in option_dict['TYPES']:
            gal_tmp = res[key].image
            g_gal = galsim.Image(gal_tmp, scale=pixel_scale)
            g_psf_mc = galsim.Image(res[key].get_psf().image, scale=pixel_scale)
            sky_var = mad(gal_tmp)**2.
            res_gal[key] = galsim.hsm.EstimateShear(g_gal,
                                                    g_psf_mc,
                                                    shear_est='KSB',
                                                    sky_var=sky_var,
                                                    weight=g_weight,
                                                    badpix=g_flag,
                                                    strict=False)
        res_gal['original_psf'] = galsim.hsm.FindAdaptiveMom(g_psf,
                                                             strict=False)
    else:
        g_gal = galsim.Image(gal, scale=pixel_scale)
        res_gal['classic'] = galsim.hsm.EstimateShear(g_gal,
                                                      g_psf,
                                                      sky_var=sky_var,
                                                      weight=g_weight,
                                                      badpix=g_flag,
                                                      strict=False)

    return res_gal


def compile_results(results, do_metacal, w_log):
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

    cat_keys = ['id',
                'gal_g1', 'gal_g2', 'gal_g1_err', 'gal_g2_err',
                'gal_uncorr_g1', 'gal_uncorr_g2',
                'gal_sigma',
                'gal_resolution',
                'gal_flag',
                'psf_g1', 'psf_g2', 'psf_sigma', 'psf_vignet', 'gal_vignet']

    if do_metacal:
        types = ['noshear', '1p', '1m', '2p', '2m']
        types += ['original_psf']
    else:
        types = ['classic']

    output_dict = {k: {kk: [] for kk in cat_keys} for k in types}

    for i in range(len(results)):
        for key in types:
            output_dict[key]['id'].append(results[i]['obj_id'])
            shapes_check = check_galsim_shapes(results[i][key],
                                               results[i]['obj_id'],
                                               w_log)
            output_dict[key]['gal_g1'].append(shapes_check[0])
            output_dict[key]['gal_g2'].append(shapes_check[1])
            output_dict[key]['gal_g1_err'].append(results[i][key].corrected_shape_err)
            output_dict[key]['gal_g2_err'].append(results[i][key].corrected_shape_err)
            output_dict[key]['gal_uncorr_g1'].append(results[i][key].observed_shape.g1)
            output_dict[key]['gal_uncorr_g2'].append(results[i][key].observed_shape.g2)
            output_dict[key]['gal_sigma'].append(results[i][key].moments_sigma)
            output_dict[key]['gal_flag'].append(shapes_check[2])
            output_dict[key]['gal_resolution'].append(results[i][key].resolution_factor)
            output_dict[key]['psf_g1'].append(results[i][key].psf_shape.g1)
            output_dict[key]['psf_g2'].append(results[i][key].psf_shape.g2)
            output_dict[key]['psf_sigma'].append(results[i][key].psf_sigma)

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

    f = io.FITSCatalog(output_name, open_mode=io.BaseCatalog.OpenMode.ReadWrite)

    for key in output_dict.keys():
        f.save_as_fits(output_dict[key], ext_name=key.upper())


def process(tile_cat_path, tile_weight_path, gal_vignet_path, bkg_vignet_path,
            psf_vignet_path, weight_vignet_path, flag_vignet_path,
            f_wcs_path, do_metacal, w_log, id_obj_min=-1, id_obj_max=-1):
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
    w_log: log file object
        log file
    id_obj_min: int, optional, default=-1
        minimum object ID to be processed if > 0
    id_obj_max: int, optional, default=-1
        maximum object ID to be processed if > 0

    Returns
    -------
    final_res : dict
        Dictionary containing the ngmix metacal results.

    """

    obj_id = fits.getdata(tile_cat_path, 2, memmap=True)['NUMBER']
    tile_vign = fits.getdata(tile_cat_path, 2, memmap=True)['VIGNET']
    tile_ra = fits.getdata(tile_cat_path, 2, memmap=True)['XWIN_WORLD']
    tile_dec = fits.getdata(tile_cat_path, 2, memmap=True)['YWIN_WORLD']
    tile_n_epoch = fits.getdata(tile_cat_path, 2, memmap=True)['N_EPOCH']
    tile_fwhm = fits.getdata(tile_cat_path, 2, memmap=True)['FWHM_IMAGE']
    tile_wcs = get_wcs_from_sexcat(fits.getdata(tile_cat_path, 1, memmap=True)[0][0])

    tile_weight = fits.getdata(tile_weight_path, 2, memmap=True)['VIGNET']
    bkg_vign_cat = SqliteDict(bkg_vignet_path)
    psf_vign_cat = SqliteDict(psf_vignet_path)
    weight_vign_cat = SqliteDict(weight_vignet_path)
    flag_vign_cat = SqliteDict(flag_vignet_path)
    f_wcs_file = SqliteDict(f_wcs_path)

    count = 0
    id_first = -1
    id_last = -1

    final_res = []
    for i_tile, id_tmp in enumerate(obj_id):

        if id_obj_min > 0 and id_tmp < id_obj_min:
            continue
        if id_obj_max > 0 and id_tmp > id_obj_max:
            continue

        if id_first == -1:
            id_first = id_tmp
        id_last = id_tmp

        count = count + 1

        res = {}
        w_log.info('{}'.format(i_tile))
        psf_vign = []
        sigma_psf = []
        weight_vign = []
        flag_vign = []
        loc_wcs_list = []
        if (psf_vign_cat[str(id_tmp)] == 'empty'):
            continue

        skip = False
        psf_expccd_name = list(psf_vign_cat[str(id_tmp)].keys())
        for expccd_name_tmp in psf_expccd_name:
            tile_vign_tmp = np.copy(tile_vign[i_tile])
            flag_vign_tmp = flag_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            flag_vign_tmp[np.where(tile_vign_tmp == -1e30)] = 2**10
            v_flag_tmp = flag_vign_tmp.ravel()
            if len(np.where(v_flag_tmp != 0)[0])/(51*51) > 1/3.:
                skip = True
                continue

            psf_vign.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            sigma_psf.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['SHAPES']['SIGMA_PSF_HSM'])
            weight_vign.append(weight_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            flag_vign.append(flag_vign_tmp)

            tile_loc_wcs = get_local_wcs(tile_wcs,
                                         tile_ra[i_tile],
                                         tile_dec[i_tile],
                                         tile_vign_tmp.shape)

            exp_name, ccd_n = re.split('-', expccd_name_tmp)
            loc_wcs_list.append(get_local_wcs(f_wcs_file[exp_name][int(ccd_n)],
                                              tile_ra[i_tile],
                                              tile_dec[i_tile],
                                              tile_vign_tmp.shape))

        if len(psf_vign) != tile_n_epoch[i_tile]:
            continue
        if skip:
            skip = False
            continue

        tile_jacob = get_jacob(tile_loc_wcs, tile_ra[i_tile], tile_dec[i_tile])

        try:
            res = do_galsim_shapes(tile_vign[i_tile],
                                   tile_weight[i_tile],
                                   tile_fwhm[i_tile]/2.335,
                                   psf_vign,
                                   tile_loc_wcs,
                                   tile_jacob,
                                   loc_wcs_list,
                                   sigma_psf,
                                   weight_vign,
                                   flag_vign,
                                   0.187,
                                   do_metacal)
        except Exception as ee:
            w_log.info('ngmix fail on object {}\nMessage : {}'.format(id_tmp, ee))
            continue

        if res == 'Error':
            w_log.info('Something went wrong with the psf on object id : {}.'.format(id_tmp))
            continue

        res['obj_id'] = id_tmp

        final_res.append(res)

    w_log.info('galsim loop over objects finished, processed {} '
               'objects, id first/last={}/{}'.format(count,
                                                     id_first,
                                                     id_last))

    bkg_vign_cat.close()
    flag_vign_cat.close()
    weight_vign_cat.close()
    psf_vign_cat.close()

    return final_res


@module_runner(input_module=['sextractor_runner', 'psfexinterp_runner', 'vignetmaker_runner'],
               version='0.0.1',
               file_pattern=['tile_sexcat', 'weight', 'image', 'exp_background', 'galaxy_psf', 'weight', 'flag'],
               file_ext=['.fits', '.fits', '.sqlite', '.sqlite', '.sqlite', '.sqlite', '.sqlite'],
               depends=['numpy', 'ngmix', 'galsim'])
def galsim_shapes_v2_runner(input_file_list, run_dirs, file_number_string,
                            config, w_log):

    output_name = run_dirs['output'] + '/' + 'galsim' + file_number_string + '.fits'

    f_wcs_path = config.getexpanded('GALSIM_SHAPES_V2_RUNNER', 'LOG_WCS')

    id_obj_min = config.getint('GALSIM_SHAPES_V2_RUNNER', 'ID_OBJ_MIN')
    id_obj_max = config.getint('GALSIM_SHAPES_V2_RUNNER', 'ID_OBJ_MAX')

    do_metacal = True

    metacal_res = process(*input_file_list, f_wcs_path, do_metacal, w_log,
                          id_obj_min=id_obj_min, id_obj_max=id_obj_max)
    res_dict = compile_results(metacal_res, do_metacal, w_log)
    save_results(res_dict, output_name)

    return None, None
