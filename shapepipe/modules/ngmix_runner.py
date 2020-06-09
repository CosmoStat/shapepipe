# -*- coding: utf-8 -*-

""" NGMIX RUNNER

This file contains methods to run ngmix for shape measurement.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io
from sqlitedict import SqliteDict

import re

import numpy as np
from numpy.random import uniform as urand

import ngmix
from ngmix.observation import Observation, ObsList, MultiBandObsList
from ngmix.fitting import LMSimple

import galsim


def get_prior():
    """ Get prior

    Return prior for the different parameters

    Return
    ------
    prior : ngmix.priors
        Priors for the different parameters.

    """

    # prior on ellipticity.  The details don't matter, as long
    # as it regularizes the fit.  This one is from Bernstein & Armstrong 2014
    g_sigma = 0.4
    g_prior = ngmix.priors.GPriorBA(g_sigma)

    # 2-d gaussian prior on the center
    # row and column center (relative to the center of the jacobian, which
    # would be zero)
    # and the sigma of the gaussians
    # units same as jacobian, probably arcsec
    row, col = 0.0, 0.0
    row_sigma, col_sigma = 0.186, 0.186  # pixel size of DES
    cen_prior = ngmix.priors.CenPrior(row, col, row_sigma, col_sigma)

    # T prior.  This one is flat, but another uninformative you might
    # try is the two-sided error function (TwoSidedErf)
    Tminval = -10.0  # arcsec squared
    Tmaxval = 1.e6
    T_prior = ngmix.priors.FlatPrior(Tminval, Tmaxval)

    # similar for flux.  Make sure the bounds make sense for
    # your images
    Fminval = -1.e4
    Fmaxval = 1.e9
    F_prior = ngmix.priors.FlatPrior(Fminval, Fmaxval)

    # now make a joint prior.  This one takes priors
    # for each parameter separately
    prior = ngmix.joint_prior.PriorSimpleSep(cen_prior,
                                             g_prior,
                                             T_prior,
                                             F_prior)

    return prior


def get_guess(img, pixel_scale=0.187,
              guess_flux_unit='img',
              guess_size_type='T', guess_size_unit='sky',
              guess_centroid=True, guess_centroid_unit='sky'):
    """Get guess

    Get the guess vector for the ngmix shape measurement
    [center_x, center_y, g1, g2, size_T, flux]
    No guess are given for the ellipticity (0., 0.)

    Parameters
    ----------
    img : numpy.ndarray
        Array containing the image
    pixel_scale : float
        Approximation of the pixel scale
    guess_flux_unit : string
        If 'img' return the flux in pixel unit
        if 'sky' return the flux in arcsec^-2
    guess_size_type : string
        if 'T' return the size in quadrupole moments definition (2 * sigma**2)
        if 'sigma' return moments sigma
    guess_size_unit : string
        If 'img' return the size in pixel unit
        if 'sky' return the size in arcsec
    guess_centroid : bool
        If True, will return a guess on the object centroid
        if False, will return the image center
    guess_centroid_unit : string
        If 'img' return the centroid in pixel unit
        if 'sky' return the centroid in arcsec

    Returns
    -------
    guess : numpy.ndarray
        Return the guess array : [center_x, center_y, g1, g2, size_T, flux]
    """

    galsim_img = galsim.Image(img, scale=pixel_scale)

    hsm_shape = galsim.hsm.FindAdaptiveMom(galsim_img, strict=False)

    error_msg = hsm_shape.error_message

    if error_msg != '':
        raise galsim.hsm.GalSimHSMError('Error in adaptive moments :\n{}'.format(error_msg))

    if guess_flux_unit == 'img':
        guess_flux = hsm_shape.moments_amp
    elif guess_flux_unit == 'sky':
        guess_flux = hsm_shape.moments_amp/pixel_scale**2.
    else:
        raise ValueError("guess_flux_unit must be in ['img', 'sky'], got : {}".format(guess_flux_unit))

    if guess_size_unit == 'img':
        size_unit = 1.
    elif guess_size_unit == 'sky':
        size_unit = pixel_scale
    else:
        raise ValueError("guess_size_unit must be in ['img', 'sky'], got : {}".format(guess_size_unit))

    if guess_size_type == 'sigma':
        guess_size = hsm_shape.moments_sigma*size_unit
    elif guess_size_type == 'T':
        guess_size = 2.*(hsm_shape.moments_sigma*size_unit)**2.

    if guess_centroid_unit == 'img':
        centroid_unit = 1.
    elif guess_centroid_unit == 'sky':
        centroid_unit = pixel_scale
    else:
        raise ValueError("guess_centroid_unit must be in ['img', 'sky'], got : {}".format(guess_centroid_unit))

    if guess_centroid:
        guess_centroid = (hsm_shape.moments_centroid-galsim_img.center) * centroid_unit
    else:
        guess_centroid = galsim_img.center * centroid_unit

    guess = np.array([guess_centroid.x,
                      guess_centroid.y,
                      0., 0.,
                      guess_size,
                      guess_flux])

    return guess


def make_galsimfit(obs, model, guess0, prior=None, lm_pars=None, ntry=5):
    """
    """

    guess = np.copy(guess0)
    fres = {}
    for it in range(ntry):
        guess[0:5] += urand(low=-0.1, high=0.1)
        guess[5:] *= (1. + urand(low=-0.1, high=0.1))
        fres['flags'] = 1
        try:
            fitter = ngmix.galsimfit.GalsimSimple(obs,
                                                  model,
                                                  prior=prior,
                                                  lm_pars=lm_pars)
            fitter.go(guess)
            fres = fitter.get_result()
        except:
            continue

        if fres['flags'] == 0:
            break

    if fres['flags'] != 0:
        raise ngmix.gexceptions.BootGalFailure("Failes to fit galaxy with galsimfit")

    fres['ntry'] = it + 1

    return fres


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


def do_ngmix_metacal(gals, psfs, psfs_sigma, weights, flags, jacob_list,
                     prior):
    """ Do ngmix metacal

    Do the metacalibration on a multi-epoch object and return the join shape
    measurement with ngmix

    Parameters
    ---------
    gals : list
        List of the galaxy vignets.
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

    pixel_scale = 0.187

    n_epoch = len(gals)

    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    # Make observation
    gal_obs_list = ObsList()
    T_guess_psf = []
    psf_res_gT = {'g_PSFo': np.array([0., 0.]),
                  'g_err_PSFo': np.array([0., 0.]),
                  'T_PSFo': 0.,
                  'T_err_PSFo': 0.}
    gal_guess = []
    gal_guess_flag = True
    wsum = 0.
    for n_e in range(n_epoch):

        psf_jacob = ngmix.Jacobian(row=(psfs[0].shape[0]-1)/2.,
                                   col=(psfs[0].shape[1]-1)/2.,
                                   wcs=jacob_list[n_e])

        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)

        # psf_T = 2. * psfs_sigma[n_e]**2.
        psf_T = psfs_sigma[n_e]*1.17741*pixel_scale

        w = np.copy(weights[n_e])
        w[np.where(flags[n_e] != 0)] = 0.

        psf_guess = np.array([0., 0., 0., 0., psf_T, 1.])
        try:
            psf_res = make_galsimfit(psf_obs, 'gauss', psf_guess, None)
        except:
            continue

        # Original PSF fit
        w_tmp = np.sum(weights[n_e])
        psf_res_gT['g_PSFo'] += psf_res['g']*w_tmp
        psf_res_gT['g_err_PSFo'] += np.array([psf_res['pars_err'][2], psf_res['pars_err'][3]])*w_tmp
        psf_res_gT['T_PSFo'] += psf_res['T']*w_tmp
        psf_res_gT['T_err_PSFo'] += psf_res['T_err']*w_tmp
        wsum += w_tmp

        # Gal guess
        try:
            gal_guess_tmp = get_guess(gals[n_e], pixel_scale=pixel_scale, guess_size_type='sigma')
        except:
            gal_guess_flag = False
            gal_guess_tmp = np.array([0., 0., 0., 0., 1, 100])

        # Recenter jacobian if necessary
        gal_jacob = ngmix.Jacobian(row=(gals[0].shape[0]-1)/2. + gal_guess_tmp[0],
                                   col=(gals[0].shape[1]-1)/2. + gal_guess_tmp[1],
                                   wcs=jacob_list[n_e])
        gal_obs = Observation(gals[n_e], weight=w, jacobian=gal_jacob,
                              psf=psf_obs)

        if gal_guess_flag:
            gal_guess_tmp[:2] = 0
            gal_guess.append(gal_guess_tmp)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)
        gal_guess_flag = True

    if wsum == 0:
        raise ZeroDivisionError('Sum of weights = 0, division by zero')

    # Normalize PSF fit output
    for key in psf_res_gT.keys():
        psf_res_gT[key] /= wsum

    # Gal guess handling
    fail_get_guess = False
    if len(gal_guess) == 0:
        fail_get_guess = True
        gal_pars = [0., 0., 0., 0., 1, 100]
    else:
        gal_pars = np.mean(gal_guess, 0)

    # boot = ngmix.bootstrap.MaxMetacalBootstrapper(gal_obs_list)

    psf_model = 'gauss'
    gal_model = 'gauss'

    # metacal specific parameters
    metacal_pars = {'types': ['noshear', '1p', '1m', '2p', '2m'],
                    'step': 0.01,
                    'psf': 'gauss',
                    'fixnoise': True,
                    'cheatnoise': False,
                    'symmetrize_psf': False}

    # maximum likelihood fitter parameters
    # parameters for the Levenberg-Marquardt fitter in scipy
    # lm_pars = {'maxfev': 2000,
    #            'xtol': 5.0e-5,
    #            'ftol': 5.0e-5}
    # max_pars = {
    #     # use scipy.leastsq for the fitting
    #     'method': 'lm',

    #     # parameters for leastsq
    #     'lm_pars': lm_pars}

    # # psf_pars = {'maxiter': 5000,
    # #             'tol': 5.0e-6}
    # psf_pars = {'maxfev': 5000,
    #            'xtol': 5.0e-6,
    #            'ftol': 5.0e-6}

    # Tguess = np.mean(T_guess_psf)*0.186**2  # size guess in arcsec

    Tguess = np.mean(T_guess_psf)

    # Tguess = 4.0*0.186**2
    ntry = 2       # retry the fit twice

    obs_dict_mcal = ngmix.metacal.get_all_metacal(gal_obs_list, **metacal_pars)
    res = {'mcal_flags': 0}

    ntry = 5

    for key in sorted(obs_dict_mcal):

        fres = make_galsimfit(obs_dict_mcal[key],
                              gal_model, gal_pars,
                              prior=prior)

        res['mcal_flags'] |= fres['flags']
        tres = {}

        for name in fres.keys():
            tres[name] = fres[name]
        tres['flags'] = fres['flags']

        wsum = 0.0
        Tpsf_sum = 0.0
        gpsf_sum = np.zeros(2)
        npsf = 0
        for obs in obs_dict_mcal[key]:

            if hasattr(obs, 'psf_nopix'):
                try:
                    psf_res = make_galsimfit(obs.psf_nopix,
                                             psf_model,
                                             np.array([0., 0., 0., 0., Tguess, 1.]),
                                             prior=None,
                                             ntry=ntry)
                except:
                    continue
                g1, g2 = psf_res['g']
                T = psf_res['T']
            else:
                try:
                    psf_res = make_galsimfit(obs.psf,
                                             psf_model,
                                             np.array([0., 0., 0., 0., Tguess, 1.]),
                                             prior=None)
                except:
                    continue
                g1, g2 = psf_res['g']
                T = psf_res['T']

            # TODO we sometimes use other weights
            twsum = obs.weight.sum()

            wsum += twsum
            gpsf_sum[0] += g1*twsum
            gpsf_sum[1] += g2*twsum
            Tpsf_sum += T*twsum
            npsf += 1

        tres['gpsf'] = gpsf_sum/wsum
        tres['Tpsf'] = Tpsf_sum/wsum

        res[key] = tres

    # result dictionary, keyed by the types in metacal_pars above
    metacal_res = res

    metacal_res.update(psf_res_gT)
    metacal_res['moments_fail'] = fail_get_guess

    return metacal_res


# def compile_results(results):
#     """ Compile results

#     Prepare the results of ngmix before saving.

#     Parameters
#     ----------
#     results : dict
#         Dictionary containing the results of ngmix metacal.

#     Returns
#     -------
#     output_dict : dict
#         Dictionary containing ready to be saved.

#     """

#     names = ['1m', '1p', '2m', '2p', 'noshear']
#     names2 = ['id', 'n_epoch_model', 'g1', 'g1_err', 'g2', 'g2_err', 'T',
#               'T_err', 'Tpsf', 's2n', 'flags', 'mcal_flags']
#     output_dict = {k: {kk: [] for kk in names2} for k in names}
#     for i in range(len(results)):
#         for name in names:
#             output_dict[name]['id'].append(results[i]['obj_id'])
#             output_dict[name]['n_epoch_model'].append(results[i]['n_epoch_model'])
#             output_dict[name]['g1'].append(results[i][name]['g'][0])
#             output_dict[name]['g1_err'].append(results[i][name]['pars_err'][2])
#             output_dict[name]['g2'].append(results[i][name]['g'][1])
#             output_dict[name]['g2_err'].append(results[i][name]['pars_err'][3])
#             output_dict[name]['T'].append(results[i][name]['T'])
#             output_dict[name]['T_err'].append(results[i][name]['T_err'])
#             output_dict[name]['Tpsf'].append(results[i][name]['Tpsf'])
#             output_dict[name]['s2n'].append(results[i][name]['s2n'])
#             output_dict[name]['flags'].append(results[i][name]['flags'])
#             output_dict[name]['mcal_flags'].append(results[i]['mcal_flags'])

#     return output_dict


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

    names = ['1m', '1p', '2m', '2p', 'noshear']
    names2 = ['id', 'n_epoch_model', 'moments_fail', 'ntry_fit',
              'g1_psfo_ngmix', 'g2_psfo_ngmix', 'T_psfo_ngmix',
              'g1_err_psfo_ngmix', 'g2_err_psfo_ngmix', 'T_err_psfo_ngmix',
              'g1', 'g1_err', 'g2', 'g2_err',
              'T', 'T_err', 'Tpsf', 'g1_psf', 'g2_psf',
              'flux', 'flux_err', 's2n',
              'flags', 'mcal_flags']
    output_dict = {k: {kk: [] for kk in names2} for k in names}
    for i in range(len(results)):
        for name in names:
            output_dict[name]['id'].append(results[i]['obj_id'])
            output_dict[name]['n_epoch_model'].append(results[i]['n_epoch_model'])
            output_dict[name]['moments_fail'].append(results[i]['moments_fail'])
            output_dict[name]['ntry_fit'].append(results[i][name]['ntry'])
            output_dict[name]['g1_psfo_ngmix'].append(results[i]['g_PSFo'][0])
            output_dict[name]['g2_psfo_ngmix'].append(results[i]['g_PSFo'][1])
            output_dict[name]['g1_err_psfo_ngmix'].append(results[i]['g_err_PSFo'][0])
            output_dict[name]['g2_err_psfo_ngmix'].append(results[i]['g_err_PSFo'][1])
            output_dict[name]['T_psfo_ngmix'].append(results[i]['T_PSFo'])
            output_dict[name]['T_err_psfo_ngmix'].append(results[i]['T_err_PSFo'])
            output_dict[name]['g1'].append(results[i][name]['g'][0])
            output_dict[name]['g1_err'].append(results[i][name]['pars_err'][2])
            output_dict[name]['g2'].append(results[i][name]['g'][1])
            output_dict[name]['g2_err'].append(results[i][name]['pars_err'][3])
            output_dict[name]['T'].append(results[i][name]['T'])
            output_dict[name]['T_err'].append(results[i][name]['T_err'])
            output_dict[name]['Tpsf'].append(results[i][name]['Tpsf'])
            output_dict[name]['g1_psf'].append(results[i][name]['gpsf'][0])
            output_dict[name]['g2_psf'].append(results[i][name]['gpsf'][1])
            output_dict[name]['flux'].append(results[i][name]['flux'])
            output_dict[name]['flux_err'].append(results[i][name]['flux_err'])

            try:
                output_dict[name]['s2n'].append(results[i][name]['s2n'])
            except:
                output_dict[name]['s2n'].append(results[i][name]['s2n_r'])
            output_dict[name]['flags'].append(results[i][name]['flags'])
            output_dict[name]['mcal_flags'].append(results[i]['mcal_flags'])

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

    for key in output_dict.keys():
        f.save_as_fits(output_dict[key], ext_name=key.upper())


def process(tile_cat_path, gal_vignet_path, bkg_vignet_path,
            psf_vignet_path, weight_vignet_path, flag_vignet_path,
            f_wcs_path, w_log, id_obj_min=-1, id_obj_max=-1):
    """ Process

    Process function.

    Parameters
    ----------
    tile_cat_path: str
        Path to the tile SExtractor catalog.
    gal_vignet_path: str
        Path to the galaxy vignets catalog.
    bkg_vignet_path: str
        Path to the background vignets catalog.
    psf_vignet_path: str
        Path to the PSF vignets catalog.
    weight_vignet_path: str
        Path to the weight vignets catalog.
    flag_vignet_path: str
        Path to the flag vignets catalog.
    f_wcs_path: str
        Path to the log file containing the WCS for each CCDs.
    w_log: log file object
        log file
    id_obj_min: int, optional, default=-1
        minimum object ID to be processed if > 0
    id_obj_max: int, optional, default=-1
        maximum object ID to be processed if > 0

    Returns
    -------
    final_res: dict
        Dictionary containing the ngmix metacal results.

    """

    tile_cat = io.FITSCatalog(tile_cat_path, SEx_catalog=True)
    tile_cat.open()
    obj_id = np.copy(tile_cat.get_data()['NUMBER'])
    tile_vign = np.copy(tile_cat.get_data()['VIGNET'])
    # tile_flag = np.copy(tile_cat.get_data()['FLAGS'])
    # tile_imaflag = np.copy(tile_cat.get_data()['IMAFLAGS_ISO'])
    tile_ra = np.copy(tile_cat.get_data()['XWIN_WORLD'])
    tile_dec = np.copy(tile_cat.get_data()['YWIN_WORLD'])
    tile_cat.close()
    # sm_cat = io.FITSCatalog(sm_cat_path, SEx_catalog=True)
    # sm_cat.open()
    # sm = np.copy(sm_cat.get_data()['SPREAD_MODEL'])
    # sm_err = np.copy(sm_cat.get_data()['SPREADERR_MODEL'])
    # sm_cat.close()
    f_wcs_file = SqliteDict(f_wcs_path)
    gal_vign_cat = SqliteDict(gal_vignet_path)
    bkg_vign_cat = SqliteDict(bkg_vignet_path)
    psf_vign_cat = SqliteDict(psf_vignet_path)
    weight_vign_cat = SqliteDict(weight_vignet_path)
    flag_vign_cat = SqliteDict(flag_vignet_path)

    final_res = []
    prior = get_prior()

    count = 0
    id_first = -1
    id_last = -1

    for i_tile, id_tmp in enumerate(obj_id):

        if id_obj_min > 0 and id_tmp < id_obj_min:
            continue
        if id_obj_max > 0 and id_tmp > id_obj_max:
            continue

        if id_first == -1:
            id_first = id_tmp
        id_last = id_tmp

        count = count + 1

        gal_vign = []
        psf_vign = []
        sigma_psf = []
        weight_vign = []
        flag_vign = []
        jacob_list = []
        if (psf_vign_cat[str(id_tmp)] == 'empty') or (gal_vign_cat[str(id_tmp)] == 'empty'):
            continue
        psf_expccd_name = list(psf_vign_cat[str(id_tmp)].keys())
        for expccd_name_tmp in psf_expccd_name:
            gal_vign_tmp = gal_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            if len(np.where(gal_vign_tmp.ravel() == 0)[0]) != 0:
                continue

            bkg_vign_tmp = bkg_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            gal_vign_sub_bkg = gal_vign_tmp - bkg_vign_tmp

            tile_vign_tmp = np.copy(tile_vign[i_tile])
            flag_vign_tmp = flag_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            flag_vign_tmp[np.where(tile_vign_tmp == -1e30)] = 2**10
            v_flag_tmp = flag_vign_tmp.ravel()
            if len(np.where(v_flag_tmp != 0)[0])/(51*51) > 1/3.:
                continue

            exp_name, ccd_n = re.split('-', expccd_name_tmp)
            jacob_tmp = get_jacob(f_wcs_file[exp_name][int(ccd_n)],
                                  tile_ra[i_tile],
                                  tile_dec[i_tile])

            gal_vign.append(gal_vign_sub_bkg)
            psf_vign.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            sigma_psf.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['SHAPES']['SIGMA_PSF_HSM'])
            weight_vign.append(weight_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            flag_vign.append(flag_vign_tmp)
            jacob_list.append(jacob_tmp)

        if len(gal_vign) == 0:
            continue
        try:
            res = do_ngmix_metacal(gal_vign,
                                   psf_vign,
                                   sigma_psf,
                                   weight_vign,
                                   flag_vign,
                                   jacob_list,
                                   prior)
        except Exception as ee:
            w_log.info('ngmix failed for object ID={}.\nMessage: {}'.format(id_tmp, ee))
            continue

        res['obj_id'] = id_tmp
        res['n_epoch_model'] = len(gal_vign)
        final_res.append(res)

    w_log.info('ngmix loop over objects finished, processed {} '
               'objects, id first/last={}/{}'.format(count,
                                                     id_first,
                                                     id_last))

    f_wcs_file.close()
    gal_vign_cat.close()
    bkg_vign_cat.close()
    flag_vign_cat.close()
    weight_vign_cat.close()
    psf_vign_cat.close()

    return final_res


@module_runner(input_module=['sextractor_runner', 'psfexinterp_runner_me',
                             'vignetmaker_runner2'],
               version='0.0.1',
               file_pattern=['tile_sexcat', 'image', 'exp_background',
                             'galaxy_psf', 'weight', 'flag'],
               file_ext=['.fits', '.sqlite', '.sqlite', '.sqlite', '.sqlite',
                         '.sqlite'],
               depends=['numpy', 'ngmix', 'galsim', 'sqlitedict'])
def ngmix_runner(input_file_list, run_dirs, file_number_string,
                 config, w_log):

    output_name = (run_dirs['output'] + '/' + 'ngmix' +
                   file_number_string + '.fits')

    f_wcs_path = config.getexpanded('NGMIX_RUNNER', 'LOG_WCS')

    id_obj_min = config.getint('NGMIX_RUNNER', 'ID_OBJ_MIN')
    id_obj_max = config.getint('NGMIX_RUNNER', 'ID_OBJ_MAX')

    metacal_res = process(*input_file_list, f_wcs_path, w_log,
                          id_obj_min=id_obj_min, id_obj_max=id_obj_max)
    res_dict = compile_results(metacal_res)
    save_results(res_dict, output_name)

    return None, None
