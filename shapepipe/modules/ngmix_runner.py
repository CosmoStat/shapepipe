# -*- coding: utf-8 -*-

""" NGMIX RUNNER

This file contains methods to run ngmix for shape measurement.

:Author: Axel Guinot

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io

import numpy as np

import ngmix
from ngmix.observation import Observation, ObsList, MultiBandObsList
from ngmix.fitting import LMSimple


def get_prior():
    """
    """

    # prior on ellipticity.  The details don't matter, as long
    # as it regularizes the fit.  This one is from Bernstein & Armstrong 2014
    g_sigma = 0.4
    g_prior = ngmix.priors.GPriorBA(g_sigma)


    # 2-d gaussian prior on the center
    # row and column center (relative to the center of the jacobian, which would be zero)
    # and the sigma of the gaussians
    # units same as jacobian, probably arcsec
    row, col = 0.0, 0.0
    row_sigma, col_sigma = 0.186, 0.186 # pixel size of DES
    cen_prior = ngmix.priors.CenPrior(row, col, row_sigma, col_sigma)

    # T prior.  This one is flat, but another uninformative you might
    # try is the two-sided error function (TwoSidedErf)
    Tminval = -10.0 # arcsec squared
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


def do_ngmix_metacal(gals, psfs, psfs_sigma, weights, flags, prior):
    """
    """

    n_epoch = len(gals)

    psf_jacob = ngmix.UnitJacobian(row=psfs[0].shape[0]/2., col=psfs[0].shape[1]/2.)
    gal_jacob = ngmix.UnitJacobian(row=gals[0].shape[0]/2., col=gals[0].shape[1]/2.)

    # Make observation
    gal_obs_list = ObsList()
    T_guess_psf = []
    for n_e in range(n_epoch):

        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)

        psf_T = 2. * psfs_sigma[n_e]**2.

        w = np.copy(weights[n_e])
        w[np.where(flags[n_e] != 0)] = 0.

        gal_obs = Observation(gals[n_e], weight=w, jacobian=gal_jacob, psf=psf_obs)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)

    boot = ngmix.bootstrap.MaxMetacalBootstrapper(gal_obs_list)
    psf_model = 'em3'
    gal_model = 'gauss'

    # metacal specific parameters
    metacal_pars={'types': ['noshear', '1p', '1m', '2p', '2m'],
                  'psf': 'gauss',
                  'fixnoise': True,
                  'cheatnoise': False,
                  'symmetrize_psf': False}

    # maximum likelihood fitter parameters
    # parameters for the Levenberg-Marquardt fitter in scipy
    lm_pars = {'maxfev': 2000,
               'xtol': 5.0e-5,
               'ftol': 5.0e-5}
    max_pars = {
        # use scipy.leastsq for the fitting
        'method': 'lm',

        # parameters for leastsq
        'lm_pars': lm_pars,}

    psf_pars = {'maxiter': 5000,
                'tol': 5.0e-6}

    Tguess = np.mean(T_guess_psf)#*pixel_scale**2  # size guess in arcsec
    ntry = 2       # retry the fit twice
    boot.fit_metacal(psf_model,
                     gal_model,
                     max_pars,
                     Tguess,
                     prior=prior,
                     ntry=ntry,
                     metacal_pars=metacal_pars,
                     psf_fit_pars=psf_pars,
                     psf_ntry=20)

    # result dictionary, keyed by the types in metacal_pars above
    metacal_res = boot.get_metacal_result()

    # out_dict = {}
    # out_dict['R11'] = (metacal_res['1p']['g'][0] - metacal_res['1m']['g'][0])/0.02
    # out_dict['R22']= (metacal_res['2p']['g'][1] - metacal_res['2m']['g'][1])/0.02
    # out_dict['R12'] = (metacal_res['2p']['g'][0] - metacal_res['2m']['g'][0])/0.02
    # out_dict['R21']= (metacal_res['1p']['g'][1] - metacal_res['1m']['g'][1])/0.02
    # out_dict['g_noshear'] = metacal_res['noshear']['g']
    # out_dict['g_psf'] = g_psf
    # out_dict['gal_snr'] = gal_snr

    return metacal_res


def compile_results(results):
    """
    """

    names = ['1m', '1p', '2m', '2p', 'noshear']
    names2 = ['id', 'g1', 'g1_err', 'g2', 'g2_err', 'T', 'Tpsf', 's2n', 'flags', 'mcal_flags']
    output_dict = {k: {kk: [] for kk in names2} for k in names}
    for i in range(len(results)):
        for name in names:
            output_dict[name]['id'].append(results[i]['obj_id'])
            output_dict[name]['n_epoch_model'].append(results[i]['n_epoch_model'])
            output_dict[name]['g1'].append(results[i][name]['g'][0])
            output_dict[name]['g1_err'].append(results[i][name]['pars_err'][2])
            output_dict[name]['g2'].append(results[i][name]['g'][1])
            output_dict[name]['g2_err'].append(results[i][name]['pars_err'][3])
            output_dict[name]['T'].append(results[i][name]['T'])
            output_dict[name]['Tpsf'].append(results[i][name]['Tpsf'])
            output_dict[name]['s2n'].append(results[i][name]['s2n'])
            output_dict[name]['flags'].append(results[i][name]['flags'])
            output_dict[name]['mcal_flags'].append(results[i]['mcal_flags'])

    return output_dict


def save_results(output_dict, output_name):
    """
    """

    f = io.FITSCatalog(output_name, open_mode=io.BaseCatalog.OpenMode.ReadWrite)

    for key in output_dict.keys():
        f.save_as_fits(output_dict[key], ext_name=key.upper())


def process(tile_cat_path, gal_vignet_path, bkg_vignet_path,
            psf_vignet_path, weight_vignet_path, flag_vignet_path):
    """
    """

    tile_cat = io.FITSCatalog(tile_cat_path, SEx_catalog=True)
    tile_cat.open()
    obj_id = np.copy(tile_cat.get_data()['NUMBER'])
    tile_vign = np.copy(tile_cat.get_data()['VIGNET'])
    tile_flag = np.copy(tile_cat.get_data()['FLAGS'])
    tile_imaflag = np.copy(tile_cat.get_data()['IMAFLAGS_ISO'])
    tile_cat.close()
    mask = (tile_flag == 0) & (tile_imaflag == 0)
    gal_vign_cat = io.FITSCatalog(gal_vignet_path)#, SEx_catalog=True)
    gal_vign_cat.open()
    bkg_vign_cat = io.FITSCatalog(bkg_vignet_path)#, SEx_catalog=True)
    bkg_vign_cat.open()
    psf_vign_cat = io.FITSCatalog(psf_vignet_path)#, SEx_catalog=True)
    psf_vign_cat.open()
    weight_vign_cat = io.FITSCatalog(weight_vignet_path)#, SEx_catalog=True)
    weight_vign_cat.open()
    flag_vignet_cat = io.FITSCatalog(flag_vignet_path)#, SEx_catalog=True)
    flag_vignet_cat.open()

    list_exp_name = gal_vign_cat.get_ext_name()
    hdu_ind = [i for i in range(len(list_exp_name)) if 'EPOCH' in list_exp_name[i]]

    skip_obj = False
    final_res = []
    prior = get_prior()
    for i_tile, id_tmp in enumerate(obj_id[mask][::2]):
        print("{}/{}".format(i_tile,len(obj_id[mask][::2])))
        gal_vign = []
        psf_vign = []
        sigma_psf = []
        weight_vign = []
        flag_vign = []
        counter = 0
        for hdu_index in hdu_ind:
            id_epoch = gal_vign_cat.get_data(hdu_index)['NUMBER']
            i_exp = np.where(i_tile == id_epoch)[0]
            if len(i_exp) == 0:
                if counter == 0:
                    skip_obj = True
                    break
                else:
                    continue
            tile_vign_tmp = tile_vign[i_tile]
            gal_vign_tmp = gal_vign_cat.get_data(hdu_index)['VIGNET'][i_exp[0]]
            bkg_vign_tmp = bkg_vign_cat.get_data(hdu_index)['VIGNET'][i_exp[0]]
            psf_vign_tmp = psf_vign_cat.get_data(hdu_index)['VIGNET'][i_exp[0]]
            weight_vign_tmp = weight_vign_cat.get_data(hdu_index)['VIGNET'][i_exp[0]]
            flag_vignet_tmp = flag_vignet_cat.get_data(hdu_index)['VIGNET'][i_exp[0]]
            gal_vign_sub_bkg = gal_vign_tmp - bkg_vign_tmp
            gal_vign.append(gal_vign_sub_bkg)
            psf_vign.append(psf_vign_tmp)
            sigma_psf.append(psf_vign_cat.get_data(hdu_index)['SIGMA_PSF_HSM'][i_exp[0]])
            weight_vign.append(weight_vign_tmp)
            flag_vign_add_mask = np.copy(flag_vignet_tmp)
            flag_vign_add_mask[np.where(tile_vign_tmp == -1e30)] = 2**10
            if len(np.where(flag_vign_add_mask.ravel() != 0)[0]) == (51*51):
                if counter == 0:
                    skip_obj = True
                    break
                else:
                    continue
            flag_vign.append(flag_vign_add_mask)
            counter += 1

        if skip_obj:
            skip_obj = False
            continue

        try:
            res = do_ngmix_metacal(gal_vign,
                                   psf_vign,
                                   sigma_psf,
                                   weight_vign,
                                   flag_vign,
                                   prior)
        except:
            continue
        res['obj_id'] = id_tmp
        res['n_epoch_model'] = len(gal_vign)
        final_res.append(res)
    gal_vign_cat.close()
    bkg_vign_cat.close()
    psf_vign_cat.close()
    weight_vign_cat.close()
    flag_vignet_cat.close()

    return final_res


@module_runner(input_module=['sextractor_runner', 'psfexinterp_runner', 'vignetmaker_runner'],
               version='1.0',
               file_pattern=['tile_sexcat', 'image', 'background', 'galaxy_psf', 'weight', 'flag'],
               file_ext=['.fits', '.fits', '.fits', '.fits', '.fits', '.fits'],
               depends=['numpy', 'ngmix'])
def ngmix_runner(input_file_list, output_dir, file_number_string,
                 config, w_log):

    output_name = output_dir + '/' + 'ngmix' + file_number_string + '.fits'

    metacal_res = process(*input_file_list)
    res_dict = compile_results(metacal_res)
    save_results(res_dict, output_name)

    return None, None
