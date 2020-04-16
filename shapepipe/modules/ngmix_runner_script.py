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

import ngmix
from ngmix.observation import Observation, ObsList, MultiBandObsList
from ngmix.fitting import LMSimple

import galsim

import sys
import os


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

    n_epoch = len(gals)

    if n_epoch == 0:
        raise ValueError("0 epoch to process")

    # Make observation
    gal_obs_list = ObsList()
    T_guess_psf = []
    for n_e in range(n_epoch):

        psf_jacob = ngmix.Jacobian(row=psfs[0].shape[0]/2.,
                                   col=psfs[0].shape[1]/2.,
                                   wcs=jacob_list[n_e])
        gal_jacob = ngmix.Jacobian(row=gals[0].shape[0]/2.,
                                   col=gals[0].shape[1]/2.,
                                   wcs=jacob_list[n_e])

        psf_obs = Observation(psfs[n_e], jacobian=psf_jacob)

        psf_T = 2. * psfs_sigma[n_e]**2.

        w = np.copy(weights[n_e])
        w[np.where(flags[n_e] != 0)] = 0.

        gal_obs = Observation(gals[n_e], weight=w, jacobian=gal_jacob,
                              psf=psf_obs)

        gal_obs_list.append(gal_obs)
        T_guess_psf.append(psf_T)

    boot = ngmix.bootstrap.MaxMetacalBootstrapper(gal_obs_list)
    psf_model = 'em3'
    gal_model = 'gauss'

    # metacal specific parameters
    metacal_pars = {'types': ['noshear', '1p', '1m', '2p', '2m'],
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
        'lm_pars': lm_pars}

    psf_pars = {'maxiter': 5000,
                'tol': 5.0e-6}

    # Tguess = np.mean(T_guess_psf)*0.186**2  # size guess in arcsec
    Tguess = 4.0*0.186**2
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

    return metacal_res


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
    names2 = ['id', 'n_epoch_model', 'g1', 'g1_err', 'g2', 'g2_err', 'T',
              'T_err', 'Tpsf', 's2n', 'flags', 'mcal_flags']
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
            output_dict[name]['T_err'].append(results[i][name]['T_err'])
            output_dict[name]['Tpsf'].append(results[i][name]['Tpsf'])
            output_dict[name]['s2n'].append(results[i][name]['s2n'])
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
            f_wcs_path):
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
    import pdb
    pdb.set_trace()
    for i_tile, id_tmp in enumerate(obj_id):

	print('MKDEBUG {}'.format(id_tmp))
        if id_tmp == 34:
		import pdb
		pdb.set_trace()
        # Preselection step
        # if (tile_flag[i_tile] > 1) or (tile_imaflag[i_tile] > 0):
        #     continue
        # if (sm[i_tile] + (5. / 3.) * sm_err[i_tile] < 0.01) and
        # (np.abs(sm[i_tile] + (5. / 3.) * sm_err[i_tile]) > 0.003):
        #     continue
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

            psf_vign.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])
            sigma_psf.append(psf_vign_cat[str(id_tmp)][expccd_name_tmp]['SHAPES']['SIGMA_PSF_HSM'])

            bkg_vign_tmp = bkg_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            gal_vign_sub_bkg = gal_vign_tmp - bkg_vign_tmp
            gal_vign.append(gal_vign_sub_bkg)

            weight_vign.append(weight_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET'])

            tile_vign_tmp = np.copy(tile_vign[i_tile])
            flag_vign_tmp = flag_vign_cat[str(id_tmp)][expccd_name_tmp]['VIGNET']
            flag_vign_tmp[np.where(tile_vign_tmp == -1e30)] = 2**10
            v_flag_tmp = flag_vign_tmp.ravel()
            if len(np.where(v_flag_tmp != 0)[0])/(51*51) > 1/3.:
                continue
            flag_vign.append(flag_vign_tmp)

            exp_name, ccd_n = re.split('-', expccd_name_tmp)
            jacob_list.append(get_jacob(f_wcs_file[exp_name][int(ccd_n)],
                                        tile_ra[i_tile],
                                        tile_dec[i_tile]))
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
        except Exception:
            print('ngmix fail on object {}'.format(id_tmp))
            continue
        res['obj_id'] = id_tmp
        res['n_epoch_model'] = len(gal_vign)
        final_res.append(res)

    f_wcs_file.close()
    gal_vign_cat.close()
    bkg_vign_cat.close()
    flag_vign_cat.close()
    weight_vign_cat.close()
    psf_vign_cat.close()

    return final_res


def ngmix_runner(input_file_list, output_name, file_number_string,
                 f_wcs_path):

    # output_name = (output_dir + '/' + 'ngmix' +
    #                file_number_string + '.fits')

    # f_wcs_path = config.getexpanded('NGMIX_RUNNER', 'LOG_WCS')

    metacal_res = process(*input_file_list, f_wcs_path)
    res_dict = compile_results(metacal_res)
    save_results(res_dict, output_name)

    return None, None


if __name__ == '__main__':

    argv = sys.argv

    list_id = int(argv[1])
    # print("list_id : {}\n{}".format(type(list_id), list_id))

    process_list = argv[2]
    # print("process_list : {}\n{}".format(type(process_list), process_list))
    f = np.load(process_list, mmap_mode='r')
    print("Input paths : \n{}\n".format(f[list_id]))

    f_wcs_path = argv[3]
    # print("wcs_path : {}\n{}".format(type(f_wcs_path), f_wcs_path))

    output_dir = argv[4]
    # print("output_dir : {}\n{}".format(type(output_dir), output_dir))
    output_name = (output_dir + '/' + 'ngmix' +
                   f[list_id][0] + '.fits')

    if not os.path.isfile(output_name):
        ngmix_runner(f[list_id][1:], output_name, f[list_id][0], f_wcs_path)
    else:
        print("File already exist : {}".format(output_name))