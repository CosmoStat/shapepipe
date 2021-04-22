# -*- coding: utf-8 -*-

"""SPREAD MODEL RU+NNER

This module compute the spread model.

:Author: Axel Guinot

"""

import numpy as np
import galsim

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io
from sqlitedict import SqliteDict


def get_sm(obj_vign, psf_vign, model_vign, weight_vign):
    """ Get Spread model

    This method compute the spread moel for an object.

    Parameters
    ----------
    obj_vign : numpy.ndarray
        Vignet of the object.
    psf_vign : numpy.ndarray
        Vignet of the gaussian model of the PSF.
    model_vign : numpy.ndarray
        Vignet of the galaxy model.
    weight_vign : numpy.ndarray
        Vignet of the weight at the object position.

    Returns
    -------
    sm : float
        Spread model value.
    sm_err : float
        Spread model error value.

    """

    m = (obj_vign > -1e29) & (weight_vign > 0.)
    w = m.astype(float)

    t_v = model_vign.ravel()
    g_v = obj_vign.ravel()
    psf_v = psf_vign.ravel()

    noise_v = (1. / weight_vign).ravel()
    noise_v[np.isinf(noise_v)] = 0.
    w_v = w.ravel()

    tg = np.sum(t_v*w_v*g_v)
    pg = np.sum(psf_v*w_v*g_v)
    tp = np.sum(t_v*w_v*psf_v)
    pp = np.sum(psf_v*w_v*psf_v)

    tnt = np.sum(t_v*noise_v*t_v*w_v)
    pnp = np.sum(psf_v*noise_v*psf_v*w_v)
    tnp = np.sum(t_v*noise_v*psf_v*w_v)
    err = tnt * pg**2. + pnp * tg**2. - 2. * tnp * pg * tg

    if pg > 0.:
        sm = (tg / pg) - (tp / pp)
    else:
        sm = 1.

    if (pg > 0.) & (err > 0.):
        sm_err = np.sqrt(err) / pg**2.
    else:
        sm_err = 1.

    return sm, sm_err


def get_model(sigma, flux, img_shape, pixel_scale=0.186):
    """ Get model

    This method compute :
     - an exponential galaxy model with scale radius = 1/16 FWHM
     - a gaussian model for the PSF

    Parameters
    ----------
    sigma : float
        Sigma of the PSF (in pixel units).
    flux : float
        Flux of the galaxy for the model.
    img_shape : list
        Size of the output vignet [xsize, ysize].
    pixel_scale : float
        Pixel scale to use for the model (in arcsec/pixel)

    Returns
    -------
    gal_vign : numpy.ndarray
        Vignet of the galaxy model.
    psf_vign : numpy.ndarray
        Vignet of the PSF model.

    """

    gal_obj = galsim.Exponential(scale_radius=1./16.*sigma*2.355*pixel_scale,
                                 flux=flux)

    psf_obj = galsim.Gaussian(sigma=sigma*pixel_scale)

    gal_obj = galsim.Convolve(gal_obj, psf_obj)

    gal_vign = gal_obj.drawImage(nx=img_shape[0], ny=img_shape[1],
                                 scale=pixel_scale).array
    psf_vign = psf_obj.drawImage(nx=img_shape[0], ny=img_shape[1],
                                 scale=pixel_scale).array

    return gal_vign, psf_vign


def save_results(sex_cat_path, output_path, sm, sm_err, mag, number,
                 mode='new'):
    """ Save results

    Save the spread model results.

    Parameters
    ----------
    sex_cat_path : str
        Path of the original SExtractor catalog.
    output_path : str
        Path of the output catalog.
    sm : numpy.ndarray
        Value of the spread model for all objects.
    sm_err : numpy.ndarray
        Value of the spread model error for all objects.
    mag : numpy.ndarray
        Magnitude of all objects (only for new catalog).
    number : numpy.ndarray
        Id of all objects (only for new catalog).
    mode : str
        Must be in ['new', 'add'].
        'new' will create a new catalog with : [number, mag, sm, sm_err]
        'add' will output a copy of the input SExtractor with the column sm
        and sm_err.

    """

    if mode == 'new':
        new_cat = io.FITSCatalog(output_path, SEx_catalog=True,
                                 open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        dict_data = {'NUMBER': number,
                     'MAG': mag,
                     'SPREAD_MODEL': sm,
                     'SPREADERR_MODEL': sm_err}
        new_cat.save_as_fits(data=dict_data, sex_cat_path=sex_cat_path)
    elif mode == 'add':
        ori_cat = io.FITSCatalog(sex_cat_path, SEx_catalog=True)
        ori_cat.open()
        new_cat = io.FITSCatalog(output_path, SEx_catalog=True,
                                 open_mode=io.BaseCatalog.OpenMode.ReadWrite)
        ori_cat.add_col('SPREAD_MODEL', sm, new_cat=True, new_cat_inst=new_cat)
        ori_cat.close()
        new_cat.open()
        new_cat.add_col('SPREADERR_MODEL', sm_err)
        new_cat.close()
    else:
        ValueError('Mode must be in [new, add].')


@module_runner(input_module=['sextractor_runner', 'psfex_interp_runner_me',
                             'vignetmaker_runner'], version='1.0',
               file_pattern=['sexcat', 'galaxy_psf', 'weight_vign'],
               file_ext=['.fits', '.sqlite', '.fits'],
               depends=['numpy', 'galsim'])
def spread_model_runner(input_file_list, run_dirs, file_number_string,
                        config, w_log):

    sex_cat_path, psf_cat_path, weight_cat_path = input_file_list

    if config.has_option('SPREAD_MODEL_RUNNER', 'SUFFIX'):
        suffix = config.get('SPREAD_MODEL_RUNNER', 'SUFFIX')
        if (suffix.lower() != 'none') & (suffix != ''):
            suffix = suffix + '_'
        else:
            suffix = ''
    else:
        suffix = ''

    pixel_scale = config.getfloat('SPREAD_MODEL_RUNNER', 'PIXEL_SCALE')
    output_mode = config.get('SPREAD_MODEL_RUNNER', 'OUTPUT_MODE')

    file_name = suffix + 'sexcat_sm' + file_number_string + '.fits'
    output_path = run_dirs['output'] + '/' + file_name

    # Get data
    sex_cat = io.FITSCatalog(sex_cat_path, SEx_catalog=True)
    sex_cat.open()
    obj_id = np.copy(sex_cat.get_data()['NUMBER'])
    obj_vign = np.copy(sex_cat.get_data()['VIGNET'])
    obj_mag = None
    if output_mode == 'new':
        obj_mag = np.copy(sex_cat.get_data()['MAG_AUTO'])
    sex_cat.close()

    # psf_cat = np.load(psf_cat_path, allow_pickle=True).item()
    psf_cat = SqliteDict(psf_cat_path)

    weight_cat = io.FITSCatalog(weight_cat_path, SEx_catalog=True)
    weight_cat.open()
    weigh_vign = weight_cat.get_data()['VIGNET']
    weight_cat.close()

    # Get spread model
    skip_obj = False
    spread_model_final = []
    spread_model_err_final = []
    for i, id_tmp in enumerate(obj_id):
        sigma_list = []

        if psf_cat[str(id_tmp)] == 'empty':
            spread_model_final.append(-1)
            spread_model_err_final.append(1)
            continue

        psf_expccd_name = list(psf_cat[str(id_tmp)].keys())

        for expccd_name_tmp in psf_expccd_name:
            sigma_list.append(psf_cat[str(id_tmp)][expccd_name_tmp]['SHAPES']['SIGMA_PSF_HSM'])

        obj_vign_tmp = obj_vign[i]
        obj_flux_tmp = 1.
        obj_sigma_tmp = np.mean(sigma_list)
        obj_weight_tmp = weigh_vign[i]
        obj_model_tmp, obj_psf_tmp = get_model(obj_sigma_tmp,
                                               obj_flux_tmp,
                                               obj_vign_tmp.shape,
                                               pixel_scale)

        obj_sm, obj_sm_err = get_sm(obj_vign_tmp,
                                    obj_psf_tmp,
                                    obj_model_tmp,
                                    obj_weight_tmp)

        spread_model_final.append(obj_sm)
        spread_model_err_final.append(obj_sm_err)

    spread_model_final = np.array(spread_model_final)
    spread_model_err_final = np.array(spread_model_err_final)

    psf_cat.close()

    save_results(sex_cat_path, output_path, spread_model_final,
                 spread_model_err_final,
                 obj_mag,
                 obj_id,
                 output_mode)

    return None, None
