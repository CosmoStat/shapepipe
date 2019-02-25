# -*- coding: utf-8 -*-

"""SPREAD MODEL RU+NNER

This module compute the spread model.

:Author: Axel Guinot

"""

import numpy as np
import galsim
from astropy import stats

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as io


def get_sm(obj_vign, psf_vign, model_vign, weight_vign):
    """
    """

    m = (obj_vign != -1e30) & (weight_vign>0.)
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
    """
    """

    gal_obj = galsim.Exponential(scale_radius=1./16.*sigma*2.634*pixel_scale, flux=flux)

    psf_obj = galsim.Gaussian(sigma=sigma*pixel_scale)

    gal_obj = galsim.Convolve(gal_obj, psf_obj)

    gal_vign = gal_obj.drawImage(nx=img_shape[0], ny=img_shape[1], scale=pixel_scale).array
    psf_vign = psf_obj.drawImage(nx=img_shape[0], ny=img_shape[1], scale=pixel_scale).array

    return gal_vign, psf_vign

def save_results(sex_cat_path, output_path, sm, sm_err):
    """
    """

    ori_cat = io.FITSCatalog(sex_cat_path, SEx_catalog=True)
    ori_cat.open()
    new_cat = io.FITSCatalog(output_path, SEx_catalog=True, open_mode=io.BaseCatalog.OpenMode.ReadWrite)
    ori_cat.add_col('SPREAD_MODEL', sm, new_cat=True, new_cat_inst=new_cat)
    ori_cat.close()
    new_cat.open()
    new_cat.add_col('SPREADERR_MODEL', sm_err)
    new_cat.close()


@module_runner(input_module=['sextractor_runner', 'psfexinterp_runner', 'vignetmaker_runner'], version='1.0',
               file_pattern=['sexcat', 'galaxy_psf', 'weight_vign'],
               file_ext=['.fits', '.fits', '.fits'],
               depends=['numpy', 'galsim'])
def spread_model_runner(input_file_list, output_dir, file_number_string,
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

    file_name = suffix + 'sexcat_sm' + file_number_string + '.fits'
    output_path = output_dir + '/' + file_name

    # Get data
    sex_cat = io.FITSCatalog(sex_cat_path, SEx_catalog=True)
    sex_cat.open()
    obj_id = np.copy(sex_cat.get_data()['NUMBER'])
    obj_vign = np.copy(sex_cat.get_data()['VIGNET'])
    sex_cat.close()

    psf_cat = io.FITSCatalog(psf_cat_path, SEx_catalog=True)
    psf_cat.open()
    ext_name = psf_cat.get_ext_name()
    hdu_ind = [i for i in range(len(ext_name)) if 'EPOCH' in ext_name[i]]
    dict_psf = []
    for i,i_h in enumerate(hdu_ind):
        dict_psf.append({})
        dict_psf[i]['id'] = psf_cat.get_data(i_h)['NUMBER']
        dict_psf[i]['sigma'] = psf_cat.get_data(i_h)['SIGMA_PSF_HSM']
    psf_cat.close()

    weight_cat = io.FITSCatalog(weight_cat_path, SEx_catalog=True)
    weight_cat.open()
    weigh_vign = weight_cat.get_data()['VIGNET']
    weight_cat.close()

    # Get spread model
    skip_obj = False
    spread_model_final = []
    spread_model_err_final = []
    for i in range(len(obj_id)):
        sigma_list = []
        obj_id_tmp = obj_id[i]
        for h in range(len(hdu_ind)):
            ind_tmp = np.where(dict_psf[h]['id'] == obj_id_tmp)[0]
            if len(ind_tmp) == 0:
                if h == 0:
                    skip_obj = True
                    break
                else:
                    continue
            sigma_list.append(dict_psf[h]['sigma'][ind_tmp])
        if skip_obj:
            skip_obj = False
            spread_model_final.append(-1)
            spread_model_err_final.append(1)
            continue
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

    save_results(sex_cat_path, output_path, spread_model_final, spread_model_err_final, )

    return None, None