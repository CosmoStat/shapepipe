# -*- coding: utf-8 -*-

"""MCCD MERGE STARCAT RUNNER

This module is used to merge the validation stars results of the MCCD runner.

:Author: Tobias Liaudat

"""

import numpy as np
import os
from astropy.io import fits
import sys
from shapepipe.pipeline import file_io as sc
from shapepipe.modules.module_decorator import module_runner


@module_runner(input_module=['mccd_runner'], version='1.0',
               file_pattern=['validation_psf'],
               file_ext=['.fits'], numbering_scheme = '-0000000',
               depends=['numpy', 'mccd_rca', 'galsim','astropy'],
               run_method='serial')
def mccd_merge_starcat_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):
    print('Merging validation results..')
    save_fullcat = True

    output_dir = run_dirs['output']

    x, y= [], []
    g1_psf, g2_psf, size_psf = [], [], []
    g1, g2, size = [], [], []
    flag_psf, flag_star = [], []
    ccd_nb = []
    pixel_mse = []
    size_mse = []
    pix_norm_mse, size_norm_mse = [] , []
    pix_filt_mse, size_filt_mse = [] , []

    bad_catalogs = 0

    for name in input_file_list:
        starcat_j = fits.open(name[0])

        # Pixel mse calculation
        pix_val = np.sum((starcat_j[2].data['VIGNET_LIST'] - starcat_j[2].data['PSF_VIGNET_LIST'])**2)

        if pix_val < 1e20:
            # Normalised pixel mse calculation
            stars_norm_vals =  np.sqrt(np.sum(starcat_j[2].data['VIGNET_LIST']**2,axis=(1,2)))
            psfs_norm_vals =  np.sqrt(np.sum(starcat_j[2].data['PSF_VIGNET_LIST']**2,axis=(1,2)))
            # Select non zero stars & psfs
            non_zero_elems = np.logical_and((psfs_norm_vals!=0),(stars_norm_vals!=0))
            # Calculate the filtered mse calculation
            pix_filt_val = np.sum((starcat_j[2].data['VIGNET_LIST'][non_zero_elems] -
                              starcat_j[2].data['PSF_VIGNET_LIST'][non_zero_elems])**2)
            # Calculate the normalized (& filtered) mse calculation
            stars_norm_vals = stars_norm_vals[non_zero_elems].reshape(-1,1,1)
            psfs_norm_vals = psfs_norm_vals[non_zero_elems].reshape(-1,1,1)
            pix_norm_val = np.sum((starcat_j[2].data['VIGNET_LIST'][non_zero_elems]/stars_norm_vals -
                                   starcat_j[2].data['PSF_VIGNET_LIST'][non_zero_elems]/psfs_norm_vals)**2)
            # Calculate sizes
            filt_size = starcat_j[2].data['VIGNET_LIST'][non_zero_elems].size
            regular_size = starcat_j[2].data['VIGNET_LIST'].size

            # Append the results to the lists
            pixel_mse.append(pix_val)
            size_mse.append(regular_size)
            pix_norm_mse.append(pix_norm_val)
            size_norm_mse.append(filt_size)
            pix_filt_mse.append(pix_filt_val)
            size_filt_mse.append(filt_size)


            # positions
            x += list(starcat_j[2].data['GLOB_POSITION_IMG_LIST'][:,0])
            y += list(starcat_j[2].data['GLOB_POSITION_IMG_LIST'][:,1])

            # shapes (convert sigmas to R^2)
            g1_psf += list(starcat_j[2].data['PSF_MOM_LIST'][:,0])
            g2_psf += list(starcat_j[2].data['PSF_MOM_LIST'][:,1])
            size_psf += list(starcat_j[2].data['PSF_MOM_LIST'][:,2]**2)
            g1 += list(starcat_j[2].data['STAR_MOM_LIST'][:,0])
            g2 += list(starcat_j[2].data['STAR_MOM_LIST'][:,1])
            size += list(starcat_j[2].data['STAR_MOM_LIST'][:,2]**2)

            # flags
            flag_psf += list(starcat_j[2].data['PSF_MOM_LIST'][:,3])
            flag_star += list(starcat_j[2].data['STAR_MOM_LIST'][:,3])

            # ccd id list
            ccd_nb += list(starcat_j[2].data['CCD_ID_LIST'])

        else:
            bad_catalogs += 1
            w_log.info('MCCD_merge_starcat: Bad catalog count = %d\n bad path = %s'%(bad_catalogs,name))


    # Regular pixel RMSE
    tot_pixel_rmse = np.sqrt(np.nansum(np.array(pixel_mse)) / np.nansum(np.array(size_mse)))
    print('Regular Total pixel RMSE = %.5e'%(tot_pixel_rmse))
    w_log.info('MCCD_merge_starcat: Regular Total pixel RMSE = %.5e'%(tot_pixel_rmse))
    # Normalized pixel RMSE
    tot_pix_norm_rmse = np.sqrt(np.nansum(np.array(pix_norm_mse)) / np.nansum(np.array(size_norm_mse)))
    print('Normalized Total pixel RMSE = %.5e'%(tot_pix_norm_rmse))
    w_log.info('MCCD_merge_starcat: Normalized Total pixel RMSE = %.5e'%(tot_pix_norm_rmse))
    # Normalized pixel RMSE
    tot_pix_filt_rmse = np.sqrt(np.nansum(np.array(pix_filt_mse)) / np.nansum(np.array(size_filt_mse)))
    print('Filtered Total pixel RMSE = %.5e'%(tot_pix_filt_rmse))
    w_log.info('MCCD_merge_starcat: Filtered Total pixel RMSE = %.5e'%(tot_pix_filt_rmse))

    output = sc.FITSCatalog(output_dir + '/full_starcat-0000000.fits',
                            open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                            SEx_catalog=True)
    # convert back to sigma for consistency
    data = {'X': x, 'Y': y,
            'E1_PSF_HSM': g1_psf, 'E2_PSF_HSM': g2_psf, 'SIGMA_PSF_HSM': np.sqrt(size_psf),
            'E1_STAR_HSM': g1, 'E2_STAR_HSM': g2, 'SIGMA_STAR_HSM': np.sqrt(size),
            'FLAG_PSF_HSM': flag_psf, 'FLAG_STAR_HSM': flag_star, 'CCD_NB': ccd_nb}
    print('Writing full catalog...')
    output.save_as_fits(data, sex_cat_path=input_file_list[0][0])
    print('... Done.')

    return None, None
