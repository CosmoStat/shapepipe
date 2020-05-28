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

    bad_catalogs = 0

    for name in input_file_list:
        starcat_j = fits.open(name[0])

        # pixel mse calculation
        pix_val = np.sum((starcat_j[2].data['VIGNET_LIST'] - starcat_j[2].data['PSF_VIGNET_LIST'])**2)

        if pix_val < 1e15:
            pixel_mse.append(np.sum((starcat_j[2].data['VIGNET_LIST'] - starcat_j[2].data['PSF_VIGNET_LIST'])**2))
            size_mse.append(starcat_j[2].data['VIGNET_LIST'].size)

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


    # Pixel RMSE
    tot_pixel_rmse = np.sqrt(np.nansum(np.array(pixel_mse)) / np.nansum(np.array(size_mse)))
    print('Total pixel RMSE = %.5e'%(tot_pixel_rmse))
    w_log.info('MCCD_merge_starcat: Total pixel RMSE = %.5e'%(tot_pixel_rmse))

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
