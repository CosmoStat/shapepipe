# -*- coding: utf-8 -*-

"""PSFEx MERGE STARCAT RUNNER

This module is used to merge the validation stars results of the PSFExinterp
runner. It is useful to calculate the pixel RMSE.

:Author: Tobias Liaudat

"""

import numpy as np
from astropy.io import fits
from shapepipe.pipeline import file_io as sc
from shapepipe.modules.module_decorator import module_runner
import re


@module_runner(input_module=['psfexinterp_runner'], version='1.0',
               file_pattern=['validation_psf'],
               file_ext=['.fits'], numbering_scheme = '-0000000-0',
               depends=['numpy', 're', 'galsim','astropy'],
               run_method='serial')
def psfex_merge_starcat_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):
    print('Merging validation results..')
    output_dir = run_dirs['output']


    x, y, ra, dec = [], [], [], []
    g1_psf, g2_psf, size_psf = [], [], []
    g1, g2, size = [], [], []
    flag_psf, flag_star = [], []
    mag, snr, psfex_acc = [], [], []
    ccd_nb = []
    pixel_mse, size_mse = [] , []
    pix_norm_mse, size_norm_mse = [] , []
    pix_filt_mse, size_filt_mse = [] , []

    bad_catalogs = 0

    for name in input_file_list:
        starcat_j = fits.open(name[0])
        stars = starcat_j[2].data['STAR_VIGNET']
        stars[stars<-1e6] = 0
        psfs = starcat_j[2].data['PSF_VIGNET']
        psfs[psfs<-1e6] = 0

        # Pixel mse calculation
        pix_val = np.sum((stars - psfs)**2)

        if pix_val < 1e20:
            # Normalised pixel mse calculation
            stars_norm_vals =  np.sqrt(np.sum(stars**2,axis=(1,2)))
            psfs_norm_vals =  np.sqrt(np.sum(psfs**2,axis=(1,2)))
            # Select non zero stars & psfs
            non_zero_elems = np.logical_and((psfs_norm_vals!=0),(stars_norm_vals!=0))
            # Calculate the filtered mse calculation
            pix_filt_val = np.sum((stars[non_zero_elems] -
                              psfs[non_zero_elems])**2)
            # Calculate the normalized (& filtered) mse calculation
            stars_norm_vals = stars_norm_vals[non_zero_elems].reshape(-1,1,1)
            psfs_norm_vals = psfs_norm_vals[non_zero_elems].reshape(-1,1,1)
            pix_norm_val = np.sum((stars[non_zero_elems]/stars_norm_vals -
                                   psfs[non_zero_elems]/psfs_norm_vals)**2)
            # Calculate sizes
            filt_size = stars[non_zero_elems].size
            regular_size = stars.size

            # Append the results to the lists
            pixel_mse.append(pix_val)
            size_mse.append(regular_size)
            pix_norm_mse.append(pix_norm_val)
            size_norm_mse.append(filt_size)
            pix_filt_mse.append(pix_filt_val)
            size_filt_mse.append(filt_size)

            # positions
            x += list(starcat_j[2].data['X'])
            y += list(starcat_j[2].data['Y'])
            try:
                ra += list(starcat_j[2].data['RA'])
                dec += list(starcat_j[2].data['DEC'])
            except:
                pass

            # shapes (convert sigmas to R^2)
            g1_psf += list(starcat_j[2].data['E1_PSF_HSM'])
            g2_psf += list(starcat_j[2].data['E2_PSF_HSM'])
            size_psf += list(starcat_j[2].data['SIGMA_PSF_HSM']**2)
            g1 += list(starcat_j[2].data['E1_STAR_HSM'])
            g2 += list(starcat_j[2].data['E2_STAR_HSM'])
            size += list(starcat_j[2].data['SIGMA_STAR_HSM']**2)

            # flags
            flag_psf += list(starcat_j[2].data['FLAG_PSF_HSM'])
            flag_star += list(starcat_j[2].data['FLAG_STAR_HSM'])

            # misc
            try:
                mag += list(starcat_j[2].data['MAG'])
                snr += list(starcat_j[2].data['SNR'])
            except:
                pass

            # CCD number
            try:
                ccd_nb += [re.split(r"\-([0-9]*)\-([0-9]+)\.", name[0])[-2]]*len(starcat_j[2].data['X'])
            except:
                ccd_nb += [99] # To know if there is some problem

        else:
            bad_catalogs += 1
            w_log.info('psfex_merge_starcat_runner: Bad catalog count = %d\n bad path = %s'%(bad_catalogs,name))



    # Regular pixel RMSE
    tot_pixel_rmse = np.sqrt(np.nansum(np.array(pixel_mse)) / np.nansum(np.array(size_mse)))
    print('Regular Total pixel RMSE = %.7e'%(tot_pixel_rmse))
    w_log.info('PSFEx_merge_starcat: Regular Total pixel RMSE = %.7e'%(tot_pixel_rmse))
    # Normalized pixel RMSE
    tot_pix_norm_rmse = np.sqrt(np.nansum(np.array(pix_norm_mse)) / np.nansum(np.array(size_norm_mse)))
    print('Normalized Total pixel RMSE = %.7e'%(tot_pix_norm_rmse))
    w_log.info('PSFEx_merge_starcat: Normalized Total pixel RMSE = %.7e'%(tot_pix_norm_rmse))
    # Normalized pixel RMSE
    tot_pix_filt_rmse = np.sqrt(np.nansum(np.array(pix_filt_mse)) / np.nansum(np.array(size_filt_mse)))
    print('Filtered Total pixel RMSE = %.7e'%(tot_pix_filt_rmse))
    w_log.info('PSFEx_merge_starcat: Filtered Total pixel RMSE = %.7e'%(tot_pix_filt_rmse))

    # Save files
    output = sc.FITSCatalog(output_dir + '/full_starcat-0000000.fits',
                            open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                            SEx_catalog=True)

    # convert back to sigma for consistency
    data = {'X': x, 'Y': y, 'RA': ra, 'DEC': dec,
            'E1_PSF_HSM': g1_psf, 'E2_PSF_HSM': g2_psf, 'SIGMA_PSF_HSM': np.sqrt(size_psf),
            'E1_STAR_HSM': g1, 'E2_STAR_HSM': g2, 'SIGMA_STAR_HSM': np.sqrt(size),
            'FLAG_PSF_HSM': flag_psf, 'FLAG_STAR_HSM': flag_star,
            'MAG': mag, 'SNR': snr, 'ACCEPTED': psfex_acc,
            'CCD_NB': ccd_nb}
    print('Writing full catalog...')
    output.save_as_fits(data, sex_cat_path=input_file_list[0][0])
    print('... Done.')

    return None, None
