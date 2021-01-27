# -*- coding: utf-8 -*-

"""CREATE LOG EXP HEADER

This module merges the PSF validation catalogs from the PSFExInterpolation module
to run the statistics on them.

:Author: Morgan Schmitz, Axel Guinot

"""


import numpy as np
from astropy.io import fits

import os
import re

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc


@module_runner(input_module='psfex_interp_runner', version='1.0',
               file_pattern=['validation_psf'],
               file_ext=['.fits'], depends=['numpy', 'sqlitedict'],
               run_method='serial')
def merge_star_cat_runner(input_file_list, run_dirs, file_number_string,
                          config, w_log):

    output_dir = run_dirs['output']
    if config.has_option('MERGE_STAR_CAT_RUNNER', 'OUTPUT_PATH'):
        output_dir = config.getexpanded('MERGE_STAR_CAT_RUNNER', 'OUTPUT_PATH')

    x, y, ra, dec = [], [], [], []
    g1_psf, g2_psf, size_psf = [], [], []
    g1, g2, size = [], [], []
    flag_psf, flag_star = [], []
    mag, snr, psfex_acc = [], [], []
    ccd_nb = []

    for name in input_file_list:
        starcat_j = fits.open(name[0])

        # positions
        x += list(starcat_j[2].data['X'])
        y += list(starcat_j[2].data['Y'])
        ra += list(starcat_j[2].data['RA'])
        dec += list(starcat_j[2].data['DEC'])

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
        mag += list(starcat_j[2].data['MAG'])
        snr += list(starcat_j[2].data['SNR'])
        psfex_acc += list(starcat_j[2].data['ACCEPTED'])

        # CCD number
        ccd_nb += [int(re.split(r"\-([0-9]*)\-([0-9]+)\.", name[0])[-2])]*len(starcat_j[2].data['RA'])

    output = sc.FITSCatalog(output_dir + '/full_starcat.fits',
                            open_mode=sc.BaseCatalog.OpenMode.ReadWrite)
    # convert back to sigma for consistency
    data = {'X': x, 'Y': y, 'RA': ra, 'DEC': dec,
            'E1_PSF_HSM': g1_psf, 'E2_PSF_HSM': g2_psf, 'SIGMA_PSF_HSM': np.sqrt(size_psf),
            'E1_STAR_HSM': g1, 'E2_STAR_HSM': g2, 'SIGMA_STAR_HSM': np.sqrt(size),
            'FLAG_PSF_HSM': flag_psf, 'FLAG_STAR_HSM': flag_star,
            'MAG': mag, 'SNR': snr, 'ACCEPTED': psfex_acc,
            'CCD_NB': ccd_nb}
    print('Writing full catalog...')
    output.save_as_fits(data)
    print('... Done.')

    return None, None
