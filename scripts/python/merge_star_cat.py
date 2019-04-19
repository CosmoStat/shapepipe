import numpy as np
import os
from astropy.io import fits
import re
from shapepipe.pipeline import file_io as sc

starcats_path = '/home/guinot/ShapePipe_dir/output_singles2/shapepipe_run_2019-03-29_15-46-39/psfexinterp_runner/output/'
stub = 'validation_psf'
save_fullcat = True

# Read starcats and save full
starcat_names = os.listdir(starcats_path)
starcat_names = [starcats_path+name for name in starcat_names if stub in name]

x, y, ra, dec = [], [], [], []
g1_psf, g2_psf, size_psf = [], [], []
g1, g2, size = [], [], []
flag_psf, flag_star = [], []
mag, snr, psfex_acc = [], [], []
ccd_nb = []


for name in starcat_names:
    starcat_j = fits.open(name)

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
    ccd_nb += [re.split(r"\-([0-9]*)\-([0-9]+)\.", name)[-2]]*len(starcat_j[2].data['RA'])

if save_fullcat:
    output = sc.FITSCatalog(starcats_path+'full_starcat.fits',
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
    output.save_as_fits(data, sex_cat_path=starcat_names[0])
    print('... Done.')
