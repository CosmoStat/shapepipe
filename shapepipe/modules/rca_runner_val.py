# -*- coding: utf-8 -*-

"""RCA RUNNER

This file is the pipeline runner for the RCA package.

:Author: Morgan Schmitz

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
from astropy.io import fits
import numpy as np
import rca
try:
    import galsim.hsm as hsm
    from galsim import Image
    import_fail = False
except ImportError:
    import_fail = True

def handle_SExtractor_mask(stars, thresh):
    """ Reads SExtracted star stamps, generates RCA-compatible masks (that is, binary weights),
    and replaces bad pixels with 0s - they will not be used by RCA, but the ridiculous numerical
    values can otherwise still lead to problems because of convolutions."""
    mask = np.ones(stars.shape)
    mask[stars < thresh] = 0
    stars[stars < thresh] = 0
    return mask

def PolynomialA(starcat, pos_params):
    xs = starcat[2].data[pos_params[0]]
    ys = starcat[2].data[pos_params[1]]
    xxs = xs - (np.max(xs) - np.min(xs))/2
    yys = ys - (np.max(ys) - np.min(ys))/2

    VT = np.ones((6,len(xs)))
    VT[1] = xxs
    VT[2] = yys
    VT[3] = xxs*yys
    VT[4] = xxs**2
    VT[5] = yys**2

    alpha = np.eye(6)

    weight_norms = np.sqrt(np.sum(VT**2,axis=1))
    VT /= weight_norms.reshape(-1,1)

    return alpha,VT



def rca_fit(starcat, pos_params, rcainst_kw, rcafit_kw, output_dir, file_number_string, sex_thresh=-1e5):
    train_stars = rca.utils.rca_format(starcat[2].data['VIGNET'])
    train_mask = handle_SExtractor_mask(train_stars, thresh=sex_thresh)
    train_pos = np.array([[x, y] for x, y in
                          zip(starcat[2].data[pos_params[0]],
                              starcat[2].data[pos_params[1]])])

    rca_instance = rca.RCA(**rcainst_kw, verbose=True)
    S, A = rca_instance.fit(train_stars, train_pos, train_mask, **rcafit_kw)
    filename = output_dir + '/fitted_model'+file_number_string
    rca_instance.quicksave(filename)
    return rca_instance

def rca_transform(rca_instance, testcat, pos_params, output_dir, file_number_string):
    test_pos = np.array([[x, y] for x, y in
                          zip(testcat[2].data[pos_params[0]],
                              testcat[2].data[pos_params[1]])])
    PSFs = rca_instance.estimate_psf(test_pos)
    np.save(output_dir+'/test_psf'+file_number_string, PSFs)
    return PSFs

def rca_degrade(rca_instance, starcat, pos_params, sex_thresh=-1e5):
    test_stars = rca.utils.rca_format(starcat[2].data['VIGNET'])
    handle_SExtractor_mask(test_stars, thresh=sex_thresh)
    test_pos = np.array([[x, y] for x, y in
                          zip(starcat[2].data[pos_params[0]],
                              starcat[2].data[pos_params[1]])])
    return rca_instance.validation_stars(test_stars, test_pos)

def rca_validation(star_cat, PSFs, pos_params, star_cat_path, output_dir, file_number_string,
                   sex_thresh=-1e5):
    if import_fail:
        raise ImportError('GalSim is required for the RCA module to run in VALIDATION mode.')
    # extract star stamps and masks
    stars = np.copy(star_cat[2].data['VIGNET'])
    mask = handle_SExtractor_mask(stars, thresh=sex_thresh)

    # and some SExtractor parameters
    star_dict = {}
    #star_dict['NUMBER'] = np.copy(star_cat[2].data['NUMBER']) # [JB] to comment when in validation mode
    star_dict['X'] = np.copy(star_cat[2].data[pos_params[0]])
    star_dict['Y'] = np.copy(star_cat[2].data[pos_params[1]])
    #star_dict['RA'] = np.copy(star_cat[2].data['XWIN_WORLD']) # [JB] to comment when in validation mode
    #star_dict['DEC'] = np.copy(star_cat[2].data['YWIN_WORLD']) # [JB] to comment when in validation mode
    #star_dict['MAG'] = np.copy(star_cat[2].data['MAG_AUTO']) # [JB] to comment when in validation mode
    #star_dict['SNR'] = np.copy(star_cat[2].data['SNR_WIN']) # [JB] to comment when in validation mode

    # compute star shapes with HSM
    badpix_mask = np.abs(mask-1) # hsm thinks 0 means good
    star_moms = [hsm.FindAdaptiveMom(Image(star), badpix=Image(bp), strict=False)
                 for star,bp in zip(stars,badpix_mask)]
    star_shapes = np.array([[moms.observed_shape.g1,
                             moms.observed_shape.g2,
                             moms.moments_sigma,
                             int(bool(moms.error_message))]
                        for moms in star_moms])

    star_dict['E1_STAR_HSM'] = star_shapes[:,0]
    star_dict['E2_STAR_HSM'] = star_shapes[:,1]
    star_dict['SIGMA_STAR_HSM'] = star_shapes[:,2]
    star_dict['FLAG_STAR_HSM'] = star_shapes[:,3]

    # and PSF shapes
    psf_moms = [hsm.FindAdaptiveMom(Image(psf), strict=False) for psf in PSFs]
    psf_shapes = np.array([[moms.observed_shape.g1,
                             moms.observed_shape.g2,
                             moms.moments_sigma,
                             int(bool(moms.error_message))]
                        for moms in psf_moms])

    star_dict['E1_PSF_HSM'] = psf_shapes[:,0]
    star_dict['E2_PSF_HSM'] = psf_shapes[:,1]
    star_dict['SIGMA_PSF_HSM'] = psf_shapes[:,2]
    star_dict['FLAG_PSF_HSM'] = psf_shapes[:,3]

    # and save in scatalog format
    filename = output_dir + '/validation_psf'+file_number_string+'.fits'
    output = sc.FITSCatalog(filename, open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                            SEx_catalog=True)
    output.save_as_fits(star_dict, sex_cat_path=star_cat_path)


@module_runner(input_module=['setools_runner'], version='1.0',
               file_pattern=['star_selection'],
               file_ext=['.fits'],
               depends=['numpy', 'rca', 'galsim'])
def rca_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):
    mode = config.get('RCA_RUNNER', 'MODE')
    sex_thresh = config.getfloat('RCA_RUNNER', 'SEXMASK_THRESH')
    pos_params = config.getlist('RCA_RUNNER', 'POSITION_PARAMS')
    star_thresh = config.getint('RCA_RUNNER', 'STAR_THRESH')
    if mode in ['FIT', 'FIT_TRANSFORM']:
        n_comp = config.getint('RCA_RUNNER', 'N_COMP')
        psf_size = config.getfloat('RCA_RUNNER', 'PSF_SIZE')
        n_eigenvects = config.getint('RCA_RUNNER', 'N_EIGENVECTS')
        ksig = config.getfloat('RCA_RUNNER', 'KSIG')
        filt_path = config.get('RCA_RUNNER', 'FILTER_PATH')
        filters = None if (filt_path == 'None') else np.load(filt_path)
        alphapath = config.get('RCA_RUNNER', 'ALPHA')

    if mode == 'FIT':
        starcat_path = input_file_list[0]
        starcat = fits.open(starcat_path)
        if len(starcat[2].data) < star_thresh:
            w_log.info('Star catalog {} rejected because it contains too few stars.'.format(starcat_path))
            return None, None
        if alphapath == 'PSFEx':
            alpha, VT = PolynomialA(starcat, pos_params)
            n_comp = 6
        else:
            alpha, VT = None, None
        rcainst_kw = {'n_comp': n_comp, 'filters': filters, 'ksig': ksig}
        rcafit_kw = {'alpha': alpha, 'VT': VT, 'psf_size': psf_size, 'n_eigenvects': n_eigenvects}
        rca_fit(starcat, pos_params, rcainst_kw, rcafit_kw, run_dirs['output'], file_number_string, sex_thresh)

    elif mode == 'TRANSFORM':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In TRANSFORM mode, both RCA outputs (as .npy) and catalogs (as .fits) are expected.')
        rca_path = input_file_list[0]
        test_path = input_file_list[1]
        testcat = fits.open(test_path)
        rca_instance = rca.quickload(rca_path)
        rca_transform(rca_instance, testcat, pos_params, run_dirs['output'], file_number_string)

    elif mode == 'FIT_TRANSFORM':
        if len(input_file_list) < 2:
            raise ValueError('In FIT_TRANSFORM mode, two catalogs (as .fits) are expected.')
        starcat_path, test_path = input_file_list[:2]
        starcat = fits.open(starcat_path)
        if len(starcat[2].data) < star_thresh:
            w_log.info('Star catalog {} rejected because it contains too few stars.'.format(starcat_path))
            return None, None
        testcat = fits.open(test_path)

        rca_instance = rca_fit(starcat, pos_params, rcainst_kw, rcafit_kw,
                               run_dirs['output'], file_number_string, sex_thresh)
        rca_transform(rca_instance, testcat, pos_params, run_dirs['output'], file_number_string)

    elif mode == 'VALIDATION':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In VALIDATION mode, both RCA outputs (as .npy) and star catalogs (as .fits) are expected.')
        apply_degradation = config.getboolean('RCA_RUNNER', 'APPLY_DEGRADATION')
        rca_path = input_file_list[0]
        test_path = input_file_list[1]
        testcat = fits.open(test_path)
        rca_instance = rca.quickload(rca_path)
        if apply_degradation:
            PSFs = rca_degrade(rca_instance, testcat, pos_params)
        else:
            PSFs = rca_transform(rca_instance, testcat, pos_params, run_dirs['output'], file_number_string)
        rca_validation(testcat, PSFs, pos_params, test_path, run_dirs['output'], file_number_string,
                       sex_thresh)

    else:
        raise ValueError('MODE should be in ["FIT", "TRANSFORM", "FIT_TRANSFORM", "VALIDATION"].')

    return None, None
