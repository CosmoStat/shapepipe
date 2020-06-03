# -*- coding: utf-8 -*-

"""MCCD RUNNER

This file is the pipeline runner for the MCCD package.

:Author: Tobias Liaudat

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
from astropy.io import fits
import numpy as np
import mccd_rca.my_mccd_rca as mccd
import mccd_rca.utils as utils
import mccd_rca.mccd_utils as mccd_utils
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

def mccd_rca_fit(starcat, rcainst_kw, rcafit_kw, output_dir, file_number_string, sex_thresh=-1e5):

    # Extract the catalog and build the list with each ccd data
    # Stars are already masked and masks are provided
    my_ccd_list = starcat[2].data['CCD_ID_LIST'].astype(int)
    ccd_unique_list = np.unique(starcat[2].data['CCD_ID_LIST']).astype(int)

    positions = starcat[2].data['GLOB_POSITION_IMG_LIST']
    stars = starcat[2].data['VIGNET_LIST']
    ccds = starcat[2].data['CCD_ID_LIST']
    masks = starcat[2].data['MASK_LIST']

    # If masks are not provided they have to be calculated
    if np.sum(masks) == 0:
        masks = handle_SExtractor_mask(stars, sex_thresh)

    # SNR treatment
    try:
        SNRs = starcat[2].data['SNR_WIN_LIST']
        SNR_weights =  SNRs/np.median(SNRs)             # Strategy N5
        # SNR_weights[SNRs<50] = SNR_weights[SNRs<50]/10. #
        SNR_weights[SNR_weights>2.] = 2.
        SNR_weights[SNR_weights<0.1] = 0.1 # [TL] 0.01
        SNR_weight_list = [SNR_weights[my_ccd_list==ccd] for ccd in ccd_unique_list]
    except:
        SNR_weight_list = None
        print('No SNR weights are being used.')

    pos_list = [positions[my_ccd_list==ccd] for ccd in ccd_unique_list]
    star_list = [utils.rca_format(stars[my_ccd_list==ccd]) for ccd in ccd_unique_list]
    mask_list = [utils.rca_format(masks[my_ccd_list==ccd]) for ccd in ccd_unique_list]

    ccd_list = [ccds[my_ccd_list==ccd].astype(int) for ccd in ccd_unique_list]
    ccd_list= [np.unique(_list)[0].astype(int) for _list in ccd_list]

    rca_instance = mccd.mccd_RCA(**rcainst_kw, verbose=True)
    S, A_loc, A_glob, alpha, pi = rca_instance.fit(star_list, pos_list, ccd_list,
                                    mask_list, SNR_weight_list, **rcafit_kw)

    filename = output_dir + '/fitted_model'+file_number_string
    rca_instance.quicksave(filename)
    return rca_instance

def mccd_validation(rca_path, test_path, apply_degradation = True, mccd_debug = False, sex_thresh=-1e5):
    # Principal validation for mccd catalogs

    # Import
    rca_instance = mccd.mccd_quickload(rca_path)
    testcat = fits.open(test_path)

    # Saving file dictionary
    star_dict = {}

    # Preparing data in ccd-list format
    # Loading
    my_ccd_list = testcat[2].data['CCD_ID_LIST'].astype(int)
    ccd_unique_list = np.unique(testcat[2].data['CCD_ID_LIST']).astype(int)

    positions = testcat[2].data['GLOB_POSITION_IMG_LIST']
    stars = testcat[2].data['VIGNET_LIST']
    masks = testcat[2].data['MASK_LIST']
    ccds = testcat[2].data['CCD_ID_LIST']

    # If masks are not provided they have to be calculated
    if masks[0] == False:
        masks = handle_SExtractor_mask(stars, sex_thresh)

    # Save test stars
    # star_dict['GLOB_POSITION_IMG_LIST'] = np.copy(positions)
    # star_dict['VIGNET_LIST'] = np.copy(stars)
    # star_dict['MASK_LIST'] = np.copy(masks)
    # star_dict['CCD_ID_LIST'] = np.copy(ccds)

    val_pos_list = [positions[my_ccd_list==ccd] for ccd in ccd_unique_list]
    val_star_list = [utils.rca_format(stars[my_ccd_list==ccd]) for ccd in ccd_unique_list]
    val_mask_list = [utils.rca_format(masks[my_ccd_list==ccd]) for ccd in ccd_unique_list]
    val_ccd_list = [ccds[my_ccd_list==ccd].astype(int) for ccd in ccd_unique_list]
    val_ccd_list_to_save = np.copy(val_ccd_list)
    val_ccd_list= [np.unique(_list)[0].astype(int) for _list in val_ccd_list]


    if apply_degradation:
        if mccd_debug:
            PSF_list = []
            PSF_glob_list = []
            PSF_loc_list = []

            for it in range(len(val_star_list)):

                deg_PSFs, deg_PSFs_glob, deg_PSFs_loc = rca_instance.validation_stars(
                    val_star_list[it], val_pos_list[it], val_mask_list[it],
                    val_ccd_list[it],mccd_debug)
                PSF_list.append(utils.rca_format(deg_PSFs))
                PSF_glob_list.append(deg_PSFs_glob)
                PSF_loc_list.append(deg_PSFs_loc)

            star_dict['PSF_GLOB_VIGNET_LIST'] = np.copy(np.concatenate(PSF_glob_list,axis=0))
            star_dict['PSF_LOC_VIGNET_LIST'] = np.copy(np.concatenate(PSF_loc_list,axis=0))
        else:
            PSF_list = [rca_instance.validation_stars(_star, _pos, _mask, _ccd_id, mccd_debug)
                        for _star,_pos,_mask,_ccd_id in
                        zip(val_star_list,val_pos_list,val_mask_list,val_ccd_list)]
            # Have the PSFs in rca format, to match the stars
            PSF_list = [utils.rca_format(psfs) for psfs in PSF_list]

    # Calculate pixel-MSE
    myMSE = [np.mean((psfs - stars)**2) for psfs,stars in zip(PSF_list,val_star_list)]

    # Calculate moments
    n_ccds = len(val_star_list)
    star_shapes_list = []
    psf_shapes_list = []
    psf_vignet_list = []
    pos_list = []
    star_list = []
    mask_list = []
    ccd_list = []

    for it in range(n_ccds):
        # for each ccd
        test_stars = val_star_list[it]
        badpix_mask = np.rint(np.abs(val_mask_list[it]-1)) # hsm thinks 0 means good
        matched_psfs = PSF_list[it]

        # Stars
        star_moms = [hsm.FindAdaptiveMom(Image(star), badpix=Image(bp), strict=False)
                     for star,bp in zip(utils.reg_format(test_stars),utils.reg_format(badpix_mask))]
        star_shapes = np.array([[moms.observed_shape.g1,
                                 moms.observed_shape.g2,
                                 moms.moments_sigma,
                                 int(bool(moms.error_message))]
                            for moms in star_moms])
        # PSFs
        psf_moms = [hsm.FindAdaptiveMom(Image(psf), strict=False) for psf in utils.reg_format(matched_psfs)]
        psf_shapes = np.array([[moms.observed_shape.g1,
                                 moms.observed_shape.g2,
                                 moms.moments_sigma,
                                 int(bool(moms.error_message))]
                            for moms in psf_moms])

        star_shapes_list.append(star_shapes)
        psf_shapes_list.append(psf_shapes)
        psf_vignet_list.append(utils.reg_format(matched_psfs))
        pos_list.append(val_pos_list[it])
        star_list.append(utils.reg_format(val_star_list[it]))
        mask_list.append(utils.reg_format(val_mask_list[it]))
        ccd_list.append(val_ccd_list_to_save[it])

    # Prepare the PSF list and the moments in an array form
    # To be able to save them in the fits format
    psf_shapes = np.concatenate(psf_shapes_list,axis=0)
    star_shapes = np.concatenate(star_shapes_list,axis=0)
    psf_vignets = np.concatenate(psf_vignet_list,axis=0)
    pos_ordered = np.concatenate(pos_list ,axis=0)
    star_ordered = np.concatenate(star_list ,axis=0)
    mask_ordered = np.concatenate(mask_list ,axis=0)
    ccd_ordered = np.concatenate(ccd_list ,axis=0)

    # Save the results and the psfs
    star_dict['PSF_VIGNET_LIST'] = np.copy(psf_vignets)
    star_dict['PSF_MOM_LIST'] = np.copy(psf_shapes)
    star_dict['STAR_MOM_LIST'] = np.copy(star_shapes)
    star_dict['GLOB_POSITION_IMG_LIST'] = np.copy(pos_ordered)
    star_dict['VIGNET_LIST'] = np.copy(star_ordered)
    star_dict['MASK_LIST'] = np.copy(mask_ordered)
    star_dict['CCD_ID_LIST'] = np.copy(ccd_ordered)

    return star_dict

def mccd_response(rca_path, grid_xy, apply_degradation = True,mccd_debug = False):
    # Response for the mccd models
    rca_instance = mccd.mccd_quickload(rca_path)

    ccd_unique_list = np.copy(rca_instance.ccd_list)

    # Saving file dictionary
    star_dict = {}

    loc2glob = mccd_utils.Loc2Glob()
    im_dim = rca_instance.S[0].shape[0]

    # Generate local generic grid
    x_lin = np.linspace(start = im_dim, stop = loc2glob.x_npix - im_dim, num=grid_xy[0])
    y_lin = np.linspace(start = im_dim, stop = loc2glob.y_npix - im_dim, num=grid_xy[1])
    xv, yv = np.meshgrid(x_lin, y_lin)
    x_coor = xv.flatten()
    y_coor = yv.flatten()

    position_list = []
    ccd_list = []

    for it in range(len(ccd_unique_list)):

        x_glob , y_glob = loc2glob.loc2glob_img_coord(ccd_n = ccd_unique_list[it],
                                                        x_coor=np.copy(x_coor), y_coor=np.copy(y_coor))
        position_list.append(np.array([x_glob,y_glob]).T)
        ccd_list.append((np.ones(len(x_glob),dtype=int)*ccd_unique_list[it]).astype(int))


    if apply_degradation:
        PSF_list = [rca_instance.validation_stars(None, _pos, test_masks=None,
                        ccd_id=_ccd_id.astype(int),mccd_debug=mccd_debug,response_flag=True)
                    for _pos,_ccd_id in zip(position_list,ccd_unique_list)]


    # Calculate moments
    n_ccds = len(ccd_unique_list)

    psf_shapes_list = []

    for it in range(n_ccds):
        # for each ccd
        matched_psfs = PSF_list[it]

        # PSFs
        psf_moms = [hsm.FindAdaptiveMom(Image(psf), strict=False) for psf in matched_psfs]
        psf_shapes = np.array([[moms.observed_shape.g1,
                                 moms.observed_shape.g2,
                                 moms.moments_sigma,
                                 int(bool(moms.error_message))]
                            for moms in psf_moms])

        psf_shapes_list.append(psf_shapes)

    # Prepare the PSF list and the moments in an array form
    # To be able to save them in the fits format
    psf_shapes = np.concatenate(psf_shapes_list,axis=0)
    psf_positions = np.concatenate(position_list,axis=0)
    psf_ccd_id = np.concatenate(ccd_list,axis=0)

    # Save the results and the psfs
    star_dict['PSF_MOM_LIST'] = np.copy(psf_shapes)
    star_dict['GLOB_POSITION_IMG_LIST'] = np.copy(psf_positions)
    star_dict['CCD_ID_LIST'] = np.copy(psf_ccd_id)

    return star_dict

def fit(input_file_list, run_dirs, file_number_string, config, w_log):
    # Get config variables
    sex_thresh = config.getfloat('MCCD', 'SEXMASK_THRESH')
    star_thresh = config.getint('MCCD', 'STAR_THRESH')
    n_comp_loc = config.getint('MCCD', 'N_COMP_LOC')
    d_comp_glob = config.getint('MCCD', 'D_COMP_GLOB')
    psf_size = config.getfloat('MCCD', 'PSF_SIZE')
    psf_size_type = config.get('MCCD', 'PSF_SIZE_TYPE')
    n_eigenvects = config.getint('MCCD', 'N_EIGENVECTS')
    if n_eigenvects == 0:
        n_eigenvects = None
    ksig_loc = config.getfloat('MCCD', 'KSIG_LOC')
    ksig_glob = config.getfloat('MCCD', 'KSIG_GLOB')
    n_iter_rca = config.getint('MCCD', 'N_ITER_RCA') # [TL] modif
    nb_iter_glob = config.getint('MCCD', 'N_ITER_GLOB') # [TL] modif
    nb_iter_loc = config.getint('MCCD', 'N_ITER_LOC') # [TL] modif
    nb_subiter_S_loc = config.getint('MCCD', 'NB_SUBITER_S_LOC') # [TL] modif
    nb_subiter_A_loc = config.getint('MCCD', 'NB_SUBITER_A_LOC') # [TL] modif
    nb_subiter_S_glob = config.getint('MCCD', 'NB_SUBITER_S_GLOB') # [TL] modif
    nb_subiter_A_glob = config.getint('MCCD', 'NB_SUBITER_A_GLOB') # [TL] modif
    loc_model = config.get('MCCD', 'LOC_MODEL')
    filt_path = config.get('MCCD', 'FILTER_PATH')
    filters = None if (filt_path == 'None') else np.load(filt_path)

    # Prepare inputs to run the main fit function
    starcat_path = input_file_list[0] # train path
    starcat = fits.open(starcat_path)
    if len(starcat[2].data) < star_thresh:
        w_log.info('Star catalog {} rejected because it contains too few stars.'.format(starcat_path))
        return None, None

    rcainst_kw = {'n_comp_loc': n_comp_loc, 'd_comp_glob': d_comp_glob, 'filters': filters,
        'ksig_loc': ksig_loc,'ksig_glob':ksig_glob}
    rcafit_kw = {'psf_size': psf_size, 'psf_size_type':psf_size_type, 'n_eigenvects': n_eigenvects,
    'nb_iter':n_iter_rca, 'nb_iter_glob':nb_iter_glob, 'nb_iter_loc':nb_iter_loc,
    'nb_subiter_S_loc':nb_subiter_S_loc, 'nb_subiter_A_loc':nb_subiter_A_loc,
    'nb_subiter_S_glob':nb_subiter_S_glob, 'nb_subiter_A_glob':nb_subiter_A_glob,
    'loc_model':loc_model} # [TL] modif

    mccd_rca_fit(starcat, rcainst_kw, rcafit_kw, run_dirs['output'],
                    file_number_string, sex_thresh)

def validate(test_star_path, rca_path, run_dirs, file_number_string,config, w_log):
    sex_thresh = config.getfloat('MCCD', 'SEXMASK_THRESH')
    apply_degradation = config.getboolean('MCCD', 'APPLY_DEGRADATION')
    try:
        mccd_debug = config.getboolean('MCCD', 'MCCD_DEBUG')
    except Exception:
        mccd_debug = False

    star_dict = mccd_validation(rca_path, test_star_path, apply_degradation, mccd_debug,sex_thresh)

    # save in catalog format
    filename = run_dirs['output'] + '/validation_psf'+file_number_string+'.fits'
    output = sc.FITSCatalog(filename, open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                            SEx_catalog=True)
    output.save_as_fits(star_dict, sex_cat_path=test_star_path)



@module_runner(input_module=['mccd_preprocessing_runner'], version='1.0',
               file_pattern=['train_star_selection','test_star_selection'],
               file_ext=['.fits','.fits'], numbering_scheme = '-0000000',
               depends=['numpy', 'mccd_rca', 'galsim'], run_method='parallel')
def mccd_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):

    mode = config.get('MCCD', 'MODE')


    if mode == 'FIT':
        fit(input_file_list, run_dirs, file_number_string,config, w_log)

    elif mode == 'FIT_VALIDATION':
        fit(input_file_list, run_dirs, file_number_string,config, w_log)
        # Validation stars are in the second position of the list
        test_star_path = input_file_list[1]
        # Fitted model is found in the output directory
        rca_path = run_dirs['output'] + '/fitted_model' + file_number_string + '.npy'
        validate(test_star_path, rca_path, run_dirs, file_number_string, config, w_log)

    elif mode == 'VALIDATION':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In VALIDATION mode, both RCA outputs (as .npy) and star catalogs (as .fits) are expected.')
        apply_degradation = config.getboolean('MCCD', 'APPLY_DEGRADATION')
        mccd_debug = config.getboolean('MCCD', 'MCCD_DEBUG')
        # Validation stars are in the second position of the list
        test_star_path = input_file_list[1]
        # Fitted model is the first position in VALIDATION mode
        rca_path = input_file_list[0]
        validate(test_star_path, rca_path, run_dirs, file_number_string, config, w_log)

    elif mode == 'RESPONSE':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In RESPONSE mode, both RCA outputs (as .npy) and catalogs (as .fits) are expected.')

        apply_degradation = config.getboolean('MCCD', 'APPLY_DEGRADATION')
        mccd_debug = config.getboolean('MCCD', 'MCCD_DEBUG')
        x_grid = config.getint('MCCD', 'X_GRID')
        y_grid = config.getint('MCCD', 'Y_GRID')

        rca_path = input_file_list[0]
        test_path = input_file_list[1]

        grid_xy = np.array([x_grid, y_grid])

        star_dict = mccd_response(rca_path, grid_xy, apply_degradation, mccd_debug)

        # save in catalog format
        filename = run_dirs['output'] + '/response_psf'+file_number_string+'.fits'
        output = sc.FITSCatalog(filename, open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                SEx_catalog=True)
        output.save_as_fits(star_dict, sex_cat_path=test_path)




    else:
        raise ValueError('MODE should be in ["FIT", "FIT_VALIDATION", "VALIDATION", "RESPONSE"].')

    return None, None
