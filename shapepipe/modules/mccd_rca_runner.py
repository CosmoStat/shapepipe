# -*- coding: utf-8 -*-

"""RCA RUNNER

This file is the pipeline runner for the RCA package.

:Author: Morgan Schmitz

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
from astropy.io import fits
import numpy as np
import mccd_rca.my_mccd_rca as rca
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
    masks = starcat[2].data['MASK_LIST']
    ccds = starcat[2].data['CCD_ID_LIST']

    # SNR treatment
    try:
        SNRs = starcat[2].data['SNR_WIN_LIST']
        SNR_weights =  SNRs/np.median(SNRs)             # Strategy N5
        # SNR_weights[SNRs<50] = SNR_weights[SNRs<50]/10. #
        SNR_weights[SNR_weights>2.] = 2.
        SNR_weights[SNR_weights<0.01] = 0.01
        SNR_weight_list = [SNR_weights[my_ccd_list==ccd] for ccd in ccd_unique_list]
    except:
        SNR_weight_list = None
        print('No SNR weights are being used.')

    pos_list = [positions[my_ccd_list==ccd] for ccd in ccd_unique_list]
    star_list = [utils.rca_format(stars[my_ccd_list==ccd]) for ccd in ccd_unique_list]
    mask_list = [utils.rca_format(masks[my_ccd_list==ccd]) for ccd in ccd_unique_list]

    ccd_list = [ccds[my_ccd_list==ccd].astype(int) for ccd in ccd_unique_list]
    ccd_list= [np.unique(_list)[0].astype(int) for _list in ccd_list]

    rca_instance = rca.mccd_RCA(**rcainst_kw, verbose=True)
    S, A_loc, A_glob, alpha, pi = rca_instance.fit(star_list, pos_list, ccd_list,
                                    mask_list, SNR_weight_list, **rcafit_kw)

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
    try:
        star_dict['NUMBER'] = np.copy(star_cat[2].data['NUMBER']) # [JB] to comment when in validation mode
        star_dict['X'] = np.copy(star_cat[2].data[pos_params[0]])
        star_dict['Y'] = np.copy(star_cat[2].data[pos_params[1]])
        star_dict['RA'] = np.copy(star_cat[2].data['XWIN_WORLD']) # [JB] to comment when in validation mode
        star_dict['DEC'] = np.copy(star_cat[2].data['YWIN_WORLD']) # [JB] to comment when in validation mode
        star_dict['MAG'] = np.copy(star_cat[2].data['MAG_AUTO']) # [JB] to comment when in validation mode
        star_dict['SNR'] = np.copy(star_cat[2].data['SNR_WIN']) # [JB] to comment when in validation mode
    except:
        #star_dict['NUMBER'] = np.copy(star_cat[2].data['NUMBER']) # [JB] to comment when in validation mode
        star_dict['X'] = np.copy(star_cat[2].data[pos_params[0]])
        star_dict['Y'] = np.copy(star_cat[2].data[pos_params[1]])
        #star_dict['RA'] = np.copy(star_cat[2].data['XWIN_WORLD']) # [JB] to comment when in validation mode
        #star_dict['DEC'] = np.copy(star_cat[2].data['YWIN_WORLD']) # [JB] to comment when in validation mode
        #star_dict['MAG'] = np.copy(star_cat[2].data['MAG_AUTO']) # [JB] to comment when in validation mode
        #star_dict['SNR'] = np.copy(star_cat[2].data['SNR_WIN']) # [JB] to comment when in validation mode

    # raw_path_PSFs = '/Users/tliaudat/Documents/PhD/codes/venv_p3/all-W3-tests/raw-data/test-RCA_hybrid_NOPSFEx/PSFs/'
    # raw_path_stars = '/Users/tliaudat/Documents/PhD/codes/venv_p3/all-W3-tests/raw-data/test-RCA_hybrid_NOPSFEx/stars/'
    # raw_path_badpix = '/Users/tliaudat/Documents/PhD/codes/venv_p3/all-W3-tests/raw-data/test-RCA_hybrid_NOPSFEx/badpixs/'
    #np.save(raw_path_PSFs + file_number_string + 'PSFs.npy',PSFs)
    #np.save(raw_path_stars + file_number_string +'stars.npy',stars)

    # compute star shapes with HSM
    badpix_mask = np.rint(np.abs(mask-1)) # hsm thinks 0 means good

    #np.save(raw_path_badpix + file_number_string + 'badpix.npy',badpix_mask)

    # Pixel MSE saving [TL]
    # raw_path_pixel_MSE = '/Users/tliaudat/Documents/PhD/codes/venv_p3/tests/MSE_pixel/test-27/'
    raw_path_pixel_MSE = '/Users/tliaudat/Documents/PhD/codes/venv_p3/sandbox_RCAv3/output/val/test-32/'
    doc_name = 'results.txt'
    try:
        f = open(raw_path_pixel_MSE + doc_name)
    except IOError:
        # If it does not exist we write the first line
        f = open(raw_path_pixel_MSE + doc_name,'a')
        f.write('catalogId\tMSE\tnStars\tDx\tDy\n')
        f.close()

    f = open(raw_path_pixel_MSE + doc_name,'a')
    myMSE = np.sum(((PSFs-stars)**2)/(PSFs.shape[0]*PSFs.shape[1]*PSFs.shape[2]))
    f.write('%s\t%.18f\t%d\t%d\t%d\n'%(file_number_string,myMSE,PSFs.shape[0],PSFs.shape[1],PSFs.shape[2]))
    f.close()

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

def mccd_validation(rca_path, test_path, apply_degradation = True, mccd_debug = False):
    # Principal validation for mccd catalogs

    # Import
    rca_instance = rca.mccd_quickload(rca_path)
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
        badpix_mask = np.abs(val_mask_list[it]-1) # hsm thinks 0 means good
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
    rca_instance = rca.mccd_quickload(rca_path)

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




@module_runner(input_module=['setools_runner'], version='1.0',
               file_pattern=['star_selection'],
               file_ext=['.fits'],
               depends=['numpy', 'mccd_rca', 'galsim'])
def mccd_rca_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):
    mode = config.get('RCA_RUNNER', 'MODE')
    sex_thresh = config.getfloat('RCA_RUNNER', 'SEXMASK_THRESH')
    # pos_params = config.getlist('RCA_RUNNER', 'POSITION_PARAMS') # [TL] Not used -> Fixed in global(XWIN_IMAGE)
    star_thresh = config.getint('RCA_RUNNER', 'STAR_THRESH')

    if mode in ['FIT', 'FIT_TRANSFORM']:
        n_comp_loc = config.getint('RCA_RUNNER', 'N_COMP_LOC')
        d_comp_glob = config.getint('RCA_RUNNER', 'D_COMP_GLOB')
        psf_size = config.getfloat('RCA_RUNNER', 'PSF_SIZE')
        psf_size_type = config.get('RCA_RUNNER', 'PSF_SIZE_TYPE')
        n_eigenvects = config.getint('RCA_RUNNER', 'N_EIGENVECTS')
        if n_eigenvects == 0:
            n_eigenvects = None
        ksig_loc = config.getfloat('RCA_RUNNER', 'KSIG_LOC')
        ksig_glob = config.getfloat('RCA_RUNNER', 'KSIG_GLOB')
        n_iter_rca = config.getint('RCA_RUNNER', 'N_ITER_RCA') # [TL] modif
        nb_iter_glob = config.getint('RCA_RUNNER', 'N_ITER_GLOB') # [TL] modif
        nb_iter_loc = config.getint('RCA_RUNNER', 'N_ITER_LOC') # [TL] modif
        nb_subiter_S_loc = config.getint('RCA_RUNNER', 'NB_SUBITER_S_LOC') # [TL] modif
        nb_subiter_A_loc = config.getint('RCA_RUNNER', 'NB_SUBITER_A_LOC') # [TL] modif
        nb_subiter_S_glob = config.getint('RCA_RUNNER', 'NB_SUBITER_S_GLOB') # [TL] modif
        nb_subiter_A_glob = config.getint('RCA_RUNNER', 'NB_SUBITER_A_GLOB') # [TL] modif
        loc_model = config.get('RCA_RUNNER', 'LOC_MODEL')
        filt_path = config.get('RCA_RUNNER', 'FILTER_PATH')
        filters = None if (filt_path == 'None') else np.load(filt_path)


    if mode == 'FIT':
        starcat_path = input_file_list[0]
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

    elif mode == 'RESPONSE':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In RESPONSE mode, both RCA outputs (as .npy) and catalogs (as .fits) are expected.')

        apply_degradation = config.getboolean('RCA_RUNNER', 'APPLY_DEGRADATION')
        mccd_debug = config.getboolean('RCA_RUNNER', 'MCCD_DEBUG')
        x_grid = config.getint('RCA_RUNNER', 'X_GRID')
        y_grid = config.getint('RCA_RUNNER', 'Y_GRID')

        rca_path = input_file_list[0]
        test_path = input_file_list[1]

        grid_xy = np.array([x_grid, y_grid])

        star_dict = mccd_response(rca_path, grid_xy, apply_degradation, mccd_debug)

        # save in catalog format
        filename = run_dirs['output'] + '/response_psf'+file_number_string+'.fits'
        output = sc.FITSCatalog(filename, open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                SEx_catalog=True)
        output.save_as_fits(star_dict, sex_cat_path=test_path)

    elif mode == 'VALIDATION':
        if len(input_file_list) < 2 or '.npy' not in input_file_list[0]:
            raise ValueError('In VALIDATION mode, both RCA outputs (as .npy) and star catalogs (as .fits) are expected.')
        apply_degradation = config.getboolean('RCA_RUNNER', 'APPLY_DEGRADATION')
        mccd_debug = config.getboolean('RCA_RUNNER', 'MCCD_DEBUG')
        rca_path = input_file_list[0]
        test_path = input_file_list[1]

        star_dict = mccd_validation(rca_path, test_path, apply_degradation, mccd_debug)

        # save in catalog format
        filename = run_dirs['output'] + '/validation_psf'+file_number_string+'.fits'
        output = sc.FITSCatalog(filename, open_mode=sc.BaseCatalog.OpenMode.ReadWrite,
                                SEx_catalog=True)
        output.save_as_fits(star_dict, sex_cat_path=test_path)


    else:
        raise ValueError('MODE should be in ["FIT", "TRANSFORM", "FIT_TRANSFORM", "VALIDATION"].')

    return None, None
