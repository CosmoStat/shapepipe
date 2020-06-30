# -*- coding: utf-8 -*-

"""PSFEx MEANSHAPES PLOT RUNNER

This module is used to generate a series of plots from the merged validation
catalogs.

:Author: Tobias Liaudat from Axel Guinot's code

"""
import sys
import numpy as np
from astropy.io import fits
from shapepipe.modules.module_decorator import module_runner
from shapepipe.pipeline import file_io as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import mccd_rca.mccd_utils as mccd_utils


# MegaCam -> plt.subplot correspondance, as given by:
'''        'COMMENT Unique detector IDs for MegaCam (North on top, East to the left)',
           'COMMENT    --------------------------',
           'COMMENT    ba ba ba ba ba ba ba ba ba',
           'COMMENT    00 01 02 03 04 05 06 07 08',
           'COMMENT --------------------------------',
           'COMMENT ba ba ba ba ba ba ba ba ba ba ba',
           'COMMENT 36 09 10 11 12 13 14 15 16 17 37',
           'COMMENT --------------------------------',
           'COMMENT 38 18 19 20 21 22 23 24 25 26 39',
           'COMMENT ab ab ab ab ab ab ab ab ab ab ab',
           'COMMENT --------------------------------',
           'COMMENT    27 28 29 30 31 32 33 34 35',
           'COMMENT    ab ab ab ab ab ab ab ab ab',
           'COMMENT    __________________________'
'''
def MegaCamPos(j):
    if j < 9:
        # first row - shift by one
        return j+1
    elif j < 18:
        # second row, non-ears
        return j+3
    elif j < 27:
        # third row non-ears
        return j+5
    elif j < 36:
        # fourth row
        return j+7
    else:
        MegaCamCoords = {36:11, 37:21, 38:22, 39:32}
        return MegaCamCoords[j]

def MegaCamFlip(xbins, ybins, ccd_nb, nb_pixel):
    if ccd_nb < 18 or ccd_nb in [36,37]:
        # swap x axis so origin is on top-right
        xbins = nb_pixel[0] - xbins + 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1
    return xbins, ybins


def MeanShapesPlot(ccd_maps, filename, title='', colorbar_ampl=1., wind=None, cmap='bwr'):
    # colorbar amplitude
    if wind is None:
        vmax = max(np.nanmax(ccd_maps), np.abs(np.nanmin(ccd_maps))) * colorbar_ampl
        vmin = -vmax * colorbar_ampl
    else:
        vmin, vmax = wind[0]*colorbar_ampl, wind[1]*colorbar_ampl

    # create full plot
    fig, axes = plt.subplots(nrows=4,ncols=11, figsize=(18,12), dpi = 400)
    # remove corner axes (above and below ears)
    for j in [0, 10, -1, -11]:
        axes.flat[j].axis('off')
    for ccd_nb,ccd_map in enumerate(ccd_maps):
        ax = axes.flat[MegaCamPos(ccd_nb)]
        im = ax.imshow(ccd_map.T,cmap=cmap,interpolation='Nearest',vmin=vmin,vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('rmse=%.3e'%(np.sqrt(np.nanmean(ccd_map**2))), size=8)
    plt.suptitle(title, size=20) #TODO: fix title
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('{}.png'.format(filename))
    plt.close()





@module_runner(input_module=['psfex_merge_starcat_runner'], version='1.0',
               file_pattern=['full_starcat'],
               file_ext=['.fits'],numbering_scheme = '-0000000',
               depends=['numpy', 'mccd_rca','astropy', 'matplotlib'],
               run_method='serial')
def psfex_meanshapes_plots_runner(input_file_list, run_dirs, file_number_string,
                       config, w_log):
    # Define the backend for matplotlib
    mpl.use('agg')

    # Get user defined parameters
    x_nb_bins = config.getint('PSFEX_MEANSHAPES_PLOTS', 'X_GRID')
    y_nb_bins = config.getint('PSFEX_MEANSHAPES_PLOTS', 'Y_GRID')
    remove_outliers = config.getboolean('PSFEX_MEANSHAPES_PLOTS', 'REMOVE_OUTLIERS')

    # Define some other paramters
    nb_pixel = x_nb_bins, y_nb_bins
    starcat_path = input_file_list[0][0]
    output_path = run_dirs['output'] + '/'
    auto_colorbar = False
    colorbar_ampl = 1.

    # MegaCam: each CCD is 2048x4612
    grid = np.linspace(0, 2048, nb_pixel[0]+1), np.linspace(0, 4612, nb_pixel[1]+1)

    # READ FULL STARCAT
    starcat = fits.open(starcat_path)[2].data

    # Flag mask
    star_flags = starcat['FLAG_STAR_HSM']
    psf_flags = starcat['FLAG_PSF_HSM']
    flagmask = np.abs(star_flags-1) * np.abs(psf_flags-1)

    # convert sigma to R^2's
    all_star_shapes = np.array([starcat['E1_STAR_HSM'],starcat['E2_STAR_HSM'],
                            2.*starcat['SIGMA_STAR_HSM']**2])
    all_psf_shapes = np.array([starcat['E1_PSF_HSM'],starcat['E2_PSF_HSM'],
                            2.*starcat['SIGMA_PSF_HSM']**2])

    all_CCDs = starcat['CCD_NB']
    all_X = starcat['X']
    all_Y = starcat['Y']

    if remove_outliers == True:
        shape_std_max = 5.
        # Outlier rejection based on the size
        R2_thresh = shape_std_max*np.std(all_psf_shapes[2,:])+np.mean(all_psf_shapes[2,:])
        bad_stars = (abs(all_psf_shapes[2,:])>R2_thresh)
        bad_stars_idx = np.nonzero(bad_stars)[0]
        print('Nb of outlier stars: %d'%(np.sum(bad_stars)))
        # Remove outlier PSFs
        all_star_shapes = all_star_shapes[:,~bad_stars]
        all_psf_shapes = all_psf_shapes[:,~bad_stars]
        all_CCDs = all_CCDs[~bad_stars]
        all_X = all_X[~bad_stars]
        all_Y = all_Y[~bad_stars]
        flagmask = flagmask[~bad_stars]

    # Transform strings into int
    all_CCDs = np.array([int(ccd) for ccd in all_CCDs])

    w_log.info('TOTAL e1 residual RMSE: %.6f\n'%(np.sqrt(np.mean((all_star_shapes[0,:] - all_psf_shapes[0,:])**2))))
    w_log.info('TOTAL e2 residual RMSE: %.6f\n'%(np.sqrt(np.mean((all_star_shapes[1,:] - all_psf_shapes[1,:])**2))))
    w_log.info('TOTAL R2 residual RMSE: %.6f\n'%(np.sqrt(np.mean((all_star_shapes[2,:] - all_psf_shapes[2,:])**2))))

    ccd_maps = np.ones((40, 2, 4)+nb_pixel) * np.nan # CCDs x star/model x (e1,e2,R2,nstars) x xpos x ypos
    for ccd_nb,ccd_map in enumerate(ccd_maps):
        # handle different scatalog versions
        try:
            ccd_mask = ((all_CCDs==ccd_nb) * flagmask).astype(bool)
        except TypeError:
            ccd_mask = ((all_CCDs==str(ccd_nb)) * flagmask).astype(bool)

        star_shapes = all_star_shapes[:,ccd_mask]
        psf_shapes = all_psf_shapes[:,ccd_mask]
        xs, ys = all_X[ccd_mask], all_Y[ccd_mask]

        # if (ccd_nb >= 18 and ccd_nb not in [36, 37]) or ccd_nb in [38, 39]:
        #     xs = 2048 - xs + 1
        # if ccd_nb >= 27 and ccd_nb <=35:
        #     ys = 4612 - ys + 1
        #
        # if ccd_nb >= 18 and ccd_nb <=26:
        #     ys = 4612 - ys + 1
        # elif ccd_nb in [38,39]:
        #     ys = 4612 - ys + 1


        xbins = np.digitize(xs, grid[0])
        ybins = np.digitize(ys, grid[1])

        # swap axes to match CCD orientation and origin convention
        xbins, ybins = MegaCamFlip(xbins, ybins, ccd_nb, nb_pixel)


        for xb in range(nb_pixel[0]):
            for yb in range(nb_pixel[1]):
                bin_star_shapes = star_shapes[:,(xbins==xb+1) * (ybins==yb+1)]
                bin_psf_shapes = psf_shapes[:,(xbins==xb+1) * (ybins==yb+1)]
                ccd_map[0,:3,xb,yb] = np.mean(bin_star_shapes,axis=1)
                ccd_map[1,:3,xb,yb] = np.mean(bin_psf_shapes,axis=1)
                ccd_map[:,3,xb,yb] = bin_star_shapes.shape[1]

    # e_1
    vmax = max(np.nanmax(ccd_maps[:,:,0]), np.abs(np.nanmin(ccd_maps[:,:,0])))
    vmax = 0.125
    vmin = -vmax
    wind = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,0], output_path+'e1s', 'e_1 (stars), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,0,0]),np.nanmax(abs(ccd_maps[:,0,0]))), wind=wind)
    MeanShapesPlot(ccd_maps[:,1,0], output_path+'e1m', 'e_1 (model), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,1,0]),np.nanmax(abs(ccd_maps[:,1,0]))), wind=wind)
    if auto_colorbar:
        wind=None
    e1_res = ccd_maps[:,0,0]-ccd_maps[:,1,0]
    e1_res = e1_res[~np.isnan(e1_res)]
    rmse_e1 = np.sqrt(np.mean((e1_res)**2))
    w_log.info('Bins: e1 residual RMSE: %.6f\n'%(rmse_e1))
    vmax = np.nanmax(abs(ccd_maps[:,0,0]-ccd_maps[:,1,0]))
    vmin = -vmax
    wind = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,0]-ccd_maps[:,1,0], output_path+'e1res', 'e_1 res, rmse=%.5e\nvmax=%.4e , std=%.5e'
    %(rmse_e1,vmax,np.nanstd(ccd_maps[:,0,0]-ccd_maps[:,1,0])), wind=wind, colorbar_ampl=colorbar_ampl)

    # e_2
    vmax = max(np.nanmax(ccd_maps[:,:,1]), np.abs(np.nanmin(ccd_maps[:,:,1])))
    vmax = 0.152
    vmin = -vmax
    wind = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,1], output_path+'e2s', 'e_2 (stars), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,0,1]),np.nanmax(abs(ccd_maps[:,0,1]))), wind=wind)
    MeanShapesPlot(ccd_maps[:,1,1], output_path+'e2m', 'e_2 (model), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,1,1]),np.nanmax(abs(ccd_maps[:,1,1]))), wind=wind)
    if auto_colorbar:
        wind=None
        colorbar_ampl = 1.
    e2_res = ccd_maps[:,0,1]-ccd_maps[:,1,1]
    e2_res = e2_res[~np.isnan(e2_res)]
    rmse_e2 = np.sqrt(np.mean((e2_res)**2))
    w_log.info('Bins: e2 residual RMSE: %.6f\n'%(rmse_e2))
    vmax = np.nanmax(abs(ccd_maps[:,0,1]-ccd_maps[:,1,1]))
    vmin = -vmax
    wind = [vmin, vmax]
    MeanShapesPlot(ccd_maps[:,0,1]-ccd_maps[:,1,1], output_path+'e2res',
        'e_2 res, rmse=%.5e\nvmax=%.4e , std=%.5e'%(rmse_e2,vmax,np.nanstd(ccd_maps[:,0,1]-ccd_maps[:,1,1])), wind=wind, colorbar_ampl=colorbar_ampl)

    # R^2
    vmax = np.nanmax(ccd_maps[:,:,2])
    wind = [0,vmax]
    colorbar_ampl = 1
    MeanShapesPlot(ccd_maps[:,0,2], output_path+'R2s', 'R_2 (stars), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,0,2]),np.nanmax(abs(ccd_maps[:,0,2]))), wind=wind, cmap='Reds')
    MeanShapesPlot(ccd_maps[:,1,2], output_path+'R2m', 'R_2 (model), std=%.5e\nvmax=%.4e'%
        (np.nanstd(ccd_maps[:,1,2]),np.nanmax(abs(ccd_maps[:,1,2]))), wind=wind, cmap='Reds')
    if auto_colorbar:
        wind=[0,np.nanmax(np.abs((ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]))]
        colorbar_ampl = 1.
    R2_res = (ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]
    R2_res = R2_res[~np.isnan(R2_res)]
    rmse_r2 = np.sqrt(np.mean((R2_res)**2))
    w_log.info('Bins: R2 residual RMSE: %.6f\n'%(rmse_r2))
    vmax = np.nanmax(abs((ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]))
    wind = [0,vmax]
    if remove_outliers == True:
        plot_title = 'Outliers removed\n∆(R_2)/R_2 res, rmse=%.5e\nvmax=%.4e , std=%.5e'%(rmse_r2,vmax,np.nanstd((ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]))
    else:
        plot_title = '∆(R_2)/R_2 res, rmse=%.5e\nvmax=%.4e , std=%.5e'%(rmse_r2,vmax,np.nanstd((ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]))

    MeanShapesPlot(np.abs((ccd_maps[:,0,2]-ccd_maps[:,1,2])/ccd_maps[:,0,2]), output_path+'R2res',plot_title,
        wind=wind, colorbar_ampl=colorbar_ampl, cmap='Reds')

    # nstars
    wind=(0,np.max(ccd_maps[:,0,3]))
    MeanShapesPlot(ccd_maps[:,0,3], output_path+'nstar', 'Number of stars\nTotal=%d'
        %(np.nansum(ccd_maps[:,0,3])), wind=wind, cmap='magma')

    # Histograms
    hist_bins = 50
    plt.figure(figsize=(12,6), dpi = 300)
    plt.hist(all_star_shapes[0,:], bins=hist_bins, range=[-0.2, 0.2], label='stars', alpha=0.5)
    plt.hist(all_psf_shapes[0,:], bins=hist_bins, range=[-0.2, 0.2], label='PSFs', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('e1',fontsize=24)
    plt.savefig(output_path+'e1_hist.png')
    plt.close()

    plt.figure(figsize=(12,6), dpi = 300)
    data_hist = all_star_shapes[0,:] - all_psf_shapes[0,:]
    wind = [np.min(data_hist), np.max(data_hist)]
    plt.hist(data_hist, bins=hist_bins, range=wind, label='err(star - psf)', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('e1 err',fontsize=24)
    plt.savefig(output_path+'err_e1_hist.png')
    plt.close()

    plt.figure(figsize=(12,6), dpi = 300)
    plt.hist(all_star_shapes[1,:], bins=hist_bins,range=[-0.2, 0.2], label='stars', alpha=0.5)
    plt.hist(all_psf_shapes[1,:], bins=hist_bins, range=[-0.2, 0.2], label='PSFs', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('e2',fontsize=24)
    plt.savefig(output_path+'e2_hist.png')
    plt.close()

    plt.figure(figsize=(12,6), dpi = 300)
    data_hist = all_star_shapes[1,:] - all_psf_shapes[1,:]
    wind = [np.min(data_hist), np.max(data_hist)]
    plt.hist(data_hist, bins=hist_bins,range=wind, label='err(star - psf)', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('e2 err',fontsize=24)
    plt.savefig(output_path+'err_e2_hist.png')
    plt.close()

    plt.figure(figsize=(12,6), dpi = 300)
    mean_R2 = np.mean(all_star_shapes[2,:])
    wind = [mean_R2 - 4, mean_R2 + 4]
    plt.hist(all_star_shapes[2,:], bins=hist_bins, range=wind, label='stars', alpha=0.5)
    plt.hist(all_psf_shapes[2,:], bins=hist_bins, range=wind, label='PSFs', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('R2',fontsize=24)
    plt.savefig(output_path+'R2_hist.png')
    plt.close()

    plt.figure(figsize=(12,6), dpi = 300)
    data_hist = (all_star_shapes[2,:] - all_psf_shapes[2,:])/all_star_shapes[2,:]
    wind = [np.min(data_hist), np.max(data_hist)]
    # wind = [-0.5,0.5]
    plt.hist(data_hist, bins=hist_bins, range=wind, label='err(star - psf)/star', alpha=0.5)
    plt.legend(loc='best',fontsize=16)
    plt.title('R2 err',fontsize=24)
    plt.savefig(output_path+'err_R2_hist.png')
    plt.close()

    return None,None
