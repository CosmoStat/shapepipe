# -*- coding: utf-8 -*-

"""MCCD PLOTS UTILITIES.

This module is used to generate a series of plots from the merged validation
catalogs.
It plots the Meanshape plots for the merged validation catalog.
It can also plot the rho statistics provided that the required packages are
installed.

:Author: Tobias Liaudat base on Axel Guinot's code

"""

import numpy as np
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
import mccd.mccd_utils as mccd_utils
import time
import stile
import stile.stile_utils
from stile.sys_tests import BaseCorrelationFunctionSysTest
import treecorr

# Define the backend for matplotlib
mpl.use('agg')

# MegaCam -> plt.subplot correspondance, as given by:
'''        'COMMENT Unique detector IDs for MegaCam',
           'COMMENT (North on top, East to the left)',
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
        return j + 1
    elif j < 18:
        # second row, non-ears
        return j + 3
    elif j < 27:
        # third row non-ears
        return j + 5
    elif j < 36:
        # fourth row
        return j + 7
    else:
        MegaCamCoords = {36: 11, 37: 21, 38: 22, 39: 32}
        return MegaCamCoords[j]


def MegaCamFlip(xbins, ybins, ccd_nb, nb_pixel):
    if ccd_nb < 18 or ccd_nb in [36, 37]:
        # swap x axis so origin is on top-right
        xbins = nb_pixel[0] - xbins + 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1
    return xbins, ybins


def MegaCamFlip_2(xbins, ybins, ccd_nb, nb_pixel):
    if ccd_nb < 18 or ccd_nb in [36, 37]:
        # swap x axis so origin is on top-right
        # xbins = nb_pixel[0] - xbins + 1
        a = 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1
    return xbins, ybins


def MeanShapesPlot(ccd_maps, filename, title='', colorbar_ampl=1., wind=None,
                   cmap='bwr'):
    r"""Plot meanshapes from ccd maps."""
    # colorbar amplitude
    if wind is None:
        vmax = max(
            np.nanmax(ccd_maps), np.abs(np.nanmin(ccd_maps))) * colorbar_ampl
        vmin = -vmax * colorbar_ampl
    else:
        vmin, vmax = wind[0] * colorbar_ampl, wind[1] * colorbar_ampl

    # create full plot
    fig, axes = plt.subplots(nrows=4, ncols=11, figsize=(18, 12), dpi=400)
    # remove corner axes (above and below ears)
    for j in [0, 10, -1, -11]:
        axes.flat[j].axis('off')
    for ccd_nb, ccd_map in enumerate(ccd_maps):
        ax = axes.flat[MegaCamPos(ccd_nb)]
        im = ax.imshow(ccd_map.T, cmap=cmap, interpolation='Nearest',
                       vmin=vmin, vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title('rmse=%.3e' % (np.sqrt(np.nanmean(ccd_map ** 2))), size=8)
    plt.suptitle(title, size=20)  # TODO: fix title
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig('{}.png'.format(filename))
    plt.close()


def plot_meanshapes(starcat_path, output_path, nb_pixel, w_log,
                    remove_outliers=False, plot_meanshapes=True,
                    plot_histograms=True):
    r"""Plot meanshapes and histograms."""
    # READ FULL STARCAT
    starcat = fits.open(starcat_path)

    auto_colorbar = False
    colorbar_ampl = 1.
    loc2glob = mccd_utils.Loc2Glob()

    # MegaCam: each CCD is 2048x4612
    grid = np.linspace(0, loc2glob.x_npix, nb_pixel[0] + 1), \
        np.linspace(0, loc2glob.y_npix, nb_pixel[1] + 1)

    # Flag mask
    star_flags = starcat[2].data['FLAG_STAR_HSM']
    psf_flags = starcat[2].data['FLAG_PSF_HSM']
    flagmask = np.abs(star_flags - 1) * np.abs(psf_flags - 1)

    # convert sigma to R^2's
    all_star_shapes = np.array(
        [starcat[2].data['E1_STAR_HSM'], starcat[2].data['E2_STAR_HSM'],
         2. * starcat[2].data['SIGMA_STAR_HSM'] ** 2])
    all_psf_shapes = np.array(
        [starcat[2].data['E1_PSF_HSM'], starcat[2].data['E2_PSF_HSM'],
         2. * starcat[2].data['SIGMA_PSF_HSM'] ** 2])
    all_CCDs = starcat[2].data['CCD_NB']
    all_X = starcat[2].data['X']
    all_Y = starcat[2].data['Y']

    if remove_outliers:
        shape_std_max = 5.
        # Outlier rejection based on the size
        R2_thresh = shape_std_max * np.std(all_psf_shapes[2, :]) + np.mean(
            all_psf_shapes[2, :])
        bad_stars = (abs(all_psf_shapes[2, :]) > R2_thresh)
        bad_stars_idx = np.nonzero(bad_stars)[0]
        w_log.info('Nb of outlier stars: %d' % (np.sum(bad_stars)))
        # Remove outlier PSFs
        all_star_shapes = all_star_shapes[:, ~bad_stars]
        all_psf_shapes = all_psf_shapes[:, ~bad_stars]
        all_CCDs = all_CCDs[~bad_stars]
        all_X = all_X[~bad_stars]
        all_Y = all_Y[~bad_stars]
        flagmask = flagmask[~bad_stars]

    w_log.info('TOTAL e1 residual RMSE: %.6e\n' % (
        np.sqrt(np.mean((all_star_shapes[0, :] - all_psf_shapes[0, :]) ** 2))))
    w_log.info('TOTAL e2 residual RMSE: %.6e\n' % (
        np.sqrt(np.mean((all_star_shapes[1, :] - all_psf_shapes[1, :]) ** 2))))
    w_log.info('TOTAL R2 residual RMSE: %.6e\n' % (
        np.sqrt(np.mean((all_star_shapes[2, :] - all_psf_shapes[2, :]) ** 2))))

    # CCDs x star/model x (e1,e2,R2,nstars) x xpos x ypos
    ccd_maps = np.ones((40, 2, 4) + nb_pixel) * np.nan

    for ccd_nb, ccd_map in enumerate(ccd_maps):
        # handle different scatalog versions
        try:
            ccd_mask = ((all_CCDs.astype(int) == ccd_nb) * flagmask).astype(
                bool)
        except TypeError:
            ccd_mask = ((all_CCDs == str(ccd_nb)) * flagmask).astype(bool)

        star_shapes = all_star_shapes[:, ccd_mask]
        psf_shapes = all_psf_shapes[:, ccd_mask]

        # Claculate shift to go from global coordinates to local coordinates
        x_shift, y_shift = loc2glob.shift_coord(ccd_nb)
        xs_loc, ys_loc = all_X[ccd_mask] - x_shift, all_Y[ccd_mask] - y_shift

        # swap axes to match CCD orientation and origin convention
        ys_loc = loc2glob.y_npix - ys_loc + 1

        # digitalize into bins
        xbins = np.digitize(xs_loc, grid[0])
        ybins = np.digitize(ys_loc, grid[1])

        for xb in range(nb_pixel[0]):
            for yb in range(nb_pixel[1]):
                bin_star_shapes = star_shapes[
                    :, (xbins == xb + 1) * (ybins == yb + 1)]
                bin_psf_shapes = psf_shapes[
                    :, (xbins == xb + 1) * (ybins == yb + 1)]
                ccd_map[0, :3, xb, yb] = np.mean(bin_star_shapes, axis=1)
                ccd_map[1, :3, xb, yb] = np.mean(bin_psf_shapes, axis=1)
                ccd_map[:, 3, xb, yb] = bin_star_shapes.shape[1]

    if plot_meanshapes:
        # e_1
        vmax = max(np.nanmax(ccd_maps[:, :, 0]),
                   np.abs(np.nanmin(ccd_maps[:, :, 0])))
        vmin = -vmax
        wind = [vmin, vmax]
        MeanShapesPlot(
            ccd_maps[:, 0, 0], output_path + 'e1s',
            "e_1 (stars), std=%.5e\nvmax=%.4e" % (
                np.nanstd(ccd_maps[:, 0, 0]),
                np.nanmax(abs(ccd_maps[:, 0, 0]))),
            wind=wind)

        MeanShapesPlot(
            ccd_maps[:, 1, 0], output_path + 'e1m',
            'e_1 (model), std=%.5e\nvmax=%.4e' %
            (np.nanstd(ccd_maps[:, 1, 0]),
             np.nanmax(abs(ccd_maps[:, 1, 0]))),
            wind=wind)

        if auto_colorbar:
            wind = None
        e1_res = ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0]
        e1_res = e1_res[~np.isnan(e1_res)]
        rmse_e1 = np.sqrt(np.mean(e1_res ** 2))
        w_log.info('Bins: e1 residual RMSE: %.6f\n' % rmse_e1)
        vmax = np.nanmax(abs(ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0]))
        vmin = -vmax
        wind = [vmin, vmax]
        MeanShapesPlot(
            ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0],
            output_path + 'e1res',
            'e_1 res, rmse=%.5e\nvmax=%.4e , std=%.5e' % (
                rmse_e1, vmax,
                np.nanstd(ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0])),
            wind=wind,
            colorbar_ampl=colorbar_ampl)

        # e_2
        vmax = max(np.nanmax(ccd_maps[:, :, 1]),
                   np.abs(np.nanmin(ccd_maps[:, :, 1])))
        vmin = -vmax
        wind = [vmin, vmax]
        MeanShapesPlot(ccd_maps[:, 0, 1], output_path + 'e2s',
                       'e_2 (stars), std=%.5e\nvmax=%.4e' %
                       (np.nanstd(ccd_maps[:, 0, 1]),
                        np.nanmax(abs(ccd_maps[:, 0, 1]))), wind=wind)
        MeanShapesPlot(ccd_maps[:, 1, 1], output_path + 'e2m',
                       'e_2 (model), std=%.5e\nvmax=%.4e' %
                       (np.nanstd(ccd_maps[:, 1, 1]),
                        np.nanmax(abs(ccd_maps[:, 1, 1]))), wind=wind)

        if auto_colorbar:
            wind = None
            colorbar_ampl = 1.
        e2_res = ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1]
        e2_res = e2_res[~np.isnan(e2_res)]
        rmse_e2 = np.sqrt(np.mean(e2_res ** 2))
        w_log.info('Bins: e2 residual RMSE: %.6f\n' % rmse_e2)
        vmax = np.nanmax(abs(ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1]))
        vmin = -vmax
        wind = [vmin, vmax]
        MeanShapesPlot(
            ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1], output_path + 'e2res',
            'e_2 res, rmse=%.5e\nvmax=%.4e , std=%.5e' % (
                rmse_e2, vmax,
                np.nanstd(ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1])),
            wind=wind, colorbar_ampl=colorbar_ampl)

        # R^2
        wind = [0, np.nanmax(ccd_maps[:, :, 2])]
        colorbar_ampl = 1
        MeanShapesPlot(
            ccd_maps[:, 0, 2], output_path + 'R2s',
            'R_2 (stars), std=%.5e\nvmax=%.4e' % (
                np.nanstd(ccd_maps[:, 0, 2]),
                np.nanmax(abs(ccd_maps[:, 0, 2]))),
            wind=wind, cmap='Reds')

        MeanShapesPlot(
            ccd_maps[:, 1, 2], output_path + 'R2m',
            'R_2 (model), std=%.5e\nvmax=%.4e' % (
                np.nanstd(ccd_maps[:, 1, 2]),
                np.nanmax(abs(ccd_maps[:, 1, 2]))),
            wind=wind, cmap='Reds')

        if auto_colorbar:
            wind = [0, np.nanmax(np.abs(
                (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]))]
            colorbar_ampl = 1.
        R2_res = (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]
        R2_res = R2_res[~np.isnan(R2_res)]
        rmse_r2 = np.sqrt(np.mean(R2_res ** 2))
        w_log.info('Bins: R2 residual RMSE: %.6f\n' % rmse_r2)
        vmax = np.nanmax(
            abs((ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]))
        wind = [0, vmax]
        if remove_outliers:
            plot_title = "Outliers removed\n∆(R_2)/R_2 res, " \
                         "rmse=%.5e\nvmax=%.4e , std=%.5e" % (
                             rmse_r2, vmax,
                             np.nanstd(
                                 (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) /
                                 ccd_maps[:, 0, 2]))
        else:
            plot_title = "∆(R_2)/R_2 res, rmse=%.5e\nvmax=%.4e , std=%.5e" % (
                rmse_r2, vmax,
                np.nanstd((ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) /
                          ccd_maps[:, 0, 2]))

        MeanShapesPlot(
            np.abs((ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) /
                   ccd_maps[:, 0, 2]),
            output_path + 'R2res',
            plot_title, wind=wind, colorbar_ampl=colorbar_ampl, cmap='Reds')

        # nstars
        wind = (0, np.max(ccd_maps[:, 0, 3]))
        MeanShapesPlot(
            ccd_maps[:, 0, 3], output_path + 'nstar',
            'Number of stars\nTotal=%d' % (np.nansum(ccd_maps[:, 0, 3])),
            wind=wind, cmap='magma')

    # Histograms
    if plot_histograms:
        hist_bins = 50
        plt.figure(figsize=(12, 6), dpi=300)
        plt.hist(all_star_shapes[0, :], bins=hist_bins, range=[-0.2, 0.2],
                 label='stars', alpha=0.5)
        plt.hist(all_psf_shapes[0, :], bins=hist_bins, range=[-0.2, 0.2],
                 label='PSFs', alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('e1', fontsize=24)
        plt.savefig(output_path + 'e1_hist.png')
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = all_star_shapes[0, :] - all_psf_shapes[0, :]
        plt.hist(data_hist, bins=hist_bins,
                 range=[np.min(data_hist), np.max(data_hist)],
                 label='err(star - psf)', alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('e1 err', fontsize=24)
        plt.savefig(output_path + 'err_e1_hist.png')
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        plt.hist(all_star_shapes[1, :], bins=hist_bins, range=[-0.2, 0.2],
                 label='stars', alpha=0.5)
        plt.hist(all_psf_shapes[1, :], bins=hist_bins, range=[-0.2, 0.2],
                 label='PSFs', alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('e2', fontsize=24)
        plt.savefig(output_path + 'e2_hist.png')
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = all_star_shapes[1, :] - all_psf_shapes[1, :]
        plt.hist(data_hist, bins=hist_bins,
                 range=[np.min(data_hist), np.max(data_hist)],
                 label='err(star - psf)', alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('e2 err', fontsize=24)
        plt.savefig(output_path + 'err_e2_hist.png')
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        mean_R2 = np.mean(all_star_shapes[2, :])
        wind = [mean_R2 - 4, mean_R2 + 4]
        plt.hist(all_star_shapes[2, :], bins=hist_bins, range=wind, label='stars',
                 alpha=0.5)
        plt.hist(all_psf_shapes[2, :], bins=hist_bins, range=wind, label='PSFs',
                 alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('R2', fontsize=24)
        plt.savefig(output_path + 'R2_hist.png')
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = (all_star_shapes[2, :] -
                     all_psf_shapes[2, :]) / all_star_shapes[2, :]
        plt.hist(data_hist, bins=hist_bins,
                 range=[np.min(data_hist), np.max(data_hist)],
                 label='err(star - psf)/star', alpha=0.5)
        plt.legend(loc='best', fontsize=16)
        plt.title('R2 err', fontsize=24)
        plt.savefig(output_path + 'err_R2_hist.png')
        plt.close()

    starcat.close()

    return None


# Rho stats functions

def NegDash(x_in, y_in, yerr_in, plot_name='', vertical_lines=True,
            xlabel='', ylabel='', rho_nb='',
            ylim=None, semilogx=False, semilogy=False,
            **kwargs):
    """ This function is for making plots with vertical errorbars,
    where negative values are shown in absolute value as dashed lines.
    The resulting plot can either be saved by specifying a file name as
    `plot_name', or be kept as a pyplot instance (for instance to combine
    several NegDashes).
    """
    x = np.copy(x_in)
    y = np.copy(y_in)
    yerr = np.copy(yerr_in)
    # catch and separate errorbar-specific keywords from Lines2D ones
    safekwargs = dict(kwargs)
    errbkwargs = dict()
    if 'linestyle' in kwargs.keys():
        print("""Warning: linestyle was provided but that would kind of defeat
        the purpose, so I'll just ignore it. Sorry.""")
        del safekwargs['linestyle']
    for errorbar_kword in ['fmt', 'ecolor', 'elinewidth', 'capsize',
                           'barsabove', 'errorevery']:
        if errorbar_kword in kwargs.keys():
            # posfmt = '-'+kwargs['fmt']
            # negfmt = '--'+kwargs['fmt']
            errbkwargs[errorbar_kword] = kwargs[errorbar_kword]
            del safekwargs[errorbar_kword]
    errbkwargs = dict(errbkwargs, **safekwargs)

    # plot up to next change of sign
    current_sign = np.sign(y[0])
    first_change = np.argmax(current_sign * y < 0)
    while first_change:
        if current_sign > 0:
            plt.errorbar(x[:first_change], y[:first_change],
                         yerr=yerr[:first_change],
                         linestyle='-', **errbkwargs)
            if vertical_lines:
                plt.vlines(x[first_change - 1], 0, y[first_change - 1],
                           linestyle='-', **safekwargs)
                plt.vlines(x[first_change], 0, np.abs(y[first_change]),
                           linestyle='--', **safekwargs)
        else:
            plt.errorbar(x[:first_change], np.abs(y[:first_change]),
                         yerr=yerr[:first_change],
                         linestyle='--', **errbkwargs)
            if vertical_lines:
                plt.vlines(x[first_change - 1], 0, np.abs(y[first_change - 1]),
                           linestyle='--',
                           **safekwargs)
                plt.vlines(x[first_change], 0, y[first_change], linestyle='-',
                           **safekwargs)
        x = x[first_change:]
        y = y[first_change:]
        yerr = yerr[first_change:]
        current_sign *= -1
        first_change = np.argmax(current_sign * y < 0)
    # one last time when `first_change'==0 ie no more changes:
    if rho_nb:
        lab = r'$\rho_{}(\theta)$'.format(rho_nb)
    else:
        lab = ''
    if current_sign > 0:
        plt.errorbar(x, y, yerr=yerr, linestyle='-', label=lab, **errbkwargs)
    else:
        plt.errorbar(x, np.abs(y), yerr=yerr, linestyle='--', label=lab,
                     **errbkwargs)
    if semilogx:
        plt.xscale('log')
    if semilogy:
        plt.yscale('log')
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if plot_name:
        plt.savefig(plot_name)
        plt.close()


class new_BaseCorrelationFunctionSysTest(BaseCorrelationFunctionSysTest):
    r"""Base function for the correlation.

    Based on style package class.

    """

    def makeCatalog(self, data, config=None, use_as_k=None,
                    use_chip_coords=False):
        if data is None or isinstance(data, treecorr.Catalog):
            return data
        catalog_kwargs = {}
        fields = data.dtype.names
        if 'ra' in fields and 'dec' in fields:
            if not use_chip_coords:
                catalog_kwargs['ra'] = data['ra']
                catalog_kwargs['dec'] = data['dec']
            elif 'x' in fields and 'y' in fields:
                catalog_kwargs['x'] = data['x']
                catalog_kwargs['y'] = data['y']
            else:
                raise ValueError(
                    "Chip coordinates requested, but 'x' and 'y' fields"
                    "not found in data")
        elif 'x' in fields and 'y' in fields:
            catalog_kwargs['x'] = data['x']
            catalog_kwargs['y'] = data['y']
        else:
            raise ValueError(
                "Data must contain (ra,dec) or (x,y) in order to do"
                "correlation function tests.")
        if 'g1' in fields and 'g2' in fields:
            catalog_kwargs['g1'] = data['g1']
            catalog_kwargs['g2'] = data['g2']
        if 'w' in fields:
            catalog_kwargs['w'] = data['w']
        if use_as_k:
            if use_as_k in fields:
                catalog_kwargs['k'] = data[use_as_k]
        elif 'k' in fields:
            catalog_kwargs['k'] = data['k']
        # Quirk of length-1 formatted arrays: the fields will be floats, not
        # arrays, which would break the Catalog init.
        try:
            len(data)
        except:
            if not hasattr(data, 'len') and isinstance(data, np.ndarray):
                for key in catalog_kwargs:
                    catalog_kwargs[key] = np.array([catalog_kwargs[key]])
        catalog_kwargs['config'] = config
        return treecorr.Catalog(**catalog_kwargs)


class Rho1SysTest(new_BaseCorrelationFunctionSysTest):
    """
    Compute the auto-correlation of residual star shapes
    (star shapes - psf shapes).
    """
    short_name = 'rho1'
    long_name = 'Rho1 statistics (Auto-correlation of star-PSF shapes)'
    objects_list = ['star PSF']
    required_quantities = [('ra', 'dec', 'g1', 'g2', 'psf_g1', 'psf_g2', 'w')]

    def __call__(self, data, data2=None, random=None, random2=None,
                 config=None, **kwargs):
        new_data = data.copy()
        new_data['g1'] = new_data['g1'] - new_data['psf_g1']
        new_data['g2'] = new_data['g2'] - new_data['psf_g2']
        if data2 is not None:
            new_data2 = data2.copy()
            new_data2['g1'] = new_data2['g1'] - new_data2['psf_g1']
            new_data2['g2'] = new_data2['g2'] - new_data2['psf_g2']
        else:
            new_data2 = data2
        if random is not None:
            new_random = random.copy()
            new_random['g1'] = new_random['g1'] - new_random['psf_g1']
            new_random['g2'] = new_random['g2'] - new_random['psf_g2']
        else:
            new_random = random
        if random2 is not None:
            new_random2 = random2.copy()
            new_random2['g1'] = new_random2['g1'] - new_random2['psf_g1']
            new_random2['g2'] = new_random2['g2'] - new_random2['psf_g2']
        else:
            new_random2 = random2
        return self.getCF('gg', new_data, new_data2, new_random, new_random2,
                          config=config, **kwargs)


class DESRho2SysTest(new_BaseCorrelationFunctionSysTest):
    """
    Compute the correlation of PSF shapes with residual star shapes
    (star shapes - psf shapes).
    """
    short_name = 'rho2des'
    long_name = 'Rho2 statistics (as defined in DES shape catalogue papers)'
    objects_list = ['star PSF']
    required_quantities = [('ra', 'dec', 'g1', 'g2', 'psf_g1', 'psf_g2', 'w')]

    def __call__(self, data, data2=None, random=None, random2=None,
                 config=None, **kwargs):
        new_data = np.rec.fromarrays([data['ra'], data['dec'], data['g1'],
                                      data['g2'], data['w']],
                                     names=['ra', 'dec', 'g1', 'g2', 'w'])
        if data2 is None:
            data2 = data
        new_data2 = np.rec.fromarrays(
            [data2['ra'], data2['dec'], data2['g1'] - data2['psf_g1'],
             data2['g2'] - data2['psf_g2'], data2['w']],
            names=['ra', 'dec', 'g1', 'g2', 'w'])
        if random is not None:
            new_random = np.rec.fromarrays(
                [random['ra'], random['dec'], random['g1'],
                 random['g2'], random['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])

        else:
            new_random = random
        if random2 is None:
            random2 = random
        if random2 is not None:
            new_random2 = np.rec.fromarrays([data2['ra'], data2['dec'],
                                             data2['g1'] - data2['psf_g1'],
                                             data2['g2'] - data2['psf_g2'],
                                             data2['w']],
                                            names=['ra', 'dec', 'g1', 'g2',
                                                   'w'])
        else:
            new_random2 = random2
        return self.getCF('gg', new_data, new_data2, new_random, new_random2,
                          config=config, **kwargs)


class DESRho3SysTest(new_BaseCorrelationFunctionSysTest):
    """
    Compute the correlation of star shapes weighted by the residual size.
    """
    short_name = 'rho3'
    long_name = "Rho3 statistics (Auto-correlation of star shapes " \
                "weighted by the residual size)"
    objects_list = ['star PSF']
    required_quantities = [('ra', 'dec', 'sigma',
                            'g1', 'g2', 'psf_sigma', 'w')]

    def __call__(self, data, data2=None, random=None, random2=None,
                 config=None, **kwargs):
        new_data = np.rec.fromarrays(
            [data['ra'], data['dec'],
             data['g1'] * (data['sigma'] - data[
                 'psf_sigma']) / data['sigma'],
             data['g2'] * (data['sigma'] - data[
                 'psf_sigma']) / data['sigma'],
             data['w']],
            names=['ra', 'dec', 'g1', 'g2', 'w'])
        if data2 is not None:
            new_data2 = np.rec.fromarrays(
                [data2['ra'], data2['dec'],
                 data2['g1'] * (data2['sigma'] - data2['psf_sigma']) /
                 data2['sigma'],
                 data2['g2'] * (data2['sigma'] - data2['psf_sigma']) /
                 data2['sigma'],
                 data2['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])

        else:
            new_data2 = data2
        if random is not None:
            new_random = np.rec.fromarrays(
                [random['ra'], random['dec'],
                 random['g1'] * (random['sigma'] - random['psf_sigma']) /
                 random['sigma'],
                 random['g2'] * (random['sigma'] - random['psf_sigma']) /
                 random['sigma'],
                 random['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])
        else:
            new_random = random

        if random2 is not None:
            new_random2 = np.rec.fromarrays(
                [random2['ra'], random2['dec'],
                 random2['g1'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['g2'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])

        else:
            new_random2 = random2
        return self.getCF('gg', new_data, new_data2, new_random, new_random2,
                          config=config, **kwargs)


class DESRho4SysTest(new_BaseCorrelationFunctionSysTest):
    """
    Compute the correlation of star shapes weighted by the residual size.
    """
    short_name = 'rho4'
    long_name = "Rho4 statistics (Correlation of residual star shapes" \
                "weighted by residual size)"
    objects_list = ['star PSF']
    required_quantities = [('ra', 'dec', 'g1', 'g2', 'sigma',
                            'psf_g1', 'psf_g2', 'psf_sigma', 'w')]

    def __call__(self, data, data2=None, random=None, random2=None,
                 config=None, **kwargs):
        new_data = np.rec.fromarrays(
            [data['ra'], data['dec'], data['g1'] - data['psf_g1'],
             data['g2'] - data['psf_g2'], data['w']],
            names=['ra', 'dec', 'g1', 'g2', 'w'])
        if data2 is None:
            data2 = data
        new_data2 = np.rec.fromarrays(
            [data2['ra'], data2['dec'],
             data2['g1'] * (data2['sigma'] - data2['psf_sigma']) / data2[
                 'sigma'],
             data2['g2'] * (data2['sigma'] - data2['psf_sigma']) / data2[
                 'sigma'],
             data2['w']],
            names=['ra', 'dec', 'g1', 'g2', 'w'])
        if random is not None:
            new_random = np.rec.fromarrays(
                [random['ra'], random['dec'],
                 random['g1'] - random['psf_g1'],
                 random['g2'] - random['psf_g2'],
                 random['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])
        else:
            new_random = random
        if random2 is None:
            random2 = random
        if random2 is not None:
            new_random2 = np.rec.fromarrays(
                [random2['ra'], random2['dec'],
                 random2['g1'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['g2'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])
        else:
            new_random2 = random2
        return self.getCF('gg', new_data, new_data2, new_random, new_random2,
                          config=config, **kwargs)


class DESRho5SysTest(new_BaseCorrelationFunctionSysTest):
    r"""Compute correlations.

    The correlation of star shapes weighted by the residual size.
    """

    short_name = 'rho5'
    long_name = "Rho5 statistics (Correlation of star and PSF shapes" \
                "weighted by residual size)"
    objects_list = ['star PSF']
    required_quantities = [('ra', 'dec', 'sigma',
                            'g1', 'g2', 'psf_sigma', 'w')]

    def __call__(self, data, data2=None, random=None, random2=None,
                 config=None, **kwargs):
        new_data = np.rec.fromarrays([data['ra'], data['dec'], data['g1'],
                                      data['g2'], data['w']],
                                     names=['ra', 'dec', 'g1', 'g2', 'w'])
        if data2 is None:
            data2 = data
        new_data2 = np.rec.fromarrays(
            [data2['ra'], data2['dec'], data2['g1'] * (data2['sigma'] - data2[
                'psf_sigma']) / data2['sigma'],
             data2['g2'] * (data2['sigma'] - data2['psf_sigma']) / data2[
                 'sigma'],
             data2['w']],
            names=['ra', 'dec', 'g1', 'g2', 'w'])

        if random is not None:
            new_random = np.rec.fromarrays(
                [random['ra'], random['dec'], random['g1'], random['g2'],
                 random['w']], names=['ra', 'dec', 'g1', 'g2', 'w'])
        else:
            new_random = random
        if random2 is None:
            random2 = random
        if random2 is not None:
            new_random2 = np.rec.fromarrays(
                [random2['ra'], random2['dec'],
                 random2['g1'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['g2'] * (random2['sigma'] - random2['psf_sigma']) /
                 random2['sigma'],
                 random2['w']],
                names=['ra', 'dec', 'g1', 'g2', 'w'])

        else:
            new_random2 = random2
        return self.getCF('gg', new_data, new_data2, new_random, new_random2,
                          config=config, **kwargs)


def rho_stats(starcat_path, output_path, rho_def='HSC',
              print_fun=lambda x: print(x)):
    """Compute and plot the five rho statistics.

    Syntax:

    > python rho_stats.py path/to/starcat rhofile_name [rho_def]
    Where path/to/starcat is the path to the full validation star
    catalog. rhofile_name is the name of the file where all rho_stats
    will be written.
    rho_def is optional and can be either 'HSC' or 'DES'
    depending on the desired definition to use for tho statistics.
    """
    # Read starcat
    starcat = fits.open(starcat_path)

    rho_stats_fun = None

    # Convert HSM flags to 0/1 weights
    hdu_no = 2
    star_flags = starcat[hdu_no].data['FLAG_STAR_HSM']
    psf_flags = starcat[hdu_no].data['FLAG_PSF_HSM']
    w = np.abs(star_flags - 1) * np.abs(psf_flags - 1)
    # w = starcat[2].data['w']

    # Convert to Stile-compatible and change sigmas to R^2 (up to constant)
    stilecat = np.rec.fromarrays(
        [w, starcat[hdu_no].data['RA'], starcat[hdu_no].data['DEC'],
         starcat[hdu_no].data['E1_STAR_HSM'],
         starcat[hdu_no].data['E2_STAR_HSM'],
         starcat[hdu_no].data['SIGMA_STAR_HSM'] ** 2,
         starcat[hdu_no].data['E1_PSF_HSM'],
         starcat[hdu_no].data['E2_PSF_HSM'],
         starcat[hdu_no].data['SIGMA_PSF_HSM'] ** 2],
        names=['w', 'ra', 'dec', 'g1', 'g2', 'sigma', 'psf_g1', 'psf_g2',
               'psf_sigma'])

    # TreeCorr config:
    TreeCorrConfig = {'ra_units': 'degrees', 'dec_units': 'degrees',
                      'max_sep': '3e2', 'min_sep': 5e-1, 'sep_units': 'arcmin',
                      'nbins': 32}

    # Ininitialize all 5 rho stats
    if rho_def == 'HSC':
        rho_stats_fun = [stile.CorrelationFunctionSysTest('Rho{}'.format(j))
                         for j in range(1, 6)]
    elif rho_def == 'DES':
        rho_stats_fun = [Rho1SysTest(), DESRho2SysTest(), DESRho3SysTest(),
                         DESRho4SysTest(), DESRho5SysTest()]

    for rho in rho_stats_fun:
        print_fun(rho.required_quantities)

    # Compute them!
    print_fun(' > Computing rho statistics...')
    start = time.time()
    rho_results = [rho_stat(stilecat, config=TreeCorrConfig) for rho_stat in
                   rho_stats_fun]
    print_fun(' > Done in {}s.'.format(time.time() - start))
    np.save(output_path + 'rho_stat_results.npy', np.array(rho_results))

    # in brackets are DES Y1's ylims
    ylim_l = None  # (5e-9,5e-6)
    ylim_r = None  # (5e-8,1e-5)

    # Plots
    ylims = [ylim_l, ylim_r, ylim_l, ylim_l, ylim_r]
    colors = ['blue', 'blue', 'green', 'orange', 'green']
    markers = ['o', 'o', 's', '^', 's']

    for j, rhores in enumerate(rho_results):
        # [TL] BUG fix on index names
        # meanR [arcmin] -> meanr // sigma_xi -> sigma_xip (?)
        NegDash(rhores['meanr'], rhores['xip'], rhores['sigma_xip'],
                output_path + 'rho_{}.png'.format(j + 1),
                semilogx=True, semilogy=True,
                color=colors[j], capsize=3, fmt=markers[j], alpha=.7,
                ylim=ylims[j],
                xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')

    for j, rhores in enumerate(rho_results):
        if j in [0, 2, 3]:
            # [TL] BUG fix on index names
            # meanR [arcmin] -> meanr // sigma_xi -> sigma_xip (?)
            NegDash(rhores['meanr'], rhores['xip'], rhores['sigma_xip'],
                    semilogx=True, semilogy=True, rho_nb=j + 1,
                    color=colors[j], capsize=3, fmt=markers[j], alpha=.7,
                    ylim=ylims[j],
                    xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')
    plt.legend()
    plt.savefig(output_path + 'lefthand_rhos.pdf')
    plt.close()

    for j, rhores in enumerate(rho_results):
        if j in [1, 4]:
            # [TL] BUG fix on index names
            # meanR [arcmin] -> meanr // sigma_xi -> sigma_xip (?)
            NegDash(rhores['meanr'], rhores['xip'], rhores['sigma_xip'],
                    semilogx=True, semilogy=True, rho_nb=j + 1,
                    color=colors[j], capsize=3, fmt=markers[j], alpha=.7,
                    ylim=ylims[j],
                    xlabel=r'$\theta$ (arcmin)', ylabel=r'$\rho_i(\theta)$')
    plt.legend()
    plt.savefig(output_path + 'righthand_rhos.pdf')
    plt.close()

    starcat.close()

    return None
