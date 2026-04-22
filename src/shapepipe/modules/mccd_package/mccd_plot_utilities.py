"""MCCD PLOTS UTILITIES.

This module is used to generate a series of plots from the merged validation
catalogues. It plots the Mean shape plots for the merged validation catalogue.

:Authors: Tobias Liaudat, Martin Kilbinger

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import mccd.mccd_utils as mccd_utils
import numpy as np
from astropy.io import fits

# Define the backend for matplotlib
mpl.use("agg")

# MegaCam -> plt.subplot correspondance, as given by:
"""        'COMMENT Unique detector IDs for MegaCam',
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
"""


def megacam_pos(index):
    """Handle MegaCam Positions.

    Parameters
    ----------
    index : int
        MegaCam position index

    Returns
    -------
    int
        Updated index

    """
    if index < 9:
        # first row - shift by one
        return index + 1
    elif index < 18:
        # second row, non-ears
        return index + 3
    elif index < 27:
        # third row non-ears
        return index + 5
    elif index < 36:
        # fourth row
        return index + 7
    else:
        MegaCamCoords = {36: 11, 37: 21, 38: 22, 39: 32}
        return MegaCamCoords[index]


def megacam_flip(xbins, ybins, ccd_nb, nb_pixel):
    """Flip MegaCam.

    Parameters
    ----------
    xbins : int
        X-axis bins
    ybins : int
        Y-axis bins
    ccd_nb : int
        CCD number
    nb_pixel : int
        Number of pixels

    Returns
    -------
    tuple
        Number of bins in the x and y axes

    """
    if ccd_nb < 18 or ccd_nb in [36, 37]:
        # swap x axis so origin is on top-right
        xbins = nb_pixel[0] - xbins + 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1

    return xbins, ybins


def megacam_flip_2(xbins, ybins, ccd_nb, nb_pixel):
    """Flip MegaCam 2.

    Parameters
    ----------
    xbins : int
        X-axis bins
    ybins : int
        Y-axis bins
    ccd_nb : int
        CCD number
    nb_pixel : int
        Number of pixels

    Returns
    -------
    tuple
        Number of bins in the x and y axes

    """
    if ccd_nb < 18 or ccd_nb in [36, 37]:
        # swap x axis so origin is on top-right
        # xbins = nb_pixel[0] - xbins + 1
        a = 1
    else:
        # swap y axis so origin is on bottom-left
        ybins = nb_pixel[1] - ybins + 1

    return xbins, ybins


def mean_shapes_plot(
    ccd_maps, filename, title="", colorbar_ampl=1.0, wind=None, cmap="bwr"
):
    r"""Mean Shapes Plot.

    Plot mean shapes from CCD maps.

    Parameters
    ----------
    ccd_maps : numpy.ndarray
        CCD maps
    filename : str
        File name
    title : str, optional
        Plot title
    colorbar_ampl : float, optional
        Colour bar amplitude; default is ``1.0``
    wind : numpy.ndarray, optional
        minimum and maximum values for color map; determined from ``ccd_maps``
        if ``None`` (default)
    cmap : str, optional
        Colour map; default is ``bwr``

    """
    # colorbar amplitude
    if wind is None:
        vmax = (
            max(np.nanmax(ccd_maps), np.abs(np.nanmin(ccd_maps)))
            * colorbar_ampl
        )
        vmin = -vmax * colorbar_ampl
    else:
        vmin, vmax = wind[0] * colorbar_ampl, wind[1] * colorbar_ampl

    # create full plot
    fig, axes = plt.subplots(nrows=4, ncols=11, figsize=(18, 12), dpi=400)
    # remove corner axes (above and below ears)
    for j in [0, 10, -1, -11]:
        axes.flat[j].axis("off")
    for ccd_nb, ccd_map in enumerate(ccd_maps):
        ax = axes.flat[megacam_pos(ccd_nb)]
        im = ax.imshow(
            ccd_map.T, cmap=cmap, interpolation="Nearest", vmin=vmin, vmax=vmax
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"rmse={np.sqrt(np.nanmean(ccd_map ** 2)):.3e}", size=8)
    plt.suptitle(title, size=20)  # TODO: fix title
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.savefig(f"{filename}.png")
    plt.close()


def plot_meanshapes(
    starcat_path,
    output_path,
    nb_pixel,
    w_log,
    hdu_no=2,
    remove_outliers=False,
    plot_meanshapes=True,
    plot_histograms=True,
    psf_model_type="mccd",
    max_e=None,
    max_de=None,
    min_r2=None,
    max_r2=None,
    max_dr2=None,
):
    """Plot Mean Shapes.

    Plot mean shapes, sizes, and histograms

    Parameters
    ----------
    starcat_path : str
        Input star and PSF catalogue
    output_path : str
        Output directory for plots
    nb_pixel : numpy.ndarray
        Number of pixels per CCD in x- and y-direction
    w_log : logging.Logger
        Logging instance
    hdu_no : int, optional
        HDU number of data in input FITS file; default is ``2``
    remove_outliers : bool, optional
        Perform outlier rejection if ``True``; default is ``False``
    plot_meanshape : bool, optional
        Plot mean focal plane ellipticities, sizes, and residuals if ``True``;
        default is ``True``
    plot_histograms : bool, optional
        Plot 1D histogram of ellipticities, sizes, and residuals if ``True``;
        default is ``True``
    psf_model_type : str, optional
        PSF model type, options are ``mccd`` or ``psfex``; defualt is ``mccd``
    max_e : float, optional
        maximum value for focal plane ellipticity plots; default is ``None``,
        set according to from data
    max_de : float, optional
        maximum value for focal plane residual ellipticity plots; default is
        ``None``, set according to from data
    min_r2 : float, optional
        minimum value for focal plane size plots, default is ``None``; set
        according to data
    max_r2 : float, optional
        maximum value for focal plane size plots, default is ``None``; set
        according to data
    max_dr2 : float, optional
        maximum value for focal plane residual size plots, default is ``None``;
        set according to data

    """
    # READ FULL STARCAT
    starcat = fits.open(starcat_path, memmap=False)

    auto_colorbar = False
    colorbar_ampl = 1.0
    loc2glob = mccd_utils.Loc2Glob()

    # MegaCam: each CCD is 2048x4612
    grid = np.linspace(0, loc2glob.x_npix, nb_pixel[0] + 1), np.linspace(
        0, loc2glob.y_npix, nb_pixel[1] + 1
    )

    # Flag mask
    star_flags = starcat[hdu_no].data["FLAG_STAR_HSM"]
    psf_flags = starcat[hdu_no].data["FLAG_PSF_HSM"]
    flagmask = np.abs(star_flags - 1) * np.abs(psf_flags - 1)

    # convert sigma to R^2's
    all_star_shapes = np.array(
        [
            starcat[hdu_no].data["E1_STAR_HSM"],
            starcat[hdu_no].data["E2_STAR_HSM"],
            2.0 * starcat[hdu_no].data["SIGMA_STAR_HSM"] ** 2,
        ]
    )
    all_psf_shapes = np.array(
        [
            starcat[hdu_no].data["E1_PSF_HSM"],
            starcat[hdu_no].data["E2_PSF_HSM"],
            2.0 * starcat[hdu_no].data["SIGMA_PSF_HSM"] ** 2,
        ]
    )
    all_CCDs = starcat[hdu_no].data["CCD_NB"]
    all_X = starcat[hdu_no].data["X"]
    all_Y = starcat[hdu_no].data["Y"]

    # Remove stars/PSFs where the measured size is zero
    # Sometimes the HSM shape measurement gives objects with measured
    # size equals to zero without an error Flag.
    bad_stars = abs(all_star_shapes[2, :]) < 0.1
    bad_psfs = abs(all_psf_shapes[2, :]) < 0.1
    size_mask = np.abs(bad_stars) * np.abs(bad_psfs)
    # Remove outlier stars/PSFs
    all_star_shapes = all_star_shapes[:, ~size_mask]
    all_psf_shapes = all_psf_shapes[:, ~size_mask]
    all_CCDs = all_CCDs[~size_mask]
    all_X = all_X[~size_mask]
    all_Y = all_Y[~size_mask]
    flagmask = flagmask[~size_mask]

    # Remove stars/PSFs where the measured size is zero
    # Sometimes the HSM shape measurement gives objects with measured
    # size equals to zero without an error Flag.
    bad_stars = abs(all_star_shapes[2, :]) < 0.1
    bad_psfs = abs(all_psf_shapes[2, :]) < 0.1
    size_mask = np.abs(bad_stars) * np.abs(bad_psfs)
    # Remove outlier stars/PSFs
    all_star_shapes = all_star_shapes[:, ~size_mask]
    all_psf_shapes = all_psf_shapes[:, ~size_mask]
    all_CCDs = all_CCDs[~size_mask]
    all_X = all_X[~size_mask]
    all_Y = all_Y[~size_mask]
    flagmask = flagmask[~size_mask]

    if remove_outliers:
        shape_std_max = 5.0
        # Outlier rejection based on the size
        R2_thresh = shape_std_max * np.std(all_star_shapes[2, :]) + np.mean(
            all_star_shapes[2, :]
        )
        bad_stars = abs(all_star_shapes[2, :]) > R2_thresh
        w_log.info(f"Nb of outlier stars: {np.sum(bad_stars):d}")
        # Remove outlier PSFs
        all_star_shapes = all_star_shapes[:, ~bad_stars]
        all_psf_shapes = all_psf_shapes[:, ~bad_stars]
        all_CCDs = all_CCDs[~bad_stars]
        all_X = all_X[~bad_stars]
        all_Y = all_Y[~bad_stars]
        flagmask = flagmask[~bad_stars]

    e1_res_rmse = np.sqrt(np.mean((all_star_shapes[0, :] - all_psf_shapes[0, :]) ** 2))
    e2_res_rmse = np.sqrt(np.mean((all_star_shapes[1, :] - all_psf_shapes[1, :]) ** 2))
    R2_res_rmse = np.sqrt(np.mean((all_star_shapes[2, :] - all_psf_shapes[2, :]) ** 2))
    w_log.info(f"TOTAL e1 residual RMSE: {e1_res_rmse:.6e}\n")
    w_log.info(f"TOTAL e2 residual RMSE: {e2_res_rmse:.6e}\n")
    w_log.info(f"TOTAL R2 residual RMSE: {R2_res_rmse:.6e}\n")

    # CCDs x star/model x (e1,e2,R2,nstars) x xpos x ypos
    ccd_maps = np.ones((40, 2, 4) + nb_pixel) * np.nan

    for ccd_nb, ccd_map in enumerate(ccd_maps):

        # handle different input catalogue types
        if psf_model_type == "mccd":

            ccd_mask = ((all_CCDs.astype(int) == ccd_nb) * flagmask).astype(
                bool
            )

            # Calculate shift to go from global coordinates to local
            # coordinates
            x_shift, y_shift = loc2glob.shift_coord(ccd_nb)

        elif psf_model_type == "psfex":

            ccd_mask = ((all_CCDs == str(ccd_nb)) * flagmask).astype(bool)

            # No shift required for PSFEx
            x_shift, y_shift = 0, 0

        else:

            raise ValueError(f"Invalid psf model type {psf_model_type}")

        star_shapes = all_star_shapes[:, ccd_mask]
        psf_shapes = all_psf_shapes[:, ccd_mask]

        xs_loc, ys_loc = all_X[ccd_mask] - x_shift, all_Y[ccd_mask] - y_shift

        if psf_model_type == "mccd":
            # swap axes to match CCD orientation and origin convention
            ys_loc = loc2glob.y_npix - ys_loc + 1

        # digitalize into bins
        xbins = np.digitize(xs_loc, grid[0])
        ybins = np.digitize(ys_loc, grid[1])

        if psf_model_type == "psfex":
            xbins, ybins = megacam_flip(xbins, ybins, ccd_nb, nb_pixel)

        for xb in range(nb_pixel[0]):
            for yb in range(nb_pixel[1]):
                bin_star_shapes = star_shapes[:, (xbins == xb + 1) * (ybins == yb + 1)]
                bin_psf_shapes = psf_shapes[:, (xbins == xb + 1) * (ybins == yb + 1)]
                ccd_map[0, :3, xb, yb] = np.mean(bin_star_shapes, axis=1)
                ccd_map[1, :3, xb, yb] = np.mean(bin_psf_shapes, axis=1)
                ccd_map[:, 3, xb, yb] = bin_star_shapes.shape[1]

    if plot_meanshapes:
        # e_1
        if max_e:
            vmax = max_e
        else:
            vmax = max(
                np.nanmax(ccd_maps[:, :, 0]),
                np.abs(np.nanmin(ccd_maps[:, :, 0])),
            )
        vmin = -vmax
        wind = [vmin, vmax]
        title = (
            f"e_1 (stars), std={np.nanstd(ccd_maps[:, 0, 0]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 0, 0])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 0, 0], output_path + "e1s", title, wind=wind
        )
        mean_shapes_plot(ccd_maps[:, 0, 0], output_path + "e1s", title, wind=wind)

        title = (
            f"e_1 (model), std={np.nanstd(ccd_maps[:, 1, 0]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 1, 0])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 1, 0], output_path + "e1m", title, wind=wind
        )
        mean_shapes_plot(ccd_maps[:, 1, 0], output_path + "e1m", title, wind=wind)

        if auto_colorbar:
            wind = None
        e1_res = ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0]
        e1_res = e1_res[~np.isnan(e1_res)]
        rmse_e1 = np.sqrt(np.mean(e1_res**2))
        w_log.info(f"Bins: e1 residual RMSE: {rmse_e1:.6f}\n")
        if max_de:
            vmax = max_de
        else:
            vmax = np.nanmax(abs(ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0]))
        vmin = -vmax
        wind = [vmin, vmax]
        title = (
            f"e_1 res, rmse={rmse_e1:.5e}\nvmax={vmax:.4e} , "
            + f"std={np.nanstd(ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0]):.5e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 0, 0] - ccd_maps[:, 1, 0],
            output_path + "e1res",
            title,
            wind=wind,
            colorbar_ampl=colorbar_ampl,
        )

        # e_2
        if max_e:
            vmax = max_e
        else:
            vmax = max(
                np.nanmax(ccd_maps[:, :, 1]),
                np.abs(np.nanmin(ccd_maps[:, :, 1])),
            )
        vmin = -vmax
        wind = [vmin, vmax]
        title = (
            f"e_2 (stars), std={np.nanstd(ccd_maps[:, 0, 1]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 0, 1])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 0, 1], output_path + "e2s", title, wind=wind
        )
        mean_shapes_plot(ccd_maps[:, 0, 1], output_path + "e2s", title, wind=wind)
        title = (
            f"e_2 (model), std={np.nanstd(ccd_maps[:, 1, 1]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 1, 1])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 1, 1], output_path + "e2m", title, wind=wind
        )
        mean_shapes_plot(ccd_maps[:, 1, 1], output_path + "e2m", title, wind=wind)

        if auto_colorbar:
            wind = None
            colorbar_ampl = 1.0

        e2_res = ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1]
        e2_res = e2_res[~np.isnan(e2_res)]
        rmse_e2 = np.sqrt(np.mean(e2_res**2))
        w_log.info(f"Bins: e2 residual RMSE: {rmse_e2:.6f}\n")
        if max_de:
            vmax = max_de
        else:
            vmax = np.nanmax(abs(ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1]))
        vmin = -vmax
        wind = [vmin, vmax]
        title = (
            f"e_2 res, rmse={rmse_e2:.5e}\nvmax={vmax:.4e} , "
            + f"std={np.nanstd(ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1]):.5e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 0, 1] - ccd_maps[:, 1, 1],
            output_path + "e2res",
            title,
            wind=wind,
            colorbar_ampl=colorbar_ampl,
        )

        # R^2
        if min_r2:
            vmin = min_r2
        else:
            vmin = 0
        if max_r2:
            vmax = max_r2
        else:
            vmax = np.nanmax(ccd_maps[:, :, 2])
        wind = [vmin, vmax]
        colorbar_ampl = 1
        title = (
            f"R_2 (stars), std={np.nanstd(ccd_maps[:, 0, 2]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 0, 2])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 0, 2],
            output_path + "R2s",
            title,
            wind=wind,
            cmap="Reds",
        )
        title = (
            f"R_2 (model), std={np.nanstd(ccd_maps[:, 1, 2]):.5e}\n"
            + f"vmax={np.nanmax(abs(ccd_maps[:, 1, 2])):.4e}"
        )
        mean_shapes_plot(
            ccd_maps[:, 1, 2],
            output_path + "R2m",
            title,
            wind=wind,
            cmap="Reds",
        )

        if auto_colorbar:
            wind = [
                0,
                np.nanmax(
                    np.abs(
                        (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2])
                        / ccd_maps[:, 0, 2]
                    )
                ),
            ]
            colorbar_ampl = 1.0
        R2_res = (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]
        R2_res = R2_res[~np.isnan(R2_res)]
        rmse_r2 = np.sqrt(np.mean(R2_res**2))
        w_log.info(f"Bins: R2 residual RMSE: {rmse_r2:.6f}\n")
        if max_dr2:
            vmax = max_dr2
        else:
            vmax = np.nanmax(
                abs((ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2])
            )
        wind = [0, vmax]
        std_title = np.nanstd(
            (ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]
        )
        title = (
            f"∆(R_2)/R_2 res, rmse={rmse_r2:.5e}\nvmax={vmax:.4e} , "
            + f"std={std_title:.5e}"
        )
        if remove_outliers:
            title = "Outliers removed\n" + title

        mean_shapes_plot(
            np.abs((ccd_maps[:, 0, 2] - ccd_maps[:, 1, 2]) / ccd_maps[:, 0, 2]),
            output_path + "R2res",
            title,
            wind=wind,
            colorbar_ampl=colorbar_ampl,
            cmap="Reds",
        )

        # nstars
        wind = (0, np.max(ccd_maps[:, 0, 3]))
        title = f"Number of stars\nTotal={np.nansum(ccd_maps[:, 0, 3]):.0f}"
        mean_shapes_plot(
            ccd_maps[:, 0, 3],
            f"{output_path}nstar",
            title,
            wind=wind,
            cmap="magma",
        )

    # Histograms
    if plot_histograms:
        hist_bins = 50
        plt.figure(figsize=(12, 6), dpi=300)
        plt.hist(
            all_star_shapes[0, :],
            bins=hist_bins,
            range=[-0.2, 0.2],
            label="stars",
            alpha=0.5,
        )
        plt.hist(
            all_psf_shapes[0, :],
            bins=hist_bins,
            range=[-0.2, 0.2],
            label="PSFs",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("e1", fontsize=24)
        plt.savefig(f"{output_path}e1_hist.png")
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = all_star_shapes[0, :] - all_psf_shapes[0, :]
        plt.hist(
            data_hist,
            bins=hist_bins,
            range=[np.min(data_hist), np.max(data_hist)],
            label="err(star - psf)",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("e1 err", fontsize=24)
        plt.savefig(output_path + "err_e1_hist.png")
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        plt.hist(
            all_star_shapes[1, :],
            bins=hist_bins,
            range=[-0.2, 0.2],
            label="stars",
            alpha=0.5,
        )
        plt.hist(
            all_psf_shapes[1, :],
            bins=hist_bins,
            range=[-0.2, 0.2],
            label="PSFs",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("e2", fontsize=24)
        plt.savefig(output_path + "e2_hist.png")
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = all_star_shapes[1, :] - all_psf_shapes[1, :]
        plt.hist(
            data_hist,
            bins=hist_bins,
            range=[np.min(data_hist), np.max(data_hist)],
            label="err(star - psf)",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("e2 err", fontsize=24)
        plt.savefig(output_path + "err_e2_hist.png")
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        mean_R2 = np.mean(all_star_shapes[2, :])
        wind = [mean_R2 - 4, mean_R2 + 4]
        plt.hist(
            all_star_shapes[2, :],
            bins=hist_bins,
            range=wind,
            label="stars",
            alpha=0.5,
        )
        plt.hist(
            all_psf_shapes[2, :],
            bins=hist_bins,
            range=wind,
            label="PSFs",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("R2", fontsize=24)
        plt.savefig(output_path + "R2_hist.png")
        plt.close()

        plt.figure(figsize=(12, 6), dpi=300)
        data_hist = (
            all_star_shapes[2, :] - all_psf_shapes[2, :]
        ) / all_star_shapes[2, :]
        plt.hist(
            data_hist,
            bins=hist_bins,
            range=[np.min(data_hist), np.max(data_hist)],
            label="err(star - psf)/star",
            alpha=0.5,
        )
        plt.legend(loc="best", fontsize=16)
        plt.title("R2 err", fontsize=24)
        plt.savefig(output_path + "err_R2_hist.png")
        plt.close()

    starcat.close()
