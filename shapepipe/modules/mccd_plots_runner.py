"""MCCD PLOTS RUNNER.

This module is used to generate a series of plots from the merged validation
catalogs.
It plots the Meanshape plots for the merged validation catalog.
It can also plot the rho statistics provided that the required packages are
installed.

:Author: Tobias Liaudat

"""

import warnings

from shapepipe.modules.mccd_package import mccd_plot_utilities as mccd_plots
from shapepipe.modules.module_decorator import module_runner

try:
    import stile
    import stile.stile_utils
    from stile.sys_tests import BaseCorrelationFunctionSysTest
    has_stile = True

except ImportError:
    has_stile = False

try:
    import treecorr
    from treecorr.corr2 import corr2_valid_params
    has_treecorr = True

except ImportError:
    has_treecorr = False

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    # Define the backend for matplotlib
    mpl.use('agg')
    has_mpl = True

except ImportError:
    has_mpl = False


@module_runner(
    version='1.1',
    input_module=['merge_starcat_runner'],
    file_pattern=['full_starcat'],
    file_ext=['.fits'],
    numbering_scheme='-0000000',
    depends=['numpy', 'mccd', 'astropy', 'matplotlib', 'stile', 'treecorr'],
    run_method='serial',
)
def mccd_plots_runner(
    input_file_list,
    run_dirs,
    file_number_string,
    config,
    module_config_sec,
    w_log,
):
    """Define The MCCD Plots Runner."""
    # Input parameters
    if config.has_option(module_config_sec, 'HDU'):
        hdu_no = config.getint(module_config_sec, 'HDU')
    else:
        hdu_no = 2

    # Get parameters for meanshapes plots
    psf_model_type = config.get(module_config_sec, 'PSF')

    if config.has_option(module_config_sec, 'MAX_E'):
        max_e = config.getfloat(module_config_sec, "MAX_E")
    else:
        max_e = None

    if config.has_option(module_config_sec, 'MAX_DE'):
        max_de = config.getfloat(module_config_sec, "MAX_DE")
    else:
        max_de = None

    if config.has_option(module_config_sec, 'MIN_R2'):
        min_r2 = config.getfloat(module_config_sec, "MIN_R2")
    else:
        min_r2 = None

    if config.has_option(module_config_sec, 'MAX_R2'):
        max_r2 = config.getfloat(module_config_sec, "MAX_R2")
    else:
        max_r2 = None

    if config.has_option(module_config_sec, 'MAX_DR2'):
        max_dr2 = config.getfloat(module_config_sec, "MAX_DR2")
    else:
        max_dr2 = None

    x_nb_bins = config.getint(module_config_sec, 'X_GRID')
    y_nb_bins = config.getint(module_config_sec, 'Y_GRID')
    remove_outliers = config.getboolean(module_config_sec, 'REMOVE_OUTLIERS')
    plot_meanshapes = config.getboolean(module_config_sec, 'PLOT_MEANSHAPES')
    plot_histograms = config.getboolean(module_config_sec, 'PLOT_HISTOGRAMS')

    # Get parameters for rho stats plots
    plot_rho_stats = config.getboolean(module_config_sec, 'PLOT_RHO_STATS')
    rho_stat_plot_style = config.get(module_config_sec, 'RHO_STATS_STYLE')

    if config.has_option(module_config_sec, 'RHO_STATS_YLIM_L'):
        str_list = config.getlist(module_config_sec, 'RHO_STATS_YLIM_L')
        ylim_l = [float(s) for s in str_list]
    else:
        ylim_l = None
    if config.has_option(module_config_sec, 'RHO_STATS_YLIM_R'):
        str_list = config.getlist(module_config_sec, 'RHO_STATS_YLIM_R')
        ylim_r = [float(s) for s in str_list]
    else:
        ylim_r = None

    nb_pixel = x_nb_bins, y_nb_bins
    starcat_path = input_file_list[0][0]
    output_path = run_dirs['output'] + '/'

    if plot_meanshapes or plot_histograms:
        if has_mpl:
            mccd_plots.plot_meanshapes(
                starcat_path,
                output_path,
                nb_pixel,
                w_log,
                hdu_no=hdu_no,
                remove_outliers=remove_outliers,
                plot_meanshapes=plot_meanshapes,
                plot_histograms=plot_histograms,
                psf_model_type=psf_model_type,
                max_e=max_e,
                max_de=max_de,
                min_r2=min_r2,
                max_r2=max_r2,
                max_dr2=max_dr2,
            )
        else:
            msg = (
                '[!] In order to plot the Meanshapes the package '
                + '_matplotlib_ has to be correctly imported. This was not'
                + ' the case, so the task is aborted. For the next time make'
                + ' sure the package is installed.'
            )
            warnings.warn(msg)
            w_log.info(msg)

    if plot_rho_stats:
        if has_stile is False or has_treecorr is False:
            msg = (
                '[!] In order to calculate the Rho stats the packages '
                + '_stile_ and _treecorr_ have to be correctly imported.'
                + ' This was not the case, so the rho stat calculation is'
                + 'aborted. For the next time make sure both of the'
                + 'packages are installed.'
            )
            warnings.warn(msg)
            w_log.info(msg)
        elif rho_stat_plot_style != 'HSC' and rho_stat_plot_style != 'DES':
            msg = (
                'The rho stat definition should be HSC or DES. An unknown'
                + ' definition was used. Rho stat calculation aborted.'
            )
            warnings.warn(msg)
            w_log.info(msg)
        else:
            mccd_plots.rho_stats(
                starcat_path,
                output_path,
                rho_def=rho_stat_plot_style,
                hdu_no=hdu_no,
                ylim_l=ylim_l,
                ylim_r=ylim_r,
                print_fun=lambda x: w_log.info(x),
            )

    # No return objects
    return None, None
