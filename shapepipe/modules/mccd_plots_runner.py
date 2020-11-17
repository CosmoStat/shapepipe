# -*- coding: utf-8 -*-

"""MCCD PLOTS RUNNER.

This module is used to generate a series of plots from the merged validation
catalogs.
It plots the Meanshape plots for the merged validation catalog.
It can also plot the rho statistics provided that the required packages are
installed.

:Author: Tobias Liaudat

"""

from shapepipe.modules.module_decorator import module_runner
from shapepipe.modules.MCCD_package import mccd_plot_utilities as mccd_plots
import warnings


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


@module_runner(input_module=['mccd_merge_starcat_runner'], version='1.0',
               file_pattern=['full_starcat'],
               file_ext=['.fits'], numbering_scheme='-0000000',
               depends=['numpy', 'mccd', 'astropy', 'matplotlib', 'stile',
                        'treecorr'],
               run_method='serial')
def mccd_plots_runner(input_file_list, run_dirs, file_number_string,
                      config, w_log):
    # Get parameters for meanshapes plots
    x_nb_bins = config.getint('MCCD_PLOTS_RUNNER', 'X_GRID')
    y_nb_bins = config.getint('MCCD_PLOTS_RUNNER', 'Y_GRID')
    remove_outliers = config.getboolean('MCCD_PLOTS_RUNNER', 'REMOVE_OUTLIERS')
    plot_meanshapes = config.getboolean('MCCD_PLOTS_RUNNER', 'PLOT_MEANSHAPES')
    plot_histograms = config.getboolean('MCCD_PLOTS_RUNNER', 'PLOT_HISTOGRAMS')

    # Get parameters for Rho stats plots
    plot_rho_stats = config.getboolean('MCCD_PLOTS_RUNNER', 'PLOT_RHO_STATS')
    rho_stat_plot_style = config.get('MCCD_PLOTS_RUNNER', 'RHO_STATS_STYLE')

    nb_pixel = x_nb_bins, y_nb_bins
    starcat_path = input_file_list[0][0]
    output_path = run_dirs['output'] + '/'

    if plot_meanshapes or plot_histograms:
        if has_mpl:
            mccd_plots.plot_meanshapes(starcat_path, output_path, nb_pixel,
                                       w_log, remove_outliers,
                                       plot_meanshapes,
                                       plot_histograms)
        else:
            msg = "[!] In order to plot the Meanshapes the package " \
                  "_matplotlib_ has to be correctly imported." \
                  "This was not the case, so the task is" \
                  "aborted. For the next time make sure the package is" \
                  "installed."
            warnings.warn(msg)
            w_log.info(msg)

    if plot_rho_stats:
        if has_stile is False or has_treecorr is False:
            msg = "[!] In order to calculate the Rho stats the packages " \
                  "_stile_ and _treecorr_ have to be correctly imported." \
                  "This was not the case, so the rho stat calculation is" \
                  "aborted. For the next time make sure both of the" \
                  "packages are installed."
            warnings.warn(msg)
            w_log.info(msg)
        elif rho_stat_plot_style != 'HSC' and rho_stat_plot_style != 'DEC':
            msg = "The rho stat definition should be HSC or DEC." \
                  "An unknown definition was used. Rho stat calculation" \
                  "aborted."
            warnings.warn(msg)
            w_log.info(msg)
        else:
            mccd_plots.rho_stats(starcat_path, output_path,
                                 rho_def=rho_stat_plot_style,
                                 print_fun=lambda x: w_log.info(x))

    return None, None
