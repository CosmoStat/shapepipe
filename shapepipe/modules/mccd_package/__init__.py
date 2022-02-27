"""MCCD AUXILIARY FUNCTIONS PACKAGE.

This package contains the modules for ``mccd``.

:Author: Tobias Liaudat


Description
===========

This package contains several modules related with the MCCD PSF model :cite:`liaudat:2021`.


MCCD_PREPROCESSING_RUNNER
=========================

:Author: Tobias Liaudat

:Parent module: ``setools_runner``

:Input: star catalogue

:Output: preprocessed star catalogue 

This module preprocesses an input star catalog as required by the MCCD PSF model.
As the MCCD PSF model builds one PSF model for the entire focal plane it requires that the different 
single-CCD files be merged into a single focal plane file. 
The coordinates are transformed in order to have a global coordinate system.
The minimum number of stars required to keep the CCD can be specified inside the MCCD config 
file with the ``min_n_stars`` parameter.

The input of this runner depends on the value of ``MODE``:
  - ``FIT``: the path to the training star catalogue which will be used to fit the PSF models. 
  - ``VALIDATION``: the path to the validation star catalogue.
  - ``FIT_VALIDATION``: the path to the training and validation star catalogues, in that order.


Module-specific config file entries:

    - CONFIG_PATH: str
        path to the MCCD config file
    - MODE: str
        mode in which the MCCD PSF model will run. 
        It can be ``FIT``, ``VALIDATION`` or ``FIT_VALIDATION``.
        With ``FIT``, the training dataset is preprocessed, with 
        ``VALIDATION`` the validation dataset is preprocessed, and with
        ``FIT_VALIDATION`` both are preprocessed
    - VERBOSE: bool
        verbose mode


MCCD_RUNNERS
============

:Author: Tobias Liaudat

:Parent module: ``mccd_preprocessing_runner``

:Input: preprocessed star catalogue

:Output: fitted PSF models and or PSF validation star catalogue

The three modules ``mccd_fit_runner``, ``mccd_val_runner`` and ``mccd_fit_val_runner`` handle the 
fitting and the validation of the MCCD PSF models.
They all share the same configuration file entries, although the inputs differ as follows
  - ``[MCCD_FIT_RUNNER]``: Has as input the path to the star catalogue that will be used to train the PSF models. 
    The star catalogue should be already preprocessed by ``[MCCD_PREPROCESSING_RUNNER]``.
  - ``[MCCD_VAL_RUNNER]``: Has two inputs. First, the path to the fitted PSF models, and second, 
  the path to the preprocessed validation star catalogue.
  - ``[MCCD_FIT_VAL_RUNNER]``: Has two inputs. First, the path to the preprocessed training star catalogue 
  that will be used to train the PSF models. Second, the path to the preprocessed validation star catalogue. 


Module-specific config file entries:

    - CONFIG_PATH: str
        path to the MCCD config file
    - MODE: str
        mode in which the MCCD PSF model will run. 
        It can be ``FIT``, ``VALIDATION`` or ``FIT_VALIDATION``, and should match the used runner
    - VERBOSE: bool
        verbose mode


MCCD_INTERP_RUNNER
==================

:Author: Tobias Liaudat

:Parent module: ``setools_runner`` and (``mccd_fit_runner`` or ``mccd_fit_val_runner``)

:Input: MCCD PSF model and a catalogue with positions to interpolate the PSF

:Output: MCCD PSF estimates at interpolated positions

This module interpolates the MCCD PSF model to target postions.


Module-specific config file entries:

    - MODE: str
        run mode for module, options are ``CLASSIC`` and ``MULTI-EPOCH``
    - POSITION_PARAMS: (list) str
        list of position parameter value names in the SExtractor output catalogue.
        When running in multi-epoch those position have to be WCS. For example:
        ``MULTI-EPOCH`` case:
        POSITION_PARAMS = XWIN_WORLD,YWIN_WORLD
        ``CLASSIC`` case:
        POSITION_PARAMS = XWIN_IMAGE,YWIN_IMAGE
    - GET_SHAPES: bool
        option to calculate PSF model shapes and save them on the output dictionary
    - PSF_MODEL_DIR: str
        input directories for the fitted MCCD PSF model files
    - PSF_MODEL_PATTERN: str
        pattern of the fitted PSF models
    - PSF_MODEL_SEPARATOR: str
        separator between the PSF model pattern and the model's ID number
    - ME_LOG_WCS: str
        path to world coordinate system log file (``*sqlite``). Only use in ``MULTI-EPOCH`` mode.


MCCD_PLOTS_RUNNER
=================

:Author: Tobias Liaudat

:Parent module: ``merge_starcat_runner``

:Input: merged star catalogue

:Output: several plots

This module generates different plots of the input merged star catalogue.
The plots include meanshapes, where the ellipticity error is plotted as a function of the focal plane positions.
The errors are binned through all the merged catalogue in the grid defined by ``X_GRID`` and ``Y_GRID``. 
The other plots are a histogram of the ellipticity errors and the rho stats plot.


Module-specific config file entries:

    - PSF: str
        PSF model type, options are ``mccd`` or ``psfex``
    - PLOT_MEANSHAPES: bool
        option to produce meanshapes plots
    - X_GRID: int
        number of bins in the X axis of one CCD for the meanshapes' plot grid
    - Y_GRID: int
        number of bins in the Y axis of one CCD for the meanshapes' plot grid
    - MAX_E: float
        max value for ellipticity in the meanshapes plot
    - MAX_DE: float
        max value for the residual ellipticity in the meanshapes plot
    - PLOT_HISTOGRAMS: bool
        option to produce histogram plots for shape errors
    - REMOVE_OUTLIERS: bool
        option to remove validated stars that are outliers in terms of shape before drawing the plots
    - PLOT_RHO_STATS: bool
        option to produce the rho stats plots
    - RHO_STATS_STYLE: str
        style of the rho stats plot, options are ``HSC`` or ``DES``
    - RHO_STATS_YLIM_L: (tuple of) float
        Y-axis limits for left-hand rho stats plot
    - RHO_STATS_YLIM_R: (tuple of) float
        Y-axis limits for right-hand rho stats plot


"""

__all__ = [
    'shapepipe_auxiliary_mccd',
    'mccd_interpolation_script',
    'mccd_plot_utilities',
]

