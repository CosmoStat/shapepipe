"""MCCD AUXILIARY FUNCTIONS PACKAGE.

This package contains the modules for ``mccd``.

:Author: Tobias Liaudat


Description
===========

This package contains several modules related to the MCCD PSF model
:cite:`liaudat:2021`.


MCCD_PREPROCESSING_RUNNER
=========================

:Author: Tobias Liaudat

:Parent module: ``setools_runner``

:Input: Star catalogue

:Output: Preprocessed star catalogue

This module preprocesses an input star catalogue as required by the MCCD PSF
model. As the MCCD PSF model builds one PSF model for the entire focal plane it
requires that the different single-CCD files be merged into a single focal
plane file. The coordinates are transformed in order to have a global
coordinate system. The minimum number of stars required to keep the CCD can be
specified inside the MCCD config file with the ``min_n_stars`` parameter.

The input of this runner depends on the value of ``MODE``.

- ``FIT``: The input star catalogue path which will be used to fit the PSF
  models.
- ``VALIDATION``: The input star catalogue path which will be used to
  validate the PSF models.
- ``FIT_VALIDATION``: Two input star catalogue paths which will be used to
  fit and validate the PSF models, respectively.


Module-specific config file entries
-----------------------------------

CONFIG_PATH: str
    Path to the MCCD config file
MODE: str
    Mode in which the MCCD PSF model will run.
    It can be ``FIT``, ``VALIDATION`` or ``FIT_VALIDATION``;
    with ``FIT``, the training dataset is preprocessed, with
    ``VALIDATION`` the validation dataset is preprocessed, and with
    ``FIT_VALIDATION`` both are preprocessed
VERBOSE: bool
    Verbose mode


MCCD_RUNNERS
============

:Author: Tobias Liaudat

:Parent module: ``mccd_preprocessing_runner``

:Input: Preprocessed star catalogue

:Output: Fitted PSF models and or PSF validation star catalogue

The three modules ``mccd_fit_runner``, ``mccd_val_runner`` and
``mccd_fit_val_runner`` handle the fitting and the validation of the MCCD PSF
models.

They all share the same configuration file entries, although the inputs differ
as follows.

- ``[MCCD_FIT_RUNNER]``: Has as input the path to the star catalogue that
  will be used to train the PSF models. The star catalogue should be already
  preprocessed by ``[MCCD_PREPROCESSING_RUNNER]``.
- ``[MCCD_VAL_RUNNER]``: Has two inputs. First, the path to the fitted PSF
  models, and second, the path to the preprocessed validation star catalogue.
- ``[MCCD_FIT_VAL_RUNNER]``: Has two inputs. First, the path to the
  preprocessed training star catalogue that will be used to train the PSF
  models. Second, the path to the preprocessed validation star catalogue.

Module-specific config file entries
-----------------------------------

CONFIG_PATH: str
    Path to the MCCD config file
MODE: str
    Mode in which the MCCD PSF model will run; it can be ``FIT``,
    ``VALIDATION`` or ``FIT_VALIDATION``, and should match the used runner
VERBOSE: bool
    Verbose mode

MCCD_INTERP_RUNNER
==================

:Author: Tobias Liaudat

:Parent module: ``setools_runner`` and
    (``mccd_fit_runner`` or ``mccd_fit_val_runner``)

:Input: MCCD PSF model and a catalogue with positions to interpolate the PSF

:Output: MCCD PSF estimates at interpolated positions

This module interpolates the MCCD PSF model to target postions.

Module-specific config file entries
-----------------------------------

MODE: str
    Run mode for module, options are ``CLASSIC`` and ``MULTI-EPOCH``
POSITION_PARAMS: (list) str
    List of position parameter value names in the SExtractor output catalogue,
    when running in multi-epoch these positions have to be WCS; e.g.

    - ``MULTI-EPOCH`` case:
      ``POSITION_PARAMS = XWIN_WORLD,YWIN_WORLD``
    - ``CLASSIC`` case:
      ``POSITION_PARAMS = XWIN_IMAGE,YWIN_IMAGE``
GET_SHAPES: bool
    Calculate PSF model shapes and save to output if ``True``
PSF_MODEL_DIR: str
    Module name of last run producing the fitted MCCD PSF model files.
    The specifier "last:" is not required
PSF_MODEL_PATTERN: str
    Pattern of the fitted PSF models
PSF_MODEL_SEPARATOR: str
    Separator between the PSF model pattern and the model's ID number
ME_LOG_WCS: str
    Path to world coordinate system log file (``*sqlite``); only use in
    ``MULTI-EPOCH`` mode

MCCD_PLOTS_RUNNER
=================

:Author: Tobias Liaudat

:Parent module: ``merge_starcat_runner``

:Input: Merged star catalogue

:Output: Several plots

This module generates different series of plots of the input merged star
catalogue.

One series of plots are the _mean shapes_, which is the PSF ellipticity and
size and their residuals as a function of the focal plane positions. Values and
residuals are combined from all merged catalogues, and binned in the focal
plane according to ``X_GRID`` and ``Y_GRID``.

Another series of plots are a histogram of the ellipticity errors.

A third series are the rho statistics plot, which are correlation functions of
various combinations of PSF ellipticity, size, and their residuals,
see :cite:`rowe:10` and :cite:`jarvis:16`.

Module-specific config file entries
-----------------------------------

PSF: str
    PSF model type, options are ``mccd`` or ``psfex``
PLOT_MEANSHAPES: bool
    Option to produce _mean shapes_ plots
X_GRID: int
    Number of bins in the X axis of one CCD for the _mean shapes_ plot grid
Y_GRID: int
    Number of bins in the Y axis of one CCD for the _mean shapes_ plot grid
MAX_E: float
    Max value for ellipticity in the _mean shapes_ plot
MAX_DE: float
    Max value for the residual ellipticity in the _mean shapes_ plot
PLOT_HISTOGRAMS: bool
    Option to produce histogram plots for shape errors
REMOVE_OUTLIERS: bool
    Option to remove validated stars that are outliers in terms of shape
    before drawing the plots
PLOT_RHO_STATS: bool
    Option to produce the rho stats plots
RHO_STATS_STYLE: str
    Style of the rho stats plot, options are ``HSC`` or ``DES``
RHO_STATS_YLIM_L: tuple
    Y-axis limits for left-hand rho stats plot
RHO_STATS_YLIM_R: tuple
    Y-axis limits for right-hand rho stats plot

"""

__all__ = [
    'shapepipe_auxiliary_mccd',
    'mccd_interpolation_script',
    'mccd_plot_utilities',
]
