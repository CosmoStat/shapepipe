"""SEXTRACTOR MODULE.

This package contains the module for ``sextractor``.

:Author: Axel Guinot

:Parent module(s): ``get_images_runner' or ``split_exp_runner``;
    ``uncompress_fits_runner``; ``mask_runner``

:Input: single-exp_single-CCD image[, weights] [, flags] [, PSF file] [, detection image] [, detection weight]

:Output: ``SExtractor`` catalogue

Description
===========

This module detects objects in input images.
The detection is done with the software ``SExtractor``, which is installed by ``ShapePipe``
by default.

In addition to the input image, configuration flags indicate additional input
file. These can be a weight image, a flag (mask) image, a PSF file, a distinct
image for detection (for ``SExtractor`` in dual-image mode) that is different
from the measurement image, and a distinct detection weight.


Module-specific config file entries
===================================

EXEC_PATH : str
    path to the ``SExtractor`` executable
DOT_SEX_FILE : str
    path to main ``SExtractor`` configuration file
DOT_PARAM_FILE : str
    path to ``SExtractor`` output parameter file
DOT_CONV_FILE : str
    path to ``SExtractor`` convolution kernel file
WEIGHT_IMAGE : bool
    use input weight image if True
FLAG_IMAGE : bool
    use input flag image if True
PSF_FILE : bool
    use input PSF file if True
DETECTION_IMAGE : str
    use distinct image for detection (run ``SExtractor`` in
    dual-image mode) if True
DETECTION_WEIGHT : str
    distinct weight image for detection (run ``SExtractor``
    in dual-image mode) if True
ZP_FROM_HEADER : bool
    photometry zero-point is read from image header if True
    use ``DOT_SEX_FILE:MAG_ZEROPOINT`` if False
ZP_KEY : str
    if ZP_FROM_HEADER is True, zero-point FITS head er key name;
BKG_FROM_HEADER : False
    background value will be read from header if True; in that
    case the ``SExtractor`` value ``BACK_TYPE`` is set to ``MANUAL``
BKG_KEY : str
    background value FITS header key name; only required if
    ``BKG_FROM_HEADER`` is True
CHECKIMAGE : str, optional
    ``SExtractor`` checkimage flag, can be a list of ``BACKGROUND``,
    ``BACKGROUND_RMS``, ``INIBACKGROUND``, ``MINIBACK_RMS``, ``-BACKGROUND``, 
    ``FILTERED``, ``OBJECTS``, ``-OBJECTS``, ``SEGMENTATION``, ``APERTURES``;
    if not given, no image check is carried out
SUFFIX : str, optional
    image output name prefix; None if not given
MAKE_POST_PROCESS : bool
    if True, perform "post-processing" steps; for stacked tile images, add
    single-exposure information to additional output FITS file HDUs; st to
    False for single-exposure images
LOG_WCS : str
    path to .sqlite file with single-exposure WCS header information; only
    only required if ``MAKE_POST_PROCESS`` is True
WORLD_POSITION : list(2) of str
    X- and Y-position World coordinate keywords (``SExtractor`` output parameters);
    only required if ``MAKE_POST_PROCESS`` is True
CCD_SIZE : list(4) of int
    xmin, xmax, ymin, ymax; pixel coordinate range (x, y) of the CCD, with
    xmin < x < xmax and ymin < y < ymax;
    only required if ``MAKE_POST_PROCESS`` is True

"""

__all__ = ['sextractor_script']
