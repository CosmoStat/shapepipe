"""SEXTRACTOR MODULE.

This package contains the module for ``sextractor``.

:Author: Axel Guinot

:Parent module: ``mask_runner``

:Input: single-exposure image, weight and flag files

:Output: SExtractor output catalogue

Description
===========

This module calls the community standard source extraction software
SExtractor :cite:`bertin:96` to identify stars and galaxies in the input
images. All of the standard SExtractor options can be used by specifying
the paths to the appropriate configuration files.

Post-processing steps can be carried out to add additional FITS HDUs to tile
catalogues for each exposure, indicating the CCD number.

Module-specific config file entries
===================================

EXEC_PATH : str, optional
    full path to the SExtractor executable on the system; if not set the
    version controlled SExtractor installation in the ShapePipe environment
    will be used
DOT_SEX_FILE : str
    full path to the ``.sex`` configuration file for SExtractor
DOT_PARAM_FILE : str
    full path to the ``.param`` configuration file for SExtractor
DOT_CONV_FILE : str
    full path to the ``.conv`` kernel file for SExtractor
WEIGHT_IMAGE : bool
    option to specify if weight files have been provided
FLAG_IMAGE : bool
    option to specify if flag files have been provided
PSF_FILE : bool
    option to specify if PSF files have been provided
DETECTION_IMAGE : bool
    option to specify if detection images have been provided
DETECTION_WEIGHT : bool
    option to specify if detection weights have been provided
ZP_FROM_HEADER : bool
    option to to use zero points defined in the image headers
ZP_KEY : str, optional
    name of header key for zero points
BKG_FROM_HEADER : bool
    option to to use background levels defined in the image headers
BKG_KEY : str, optional
    name of header key for background levels
CHECKIMAGE : list, optional
    list of header key names corresponding to check images
PREFIX : str, optional
    output file name prefix
MAKE_POST_PROCESS : bool
    option to run post-processing steps
LOG_WCS : str, optional
    path to world coordinate system log file (``*sqlite``)
WORLD_POSITION : list, optional
    list of world coordinates to use to match objects
CCD_SIZE : list, optional
    size of a CCD in pixels [nx, ny]

"""

__all__ = ['sextractor_script']
