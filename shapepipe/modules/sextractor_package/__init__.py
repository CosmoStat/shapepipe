"""SEXTRACTOR MODULE.

This package contains the module for ``sextractor``.

:Author: Axel Guinot

:Parent module: ``mask_runner``

:Input: Single-exposure single-CCD image, weight and flag files

:Output: SExtractor output catalogue

Description
===========

This module calls the community standard source extraction software
SExtractor :cite:`bertin:96` to detect objects in the input
images. All of the standard SExtractor options can be used by specifying
the paths to the appropriate configuration files.

Post-processing steps can be carried out to add additional FITS HDUs to tile
catalogues for each exposure, indicating the CCD number.

Module-specific config file entries
===================================

EXEC_PATH : str, optional
    Full path to the SExtractor executable (``sex``) on the system; if not set
    the version controlled SExtractor installation in the ShapePipe environment
    will be used
DOT_SEX_FILE : str
    Full path to the ``.sex`` configuration file for SExtractor
DOT_PARAM_FILE : str
    Full path to the ``.param`` configuration file for SExtractor
DOT_CONV_FILE : str
    Full path to the ``.conv`` kernel file for SExtractor
WEIGHT_IMAGE : bool
    Option to specify if weight files have been provided
FLAG_IMAGE : bool
    Option to specify if flag files have been provided
PSF_FILE : bool
    Option to specify if PSF files have been provided
DETECTION_IMAGE : bool
    Option to specify if detection images have been provided
DETECTION_WEIGHT : bool
    Option to specify if detection weights have been provided
ZP_FROM_HEADER : bool
    Option to to use magnitude zero points defined in the image headers
ZP_KEY : str, optional
    Name of header key for magnitude zero points
BKG_FROM_HEADER : bool
    Option to to use background levels defined in the image headers
BKG_KEY : str, optional
    Name of header key for background levels
CHECKIMAGE : list, optional
    List of header key names corresponding to check images
PREFIX : str, optional
    Output file name prefix
MAKE_POST_PROCESS : bool
    Option to run post-processing steps
LOG_WCS : str, optional
    Path to world coordinate system log file (``*sqlite``)
WORLD_POSITION : list, optional
    List of world coordinates to use to match objects
CCD_SIZE : list, optional
    Size of a CCD in pixels ``[nx, ny]``

"""

__all__ = ["sextractor_script"]
