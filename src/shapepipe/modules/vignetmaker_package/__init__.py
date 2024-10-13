"""VIGNET MAKER PACKAGE.

This package contains the module(s) for ``vignetmaker``.

:Author: Axel Guinot

:Parent module: ``sextractor_runner``

:Input: SExtractor output files

:Output: Vignet FITS files

Description
===========

This module generates vignets (or postage stamps) around galaxy positions from
SExtractor output files and saves the results to FITS files.

Module-specific config file entries
===================================

MASKING : bool, optional
    Option to modify the SExtractor masking value (``-1e29``) in the vignets;
    the default value is ``False``
MASK_VALUE : float
    Value to use for masking in the vignets
STAMP_SIZE : int
    Size of the vignet in pixels, note that this must be an odd integer
COORD : str
    Coordinate convention, options are ``PIX`` for pixel coordinates
    or ``SPHE`` for world coordinates
POSITION_PARAMS : list
    List of position parameter value names in the SExtractor output catalogue
MODE : str
    Run mode for module, options are ``CLASSIC`` or ``MULTI-EPOCH``
PREFIX : str or list
    Output file name prefix(es)
ME_IMAGE_DIR : list
    Module names of last run producing single-exposure flags, images, weights,
    and SExtractor background images, for multi-epoch processing. The specifier
    "last:" is not required
ME_IMAGE_PATTERN : list
    Input file name patterns for flag, image, weight, and SExtractor background
    files, for multi-epoch processing
ME_LOG_WCS : str
    Path to world coordinate system log file (``*sqlite``)

"""

__all__ = ['vignetmaker']
