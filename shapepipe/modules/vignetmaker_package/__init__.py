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

MASKING : bool, defalut=``False``
    option to modify the SExtractor masking value (``-1e29``) in the vignets
MASK_VALUE : float
    value to use for masking in the vignets
STAMP_SIZE : int
    size of the vignet in pixels, note that this must be an odd integer
COORD : str
    coordinate convention, options are ``PIX`` for pixel coordinates
    or ``SPHE`` for world coordinates
POSITION_PARAMS : list
    list of position parameter value names in the SExtractor output catalogue
MODE : str
    run mode for module, options are ``CLASSIC`` or ``MULTI-EPOCH``
SUFFIX : str or list
    output file name suffix(es)
ME_IMAGE_DIR : list
    input image directories for multi-epoch processing
ME_IMAGE_PATTERN : list
    input file name patterns for multi-epoch processing
ME_LOG_WCS : str
    path to world coordinate system log file (``*sqlite``)

"""

__all__ = ['vignetmaker']
