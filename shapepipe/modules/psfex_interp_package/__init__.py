r"""PSFEX INTERPOLATION PACKAGE.

This package contains the module for ``psfex_interp``.

:Author: Axel Guinot

:Parent module: ``psfex_runner`` and ``setools_runner``

:Input: PSFEx PSF model

:Output: PSFEx PSF estimates at galaxy positions

Description
===========

This module interpolates the PSFEx PSF model to galaxy positions.

Module-specific config file entries
===================================

MODE : str
    run mode for module, options are ``CLASSIC``, ``MULTI-EPOCH``
    or ``VALIDATION``
POSITION_PARAMS : list
    list of position parameter value names in the SExtractor output catalogue
GET_SHAPES : bool
    option to compute shapes for the PSF model
STAR_THRESH : int
    threshold of stars under which the PSF is not interpolated
CHI2_THRESH : int
    threshold for chi squared (:math:`$\Chi^2$`)
ME_DOT_PSF_DIR : list
    input directories for PSFEx PSF model files, for multi-epoch processing
ME_DOT_PSF_PATTERN : list
    input file name patterns for PSFEx PSF model files, for multi-epoch
    processing
ME_LOG_WCS : str
    path to world coordinate system log file (``*sqlite``)

"""

__all__ = ['psfex_interp']
