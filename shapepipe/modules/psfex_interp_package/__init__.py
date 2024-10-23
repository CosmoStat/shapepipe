r"""PSFEX INTERPOLATION PACKAGE.

This package contains the module for ``psfex_interp``.

:Author: Axel Guinot

:Parent modules:

- ``psfex_runner``
- ``setools_runner``

:Input: PSFEx PSF model

:Output: PSFEx PSF estimates at interpolated positions

Description
===========

This module interpolates the PSFEx PSF model.

Module-specific config file entries
===================================

MODE : str
    Run mode for module, options are ``CLASSIC``, ``MULTI-EPOCH``
    or ``VALIDATION``
POSITION_PARAMS : list
    List of position parameter value names in the SExtractor output catalogue
GET_SHAPES : bool
    Option to compute shapes for the PSF model
STAR_THRESH : int
    Threshold of stars under which the PSF is not interpolated
CHI2_THRESH : int
    Threshold for chi squared (:math:`\chi^2`)
ME_DOT_PSF_DIR : str
    Module name of last run producing PSFEx PSF model files, for multi-epoch
    processing. The specifier "last:" is not required
ME_DOT_PSF_PATTERN : str
    Input file name pattern for PSFEx PSF model files, for multi-epoch
    processing
ME_LOG_WCS : str
    Path to world coordinate system log file (``*sqlite``)

"""

__all__ = ['psfex_interp']
