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
ME_DOT_PSF_EXP_DIR : str
    Root of the per-exposure work directories, for multi-epoch processing;
    the ``psfex_runner`` output directory of every exposure listed in the
    input ``exp_numbers`` file is resolved beneath it
ME_DOT_PSF_PATTERN : str
    Input file name pattern for PSFEx PSF model files, for multi-epoch
    processing

In ``MULTI-EPOCH`` mode the input file list is positional and carries three
entries -- the galaxy catalogue, the world-coordinate-system log
(``*sqlite``) and the ``exp_numbers`` file -- set through ``FILE_PATTERN``
and ``FILE_EXT``, not through module-specific keys.

"""

__all__ = ["psfex_interp"]
