r"""MAKE CATALOGUE PACKAGE.

This package contains the module for ``make_cat``.

:Author: Axel Guinot

:Parent modules:

- ``sextractor_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
- ``ngmix_runner``

:Input: SExtractor catalogues, sqlite catalogue

:Output: SExtractor catalogue

Description
===========

This module creates a *final* catalogue combining the output of various
previous module runs. This gathers all relevant information on the measured
galaxies for weak-lensing post-processing. This includes galaxy detection and
basic measurement parameters, the PSF model at galaxy positions, and the
shape measurement.

Module-specific config file entries
===================================

SHAPE_MEASUREMENT_TYPE : list
    Shape measurement method; the only supported option is ``ngmix``
SAVE_PSF_DATA : bool, optional
    Save PSF information if ``True``; default value is ``False``
TILE_LIST : str, optional
    Path to list of all tile IDs, used to flag objects in areas of overlap
    between tiles

"""

__all__ = ["make_cat"]
