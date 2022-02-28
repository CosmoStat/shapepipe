"""MAKE CATALOGUE PACKAGE.

This package contains the module for ``make_cat``.

:Author: Axel Guinot

:Parent modules:
- ``sextractor_runner``
- ``spread_model_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
- ``ngmix_runner``

:Input: ``SExtractor`` catalogues, ``sqlite`` catalogue

:Output: ``Sextractor`` catalogue

Description
===========

This module creates a 'final' catalogue combining the output of various previous
module runs. This gathers all relevant information on the measured galaxies
for weak-lensing post-processing. This includes galaxy detection and basic
measurement parameters, the PSF model at galaxy positions, the spread-model
classification, and the shape measurement.

Module-specific config file entries
===================================

SM_DO_CLASSIFICATION : bool, optional
    Adds spread-model star/galaxy classification flag as column `SPREAD_CLASS`
    to output if ``True``
SM_STAR_STRESH : float, optional
    Threshold for star selection; object is classified as star if
    |sm + 2 * sm_err < SM_STAR_THRESH where sm (sm_err) is the spread
    model (error); default is 0.003
SM_GAL_STRESH : float, optional
    Threshold for galaxy selection; object is classified as galaxy if
    sm + 2 * sm_err > SM_GAL_THRESH where sm (sm_err) is the spread
    model (error); default is 0.01
SHAPE_MEASUREMENT_TYPE : list
    Shape measurement method, valid is ``ngmix``, ``galsim``
SAVE_PSF_DATA : bool, optional
    Save PSF information if ``True``; default is ``False``
TILE_LIST : str, optional
    Path to list of all tile IDs, used to flag objects in areas of overlap
    between tiles; default is ``None``

"""

__all__ = ['make_cat']
