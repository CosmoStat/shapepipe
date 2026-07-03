r"""MAKE CATALOGUE PACKAGE.

This package contains the module for ``make_cat``.

:Author: Axel Guinot

:Parent modules:

- ``sextractor_runner``
- ``spread_model_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
- ``ngmix_runner``

:Input: SExtractor catalogues, sqlite catalogue

:Output: SExtractor catalogue

Description
===========

This module creates a *final* catalogue combining the output of various
previous module runs. This gathers all relevant information on the measured
galaxies for weak-lensing post-processing. This includes galaxy detection and
basic measurement parameters, the PSF model at galaxy positions, the
spread-model classification, and the shape measurement.

Module-specific config file entries
===================================

SM_DO_CLASSIFICATION : bool, optional
    Adds spread-model star/galaxy classification flag as the column
    ``SPREAD_CLASS`` to the output if ``True``
SM_STAR_THRESH : float, optional
    Threshold :math:`s_{\rm star, thresh}` for star selection; object is
    classified as a star if
    :math:`|x s + 2 \sigma_s | < s_{\textrm{star, thresh}}`
    where :math:`s` is the spread model and :math:`\sigma_s` is the spread
    model error; default value is ``0.003``
SM_GAL_THRESH : float, optional
    Threshold :math:`s_{\rm gal, thresh}` for galaxy selection; object is
    classified as a galaxy if
    :math:`s + 2 \sigma_s  > s_{\textrm{gal, thresh}}` where :math:`s` is the
    spread model and :math:`\sigma_s` is the spread model error; default value
    is ``0.01``
SHAPE_MEASUREMENT_TYPE : list
    Shape measurement method, valid options are ``ngmix`` and/or ``galsim``
SAVE_PSF_DATA : bool, optional
    Save PSF information if ``True``; default value is ``False``
TILE_LIST : str, optional
    Path to list of all tile IDs, used to flag objects in areas of overlap
    between tiles

"""

__all__ = ["make_cat"]
