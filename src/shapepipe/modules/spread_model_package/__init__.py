"""SPREAD MODEL PACKAGE.

This package contains the module for ``spread_model``.

:Author: Axel Guinot

:Parent modules:

- ``sextractor_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
  and ``vignetmaker_runner``

:Input: SExtractor output catalogue, PSFEx output catalogue and vignet files

:Output: Updated SExtractor or new FITS file

Description
===========

This module refines the galaxy sample that will be used for shape measurement
using the spread model method of :cite:`desai:12` and :cite:`mohr:12`.

Module-specific config file entries
===================================

PREFIX : str, optional
    Ouput file prefix
PIXEL_SCALE : float
    Pixel scale in arcseconds
OUTPUT_MODE : str
    Output mode, options are ``add`` to add the outputs to the input
    SExtractor catalogue or ``new`` to generte a new output catalogue

"""

__all__ = ["spread_model"]
