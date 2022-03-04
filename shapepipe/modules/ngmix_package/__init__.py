"""NGMIX PACKAGE.

This package contains the module for ``ngmix``.

:Author: Axel Guinot

:Parent modules:
- ``sextractor_runner``
- ``psfex_interp_runner``
- ``vignetmaker_runner``

:Input: galaxy image vignets

:Output: shape catalogue

Description
===========

This module calls the weak-lensing shape measurment and metacalibration
software NGMIX :cite:`sheldon:15`. In particular the GalSim :cite:`rowe:15`
fitting method of NGXMIX is called to obtain galaxy shape measurements.
The metacalibration routines of NGMIX are also called to provide all of the
measurements required to calibrate the shear values.

Module-specific config file entries
===================================

MAG_ZP : float
    Photometric zero point
PIXEL_SCALE : float
    Pixel scale in arcseconds
LOG_WCS : str
    Path to world coordinate system log file (``*sqlite``)
ID_OBJ_MIN : int
    ID of first galaxy object to be processed
ID_OBJ_MAX : int
    ID of last galaxy object to be processed

"""

__all__ = ['ngmix']
