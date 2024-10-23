"""RANDOM CATALOGUE PACKAGE.

This package contains the module for ``random_cat``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: None

:Input: Images and masks

:Output: Random catalogue FITS file

Description
===========

This module creates a random catalogue, and computes the tile area accounting
for overlapping and masked regions.

Module-specific config file entries
===================================

N_RANDOM : float
    The number of random objects requested on output
DENSITY : bool, optional
    Option to interpret the number of random objects per square degree; the
    default is ``False``
SAVE_MASK_AS_HEALPIX : bool
    Output healpix mask if ``True``
HEALPIX_OUT_FILE_BASE : str, optional
    Output halpix mask file base name; used only if SAVE_MASK_AS_HEALPIX is
    ``True``
HEALPIX_OUT_NSIDE : int, optional
    Output healpix mask nside; used only if SAVE_MASK_AS_HEALPIX is ``True``

"""

__all__ = ["random_cat.py"]
