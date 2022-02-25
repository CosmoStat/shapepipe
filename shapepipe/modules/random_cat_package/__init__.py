"""RANDOM CATALOGUE PACKAGE.

This package contains the module for ``random_cat``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: None

:Input: images and masks

:Output: random catalogue FITS file

Description
===========

This module creates a random catalogue, and computes the tile area accounting
for overlapping and masked regions.

Module-specific config file entries
===================================

N_RANDOM : float
    the number of random objects requested on output
DENSITY : bool, optional
    option to interpret the number of random objects per square degree, the
    default is ``False``
TILE_LIST : str, optional
    path to tile IDs for overlap flagging

"""

__all__ = ['random_cat.py']
