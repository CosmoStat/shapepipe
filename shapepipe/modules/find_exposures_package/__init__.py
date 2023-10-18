"""FIND EXPOSURES PACKAGE.

This package contains the module for ``find_exposures``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: ``get_images_runner``

:Input: Tile image

:Output: Single-exposure ID list

Description
===========

Identify the exposure images that were co-added to produce the tiles
(stacked image). The image names are listed in the tile FITS header,
which is read by this module to extract the names.

The output ASCII file contains the image base names (without file extension).

Note that this module is specific for CFIS, in particular the FITS keyword for
identification (``HISTORY``), and the exposure file patterns.

Module-specific config file entries
===================================

COLNUM : int
   Column number to find exposure in fits header of tile image for the HISTORY
   string
EXP_PREFIX: str
   Prefix of exposures 
"""

__all__ = ['find_exposures.py']
