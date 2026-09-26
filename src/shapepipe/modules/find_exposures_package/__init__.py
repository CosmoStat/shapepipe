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
   Prefix to strip from the exposure filename, e.g. ``simu_image-`` for
   simulated exposures. Leave empty for CFIS exposures, which carry no
   prefix -- the trailing epoch letter (``p``) is a suffix, kept as part
   of the exposure name and unaffected by this key.
"""

__all__ = ["find_exposures.py"]
