# -*- coding: utf-8 -*-

"""FIND EXPOSURES PACKAGE

This package contains the module for ``find_exposures``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: get_image_runner

:Input: tile image

:Output: single-exposure name list ASCII file

Description
===========

Identify the exposure images that where co-added to produce the tiles. The image
names are listed in the tile FITS header, which is read by this module to
extract the names.

The output ASCII file contains the image base names (without file extension).

Note that this module is specific for CFIS, in particular the FITS keyword for
identification (``HISTORY``), and the exposure file patterns.

Module-specific entries
=======================
None
"""

__all__ = ['find_exposures.py']
