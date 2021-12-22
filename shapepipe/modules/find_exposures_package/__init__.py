# -*- coding: utf-8 -*-

"""FIND EXPOSURES PACKAGE

This package contains the module for ``find_exposures``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Description: This module identifies the single-exposure image IDs that were co-added to produce a stacked image (tile). The image IDs are extracted from the tile FITS header.

:Parent module: ``get_images_runner``

:Input: tile image

:Output: single-exposure ID list
"""

__all__ = ['find_exposures.py']
