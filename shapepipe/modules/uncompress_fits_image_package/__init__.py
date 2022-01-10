# -*- coding: utf-8 -*-

"""UNCOMRESS FITS

This package contains the module for ``get_uncompress``.

:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: get_images_runner

:Input: compressed FITS file

:Ouput: uncompressed FITS file

:Date: 2020

Description
===========

This module uncompresses FITS files (e.g. images, weights, flags), and ignores
dummy HDUs with no image data. (E.g. the compressed FITS weight images contain
a primary empty HDU #0 that will interfer with further processing by
ShapePipe). The uncompressed image is save as FITS file
with a single image hdu.

Module-specific config file entries
===================================
HDU_DATA : int, optional, default=0
    HDU number of input data
OUTPUT_PATTERN : str
    output file pattern
"""

__all__ = ['uncompress']
