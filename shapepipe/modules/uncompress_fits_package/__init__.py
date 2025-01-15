"""UNCOMPRESS FITS PACKAGE.

This package contains the module for ``uncompress_fits``.

:Author: Axel Guinot, Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: ``get_images_runner``

:Input: Compressed FITS file

:Ouput: Uncompressed FITS file

Description
===========

This module uncompresses FITS files (e.g. images, weights, flags), and ignores
dummy HDUs with no image data. The compressed FITS weight images contain
a primary empty HDU #0 that will interfere with further processing with
ShapePipe. The uncompressed image is saved as FITS file with a single image
HDU.

Module-specific config file entries
===================================

HDU_DATA : int, optional
    HDU number of input data; the default value is ``0``
OUTPUT_PATTERN : str
    Output file pattern

"""

__all__ = ["uncompress_fits"]
