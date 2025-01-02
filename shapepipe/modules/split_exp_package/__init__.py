"""SPLIT EXPOSURE PACKAGE.

This package contains the module for ``split_exp``.

:Author: Axel Guinot

:Parent module: ``get_images_runner``

:Input: Single-exposure image, weight, or flag file

:Output: Single-exposure single-HDU file; header numpy binary file (``.npy``)
    if input is an image

Description
===========

This module splits up single-exposure FITS files.
Each CCD of a single-exposure mosaic image is stored in a HDU. This module
saves each HDU in a separate single-exposure single-HDU FITS file, such
that each resulting file only contains data from one CCD.

Module-specific config file entries
===================================

OUTPUT_SUFFIX : list
    Output file name prefixes. Special strings are:

    ``image``:

    1. Header is saved to numpy binary file (``.npy``);
    2. Header WCS is saved in output FITS file header;
    3. Geader WCS coordinates are transformed from pv to sip using the
       ``sip_tpv`` package

    ``flag``: data is save as ``int16``

N_HDU : int
    Number of HDUs (CCDs) of the input mosaic FITS file

"""

__all__ = ["split_exp.py"]
