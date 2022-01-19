# -*- coding: utf-8 -*-

"""GET IMAGES PACKAGE

This package contains the module for ``get_images``

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: None or find_exposure_runner

:Input: image ID list

:Output: images

Description
===========

Retrieve images. Download using the ``vos`` python library from the ``VOSpace``
(http://www.ivoa.net/documents/VOSpace) software, or by creating symbolic links
to exising images.

Note that the input for this module are (1) a ASCII file specifying the parent
image IDs, and (2) the images to retrieve.
For (1) the standard config entries ``FILE_PATTERN``, ``FILE_EXT``, ``NUMBERING_SCHEME``
are used. For (2), additional config entries are to be specified, see below.

Module-specific config file entries
===================================
RETRIEVE : str
    retrieval method, ``vos`` or ``symlink``
RETRIEVE_OPTIONS : str, optional
    options to pass to retrieval method
N_TRY : int, optional, default=3
    if RETRIEVE=``vos``, number of attempts to download image in case of vos
    failure
INPUT_PATH : (list of) str
    input path(s) of images
INPUT_FILE_PATTERN : (list of) str
    input file pattern(s) with image ID as dummy template
INPUT_FILE_EXT : (list of) str
    input file extension(s)
    input path(s) of images
OUTPUT_FILE_PATTERN : (list of) str
    output file pattern
INPUT_NUMBERING : str
    input numbering scheme, python regexp
CHECK_EXISTING_DIR : str, optional
    if given, check this directory for existing images, which will then not be
    downloaded
N_EXPECTED : int, optional, default=1
    if CHECK_EXISTING_DIR is given, number of expected images

"""

__all__ = ['get_images.py']
