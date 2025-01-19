"""GET IMAGES PACKAGE.

This package contains the module for ``get_images``.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: None or ``find_exposure_runner``

:Input: Image ID list

:Output: Images

Description
===========

Retrieve images. Download using the ``vos`` python library from the
`VOSpace <http://www.ivoa.net/documents/VOSpace>`_ software, or by creating
symbolic links to exising images.

Note that the input for this module are (1) an ASCII file specifying the parent
image IDs, and (2) the images to retrieve.
For (1) the standard config entries ``FILE_PATTERN``, ``FILE_EXT``,
``NUMBERING_SCHEME`` are used. For (2), additional config entries are to be
specified, see below.

Module-specific config file entries
===================================

RETRIEVE : str
    Retrieval method, ``vos`` or ``symlink``
RETRIEVE_OPTIONS : str, optional
    Options to pass to retrieval method
N_TRY : int, optional
    Ff ``RETRIEVE = vos``, this option specifies the number of attempts to
    download an image in case of vos failure; the default value is ``3``
INPUT_PATH : list
    Input path(s) of images
INPUT_FILE_PATTERN : list
    Input file pattern(s) with image ID as dummy template
INPUT_FILE_EXT : list
    Input file extension(s)
OUTPUT_FILE_PATTERN : list
    Output file pattern(s)
INPUT_NUMBERING : str
    Input numbering scheme, python regexp
CHECK_EXISTING_DIR : str, optional
    If given, search this directory (recursively) for existing images, which
    will then not be downloaded
N_EXPECTED : int, optional
    If ``CHECK_EXISTING_DIR`` is given, this option specifies the number of
    expected images; the default value is ``1``

"""

__all__ = ["get_images.py"]
