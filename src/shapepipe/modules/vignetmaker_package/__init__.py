"""VIGNET MAKER PACKAGE.

This package contains the module(s) for ``vignetmaker``.

:Author: Axel Guinot

:Parent module: ``sextractor_runner``

:Input: SExtractor output files

:Output: Vignet FITS files

Description
===========

This module generates vignets (or postage stamps) around galaxy positions from
SExtractor output files and saves the results to FITS files.

Module-specific config file entries
===================================

MASKING : bool, optional
    Option to modify the SExtractor masking value (``-1e29``) in the vignets;
    the default value is ``False``
MASK_VALUE : float
    Value to use for masking in the vignets
STAMP_SIZE : int
    Size of the vignet in pixels, note that this must be an odd integer
COORD : str
    Coordinate convention, options are ``PIX`` for pixel coordinates
    or ``SPHE`` for world coordinates
POSITION_PARAMS : list
    List of position parameter value names in the SExtractor output catalogue
MODE : str
    Run mode for module, options are ``CLASSIC`` or ``MULTI-EPOCH``
PREFIX : str or list
    Output file name prefix(es)
ME_IMAGE_EXP_DIR : str
    Root of the per-exposure work directories, for multi-epoch processing;
    each runner output directory is resolved beneath it for every exposure
    listed in the input ``exp_numbers`` file
ME_IMAGE_EXP_RUNNERS : list
    Names of the runners producing the single-exposure flags, images, weights
    and SExtractor background images, in the same order as
    ``ME_IMAGE_PATTERN``
ME_IMAGE_PATTERN : list
    Input file name patterns for flag, image, weight, and SExtractor background
    files, for multi-epoch processing

In ``MULTI-EPOCH`` mode the input file list is positional and carries three
entries -- the galaxy catalogue, the world-coordinate-system log
(``*sqlite``) and the ``exp_numbers`` file -- set through ``FILE_PATTERN``
and ``FILE_EXT``, not through module-specific keys.

"""

__all__ = ["vignetmaker"]
