"""MASK EXT MODULE.

This package contains the module for ``mask_ext``.

:Author: Cail Daley

:Parent module: ``split_exp_runner`` (exposures) or ``get_images_runner`` /
                ``uncompress_fits_runner`` (tiles), or None

:Input: Single image (tile or single-exposure single-CCD) and, optionally, an
        external instrument flag file

:Output: Per-image pixel flag file

Description
===========

This module produces the ShapePipe per-image pixel flag file
``<PREFIX>_flag<num>.fits`` by rasterizing an **external healsparse mask** onto
the image pixel grid, rather than generating masks internally (the job of the
``mask`` module). It is the ShapePipe consumer of the unified UNIONS healsparse
mask products (PhotoPipe footprint + bright stars, MaxiMask, manual galaxy
masks), merged and stored as ``healsparse.HealSparseMap`` files.

For each image the module:

1. Reads the image header into an ``astropy.wcs.WCS``.
2. Evaluates the pixel grid to (RA, Dec) in row chunks (bounded memory).
3. Queries ``HealSparseMap.get_values_pos`` at every pixel centre.
4. Maps healsparse bit values to ShapePipe flag values through the
   config-driven ``BIT_FLAG_MAP``, producing an ``int16`` flag image.
5. Optionally sums in an external instrument flag image (``USE_EXT_FLAG``),
   matching the combination semantics of the ``mask`` module.
6. Writes ``<PREFIX>_flag<num>.fits`` carrying the image WCS in its header.

Everything downstream (SExtractor ``IMAFLAGS_ISO``, setools star selection,
vignetmaker, ngmix) consumes this artifact unchanged; the module is a drop-in
replacement for ``mask`` at the flag-image contract.

The same module serves tiles and exposure CCDs: only the input image (hence
WCS) differs. The healsparse file, its bit meanings, and all flag values live
**only in config** — the mask products evolve and are swapped at zero code cost.

Module-specific config file entries
===================================

MASK_PATH : str
    Path to the healsparse mask file (``.hsp``/``.fits``); band-agnostic
BIT_FLAG_MAP : str
    Mapping from healsparse pixel bit value to output flag value, formatted as
    ``"<bit>:<flag>, <bit>:<flag>, ..."`` (e.g. ``"64:1, 2048:2"``). A pixel
    carrying several bits receives the bitwise-OR of the mapped flag values.
OFF_MAP_FLAG : int, optional
    Flag value assigned to pixels that fall outside the healsparse map
    footprint (sentinel pixels); default is ``0``. Off-footprint pixels are
    usually unobserved, so a non-zero value flags them.
USE_EXT_FLAG : bool, optional
    If ``True``, sum an external instrument flag file (given as the second
    input) into the rasterized mask; default is ``False``. Needed for
    exposures, where saturation / bleeding / bad columns arrive with the data.
HDU : int, optional
    HDU of the external instrument flag FITS file; default is ``0``
PREFIX : str, optional
    Prefix prepended to the output file name base ``flag``; default is ``""``
CHUNK_SIZE : int, optional
    Number of image rows evaluated per chunk; default is chosen from the image
    size to bound memory (``1e8`` px for tiles, ``1e7`` px for exposure CCDs)

"""

__all__ = ["mask_ext"]
