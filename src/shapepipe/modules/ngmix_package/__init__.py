"""NGMIX PACKAGE.

This package contains the module for ``ngmix``.

:Author: Lucie Baumont, Axel Guinot

:Parent modules:

- ``sextractor_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
- ``vignetmaker_runner``
- ``merge_headers_runner``

:Input: Galaxy image vignets

:Output: Shape catalogue

Description
===========

This module calls the weak-lensing shape measurment and metacalibration
software NGMIX :cite:`sheldon:15`. In particular the GalSim :cite:`rowe:15`
fitting method of NGXMIX is called to obtain galaxy shape measurements.
The metacalibration routines of NGMIX are also called to provide all of the
measurements required to calibrate the shear values.

The merged single-exposure WCS header log (``log_exp_headers.sqlite`` from
``merge_headers_runner``) is passed as the last entry of
``FILE_PATTERN``/``FILE_EXT``, not via a config option.

Module-specific config file entries
===================================

MAG_ZP : float
    Photometric zero point
PIXEL_SCALE : float
    Pixel scale in arcseconds
SAVE_BATCH : int, optional
    Save the output catalogue in batches of this size; default is ``-1``
    (no batch saving)
ID_OBJ_MIN : int
    ID of first galaxy object to be processed; not used if set to ``-1``
    (default). Environment variables are expanded, so an orchestrator can
    set the object range per chunk, for example
    ``ID_OBJ_MIN = $SP_NGMIX_ID_OBJ_MIN``.
ID_OBJ_MAX : int
    ID of last galaxy object to be processed; not used if set to ``-1``
    (default). Environment variables are expanded, as for ``ID_OBJ_MIN``.
BKG_RMS_VIGNET_PATH : str, optional
    Path to a ``background_rms_vignet*.sqlite`` file produced by
    ``vignetmaker_runner``. The string may contain
    ``{file_number_string}``, which is replaced by the current tile ID.

Random number generation
========================

Each object gets its own random number stream, seeded from its sky position
and CCD (``position_seed``, ngmix#796). The output is therefore identical
whether a tile is processed in one go or split into object chunks with
``ID_OBJ_MIN``/``ID_OBJ_MAX``. The older tile-seeded mode is retired; the
config option ``SEED_FROM_POSITION`` is obsolete. Setting it to ``False``
raises an error, so a stale config cannot silently change the RNG.

"""

__all__ = ["ngmix"]
