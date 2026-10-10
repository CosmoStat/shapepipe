r"""MAKE CATALOGUE PACKAGE.

This package contains the module for ``make_cat``.

:Author: Axel Guinot

:Parent modules:

- ``sextractor_runner``
- ``psfex_interp_runner`` or ``mccd_interp_runner``
- ``ngmix_runner``

:Input: SExtractor catalogues, sqlite catalogue

:Output: SExtractor catalogue

Description
===========

This module creates a *final* catalogue combining the output of various
previous module runs. This gathers all relevant information on the measured
galaxies for weak-lensing post-processing. This includes galaxy detection and
basic measurement parameters, the PSF model at galaxy positions, and the
shape measurement. ``N_EPOCH`` counts the object's epochs with a validated
PSF model, the epochs ngmix can fit; ``N_EPOCH_OVERLAP`` counts the exposure
CCDs whose footprint holds the object, including CCDs whose PSF model failed
validation. Every detected object is kept: the catalogue carries no
star/galaxy classification, which is done downstream. Each object is
tagged with its source tile via the ``TILE_ID`` column and carries the
survey-wide object ID ``TILE_UNIQUE_ID = tile_id * 10**6 + NUMBER``
(``tile_id = RRR * 1000 + DDD``, see
:func:`shapepipe.utilities.cfis.get_tile_unique_id`); objects duplicated
across overlapping tiles are not deduplicated here, and are left to
downstream selection.

Module-specific config file entries
===================================

SHAPE_MEASUREMENT_TYPE : list
    Shape measurement method; the only valid option is ``ngmix`` (the knob is
    retained as the extension point for a future estimator family)
SAVE_PSF_DATA : bool, optional
    Save PSF information if ``True``; default value is ``False``
N_EPOCH_SLOTS : int, optional
    Number of slots written for each per-epoch PSF column family
    (``HSM_*_PSF_n``, ``EXP_ID_n``, ``CCD_n``) when ``SAVE_PSF_DATA`` is
    ``True``; unfilled slots hold the family's sentinel, and an object with
    more epochs than slots raises an error. A fixed count gives every tile
    the same schema. Default is the tile's maximum ``N_EPOCH_OVERLAP`` plus
    one

"""

__all__ = ["make_cat"]
