"""FAKE PSF PACKAGE.

This package creates fake PSF postage-stamp dictionaries for image simulations.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Parent module: ``sextractor_runner``

:Input: Tile SExtractor catalogue (multi-epoch FITS format)

:Output: Per-galaxy PSF sqlite dictionary (``galaxy_psf-XXX-XXX.sqlite``)

Description
===========

For each galaxy in the tile catalogue the module looks up the pre-computed
PSF stamp for every contributing exposure-CCD combination from a pickled
PSF dictionary and writes the result to a SqliteDict file in the same format
produced by ``psfex_interp_runner``.

Module-specific config file entries
====================================

PSF_DICT_PATH : str
    Path to the pickled PSF dictionary (``Full_psf_dict.pickle``).
"""

__all__ = ["fake_psf.py"]
