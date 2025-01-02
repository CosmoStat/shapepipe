"""MERGE STAR CATALOGUES PACKAGE.

This package contains the module for ``merge_starcat``.

:Author: Tobias Liaudat, Martin Kilbinger

:Parent module: ``mccd_fit_val_runner``, ``mccd_interp_runner``,
  ``mccd_val_runner`` or ``psfex_interp_runner``

:Input: PSFEx or MCCD star catalogues

:Output: Merged star catalogue

Description
===========

This module merges the star catalogues used for building either PSFEx or MCCD
PSF models.

Module-specific config file entries
===================================

PSF_MODEL : str
    PSF model used; options are ``psfex`` or ``mccd``, ``setools''
HDU : int, optional
    HDU number of input catalogue table, default is ``1`` for ``mccd``
    and ``2`` for ``psfex`` and ``setools``

"""

__all__ = ["merge_starcat"]
