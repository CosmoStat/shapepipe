"""MERGE STAR CATALOGUES PACKAGE.

This package contains the module for ``merge_starcat``.

:Author: Tobias Liaudat, Martin Kilbinger

:Parent module: ``mccd_fit_val_runner``

:Input: PSFEx or MCCD star catalogues

:Output: merged star catalogue

Description
===========

This module merged the star catalogues used for building either PSFEx or MCCD
PSF models.

Module-specific config file entries
===================================

PSF_MODEL : str
    PSF model used; options are ``psfex`` or ``mccd``

"""

__all__ = ['merge_starcat']
