"""PSFEX PACKAGE.

This package contains the module for ``psfex``.

:Authors: Axel Guinot, Tobias Liaudat

:Parent modules:
- ``sextractor_runner``
- ``setools_runner``

:Input: SExtractor catalogue

:Output: PSFEx PSF model

Description
===========

This module calls the community standard point spread function (PSF) modelling
software PSFEx :cite:`bertin:11` to build a PSF model.

This module requires a previous run of SExtractor.

Module-specific config file entries
===================================

EXEC_PATH : str
    full path to the PSFEx executable on the system; if not set the
    version controlled PSFEx installation in the ShapePipe environment
    will be used
DOT_PSFEX_FILE : str
     full path to the ``.psfex`` configuration file for PSFEx
CHECKIMAGE : list, optional
     list of header key names corresponding to check images

"""

__all__ = ['psfex_script']
