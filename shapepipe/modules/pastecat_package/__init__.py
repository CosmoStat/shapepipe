"""PASTE CATALOGUES PACKAGE.

This package contains the module for ``pastecat``.

:Authors: Martin Kilbinger, Axel Guinot

:Parent module: ``sextractor_runner``

:Input: SExtractor catalogues

:Output: Pasted catalogue

Description
===========

This module pastes (or concatenates) a series of SExtractor catalogues into
as single output catalogue.

Module-specific config file entries
===================================

CHECK_COL_NAME : str, optional
    SExtractor column name to use to ensure consistency between the number of
    rows for the catalogues to be pasted
HDU : list, optional
    HDU numbers of the input catalogues; default is set to ``2`` for all
    catalogues
PREFIX : str, optional
    Ouput file name prefix; default is ``cat_pasted``
EXT_NAME : list, optional
    List of HDU extension names; default is to use the input file names

"""

__all__ = ["pastecat"]
