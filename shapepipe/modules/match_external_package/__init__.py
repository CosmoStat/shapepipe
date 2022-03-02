"""MATCH EXTERNAL PACKAGE.

This package contains the module for ``match_external``.

:Author: Martin Kilbinger, Xavier Jimenez

:Parent module: ``sextractor_runner``

:Input: SExtractor catalogue

:Output: matched catalogue

Description
===========

This module matches an external catalogue to a ShapePipe (SExtractor)
catalogue. Matching is performed by checking if the 2D angular separation
between objects is less than a specified tolerance.

Module-specific config file entries
===================================

TOLERANCE : float
    tolerance of matching distance in arcseconds
COL_MATCH : list
    internal SExtractor column names to be matched
HDU : int, optional
    internal SExtractor catalogue HDU number; default is ``2``
MODE : str
    run mode for module, options are ``CLASSIC`` or ``MULTI-EPOCH``
EXTERNAL_CAT_PATH : str
    path to the external catalogue
EXTERNAL_COL_MATCH : list
    external SExtractor column names to be matched
EXTERNAL_COL_COPY : list
    external SExtractor column names to be copied to the matched catalogue
EXTERNAL_HDU : int, optional
    external SExtractor catalogue HDU number; default is ``1``
PREFIX : str, optional
    output file name prefix; default is ``cat_matched``
MARK_NON_MATCHED : float, optional
    if set all objects will be included in the output catalogue and unmatched
    objects will be marked with the value specified
OUTPUT_DISTANCE : bool, optional
    option to output the matching distance; default is ``False``

"""

__all__ = ['match_external']
