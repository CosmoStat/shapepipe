"""MATCH EXTERNAL PACKAGE.

This package contains the module for ``match_external``.

:Authors: Martin Kilbinger, Xavier Jimenez

:Parent module: ``sextractor_runner``

:Input: SExtractor catalogue

:Output: Matched catalogue

Description
===========

This module matches an external catalogue to a ShapePipe (SExtractor)
catalogue. Matching is performed by checking if the 2D angular separation
between objects is less than a specified tolerance.

Module-specific config file entries
===================================

TOLERANCE : float
    Tolerance of matching distance in arcseconds
COL_MATCH : list
    Internal SExtractor column names to be matched
HDU : int, optional
    Internal SExtractor catalogue HDU number; default is ``2``
MODE : str
    Run mode for module, options are ``CLASSIC`` or ``MULTI-EPOCH``
EXTERNAL_CAT_PATH : str
    Path to the external catalogue
EXTERNAL_COL_MATCH : list
    External SExtractor column names to be matched
EXTERNAL_COL_COPY : list
    External SExtractor column names to be copied to the matched catalogue
EXTERNAL_HDU : int, optional
    External SExtractor catalogue HDU number; default is ``1``
PREFIX : str, optional
    Output file name prefix; default is ``cat_matched``
MARK_NON_MATCHED : float, optional
    If set all objects will be included in the output catalogue and unmatched
    objects will be marked with the value specified
OUTPUT_DISTANCE : bool, optional
    Option to output the matching distance; default is ``False``

"""

__all__ = ["match_external"]
