"""MASK QUERY PACKAGE.

This package contains the module for ``mask_query``.

:Author: Claude Fable 5, for PR #847

:Parent module: ``sextractor_runner``

:Input: Single-exposure single-CCD SExtractor catalogue

:Output: The same catalogue with an added ``FLAG_EXT`` column

Description
===========

ShapePipe consumes sky-fixed masks by querying them, not by rasterizing them.
This module sits between SExtractor and ``setools`` on the exposure chain: it
reads each detection's windowed world position (``XWIN_WORLD``, ``YWIN_WORLD``
in the ``LDAC_OBJECTS`` extension) and looks it up in the configured healsparse
maps, writing one integer column:

``FLAG_EXT``
    ``0`` for an object no configured map flags, nonzero otherwise. The nonzero
    value is the bitwise OR of the contributing map values, so it says *which*
    bits fired, but nothing downstream is required to read it that way.

The single column exists because ``setools`` expressions support only
``< > <= >= == !=`` — no bitwise operators — so the bit selection has to
happen here. ``star_selection.setools`` cuts on ``FLAG_EXT == 0`` beside its
existing ``IMAFLAGS_ISO == 0``: instrument flags reject pixels, the queried
masks reject objects.

The lookup itself lives in :mod:`shapepipe.utilities.mask_query`, shared with
``make_cat``'s per-band ``MASK_<band>`` columns, so the healsparse primitive is
written once. That module's docstring documents the off-coverage convention.

The diet is deliberately narrow
===============================

``MASK_PATHS`` is a *list of maps to reject PSF stars on*, not a list of every
mask that exists. The committed configs name exactly one map — the UNIONS
star-body product (bit 2) — beside the instrument flags SExtractor already
reads. Halo bits 0 and 1 are excluded on purpose: halos flag objects for the
final catalogue, they do not reject PSF stars (mask-force telecon, 2026-07-21).
MaxiMask is not in the diet either.

Widening it costs a config edit and no code — add a path. That is why the
contract is a path list rather than a bit mask: the UNIONS products are one
boolean map per bit, so choosing bits *is* choosing files, and ``MASK_BITS``
exists only for integer maps that pack several bits into one file.

Module-specific config file entries
===================================

MASK_PATHS : str
    Comma-separated healsparse map paths to query
MASK_BITS : int, optional
    Bit mask applied to integer maps (``value & MASK_BITS`` flags); default is
    to flag on any nonzero value. Ignored for boolean maps, which flag on
    ``True``
PREFIX : str, optional
    Output file prefix

"""

__all__ = ["mask_query"]
