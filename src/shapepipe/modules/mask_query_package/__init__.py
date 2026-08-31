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
    ``0`` for an object no configured map flags, nonzero otherwise. What the
    nonzero value *is* depends on the maps: a boolean map — which is what the
    UNIONS per-bit products are, and what the shipped config names — can only
    contribute ``1``, so with those the column is 0/1 and says nothing about
    which map fired. An integer map contributes its own value (optionally
    ``& MASK_BITS``), and contributions are OR-ed, so a bit-packed map does
    carry its bits through. Nothing downstream reads more than ``== 0``.

The single column exists because ``setools`` expressions support only
``< > <= >= == !=`` — no bitwise operators — so any bit selection has to happen
here, leaving the config a plain test for zero.

Nothing cuts on it by default
=============================

The shipped ``star_selection.setools`` does NOT cut on ``FLAG_EXT``. The PSF
star selection is permissive: it rejects on the instrument flags
(``IMAFLAGS_ISO == 0``) and leans on outlier rejection for the rest. The
column is written for transparency and measurement — so the effect of the
external masks on the star sample can be *measured* before it is imposed —
which is the same principle as ``make_cat``'s unfiltered ``MASK_<band>``
columns: flag here, cut downstream.

Feeding the external masks to the star selection is a one-line change if
outlier rejection turns out not to be robust enough: add ``FLAG_EXT == 0`` by
each
``IMAFLAGS_ISO == 0`` in ``star_selection.setools``. The escape hatch is
deliberate, and the file's header says so.

The lookup itself lives in :mod:`shapepipe.utilities.mask_query`, shared with
``make_cat``'s per-band ``MASK_<band>`` columns, so the healsparse primitive is
written once. That module's docstring documents the off-coverage convention.

What gets queried is deliberately narrow
========================================

``MASK_PATHS`` is a *list of maps to record against each PSF-star candidate*,
not a list of every mask that exists. The committed configs name exactly one
map — the UNIONS star-body product (bit 2). Halo bits 0 and 1 are excluded on
purpose: halos flag objects for the final catalogue, they say nothing about
whether a star is a good PSF sample (mask-force telecon, 2026-07-21). MaxiMask
is not queried here either.

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
