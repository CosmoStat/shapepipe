"""MASK QUERY PACKAGE.

This package contains the module for ``mask_query``.

:Author: Claude Fable 5

:Parent module: ``sextractor_runner``

:Input: Single-exposure single-CCD SExtractor catalogue

:Output: The same catalogue with an added ``MASK_EXT`` column

Two classes of mask
===================

ShapePipe distinguishes them, and they are not interchangeable:

**Instrument flags** mark a CORRUPTED MEASUREMENT — bad columns, saturation —
in the flag image delivered with each exposure. A flagged pixel carries no
usable signal, so these are the only masks that reject anything inside the
pipeline: ``setools`` drops flagged stars via ``IMAFLAGS_ISO``, and ``ngmix``
zero-weights flagged pixels and drops epochs that lose too many.

**Healsparse masks** are sky-fixed LOCATION flags — a star halo, a manual
region, a band with no data. They say where an object is, not that its pixels
are broken, so what to do about one is an analysis decision. They are queried
into columns: ``MASK_EXT`` here, ``MASK_<band>`` in ``make_cat``. The one
in-pipeline cut on them is PSF-star selection, which rejects candidates with
``MASK_EXT != 0``; the ``MASK_<band>`` columns are never cut on.

Description
===========

This module sits between SExtractor and ``setools`` on the exposure chain: it
reads each detection's windowed world position (``XWIN_WORLD``, ``YWIN_WORLD``
in the ``LDAC_OBJECTS`` extension) and looks it up in the configured healsparse
maps, writing one integer column:

``MASK_EXT``
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

A no-op without maps
====================

With no ``MASK_PATHS`` the module is a strict no-op: the catalogue is copied
through unchanged, with no ``MASK_EXT`` column at all — the same gating
``make_cat`` gives ``MASK_EXT_PATHS``. It stays in the ``MODULE`` chain either
way. The shipped data configs name the UNIONS manual mask (bit 8), and
``star_selection.setools`` rejects ``MASK_EXT != 0`` in its ``preselect``,
``flag`` and ``star_selection`` blocks; a config without the map must drop those
lines too, because ``setools`` cannot cut on a missing column.

The lookup itself lives in :mod:`shapepipe.utilities.mask_query`, shared with
``make_cat``'s per-band ``MASK_<band>`` columns, so the healsparse primitive is
written once. That module's docstring documents the off-coverage convention.

What gets queried is deliberately narrow
========================================

``MASK_PATHS`` is a *list of maps that veto PSF-star candidates*, not a list
of every mask that exists. The committed configs name exactly one map — the
UNIONS manual mask (bit 8), which covers large galaxies, globular clusters and
other extended objects. Halo bits 0 and 1 are excluded on purpose: halos flag
objects for the final catalogue, they say nothing about whether a star is a
good PSF sample (mask-force telecon, 2026-07-21). The star-body map (bit 2) is
excluded because it flags the stars themselves. MaxiMask is not queried here
either.

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
