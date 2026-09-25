"""MASK QUERY.

The one healsparse lookup in ShapePipe.

ShapePipe does not generate or rasterize masks. Sky-fixed masks are supplied
as healsparse maps and are consumed by *querying them at object positions*:
every object gets its mask value(s) as catalogue columns, and rejection happens
at the catalogue level, never at the pixel level. The only mask that still
reaches pixels is the per-exposure instrument flag image delivered with it.

Two callers share the primitive defined here:

* ``make_cat`` writes one ``MASK_<band>`` column per configured map, carrying
  the map value verbatim (no interpretation, no filtering);
* ``mask_query`` writes a single integer ``MASK_EXT`` column onto the exposure
  SExtractor catalogue, combining the configured maps into "clean (0) or
  flagged (nonzero)" so that ``setools`` — whose expression language has no
  bitwise operators — *could* cut on ``MASK_EXT == 0``. The shipped selection
  does not; the column is carried for measurement (see that module).

Partial reads
-------------
Never read a whole map. The UNIONS products are ~550 MB each and ``mask_query``
runs PER CCD, so a full read would cost ~22 GB of I/O per 40-CCD exposure to
answer questions about a 0.06 deg² footprint. Instead the coverage index is
read first (``HealSparseCoverage``, a few hundred kB), the coverage pixels the
catalogue actually touches are computed with ``hpgeom.angle_to_pixel`` at the
map's ``nside_coverage``, and only those get loaded.

Every queried position falls inside exactly one coverage pixel and all of them
are requested, so the padding by ``hpgeom.neighbors`` is insurance rather
than necessity — at ``nside_coverage=128`` a coverage pixel of the star map
is ~23 kB, so the padding costs under a megabyte and buys immunity to any
edge convention we did not think of. ``test_partial_read_matches_full`` is
what actually holds the two paths equal.

Measured on the DR6 star-body map — the bit-2 rung of the ladder the configs
name ``mask_ugriz_nside131072_n4.hsp``, run against the staged single-band copy
``mask_r_nside131072_n4.hsp``, 583 MB, ``nside_coverage=128`` — for 2000
positions in one CCD-sized box, one process each: partial 0.102 s / 157 MiB
peak RSS, full 12.8 s / 3364 MiB, identical values. Per 40-CCD exposure that
is ~4 s against ~8.5 min of map reading, and the query adds ~0.15 GB to a rule
already asking for 16 GB.

Coverage
--------
``healsparse.HealSparseMap.get_values_pos`` returns a map's *sentinel* for
positions outside its coverage: ``False`` for boolean maps, and typically
``-1`` for integer maps. ``make_cat`` passes that sentinel through verbatim,
which is the documented off-map flag for the final catalogue.

Coverage is reported from the COVERAGE MASK, not ``valid_mask=True``. For a
boolean map — which is what the UNIONS per-bit products are — healsparse
stores only the ``True`` pixels, so ``valid_mask`` returns the value itself
and cannot tell "inside the footprint and clean" from "outside it entirely".
The coverage mask can, at ``nside_coverage`` resolution, which is the scale the
question is asked at anyway: does this map reach this exposure at all?

``flag_positions`` treats off-coverage as **not flagged**, for both map kinds.
This makes the integer case agree with the boolean case (whose sentinel is
literally ``False``) rather than diverge from it, and it keeps a map whose
coverage does not reach an exposure from silently rejecting every star on it.
That is also why the all-off-coverage case is logged as a WARNING: it is
indistinguishable from "nothing is masked here" in the output column, so it has
to be distinguishable in the log.

:Author: Claude Fable 5, for PR #847

"""

import numpy as np


def parse_map_paths(paths_str):
    """Parse Map Paths.

    Parse a comma-separated list of healsparse map paths.

    Parameters
    ----------
    paths_str : str
        Comma-separated map paths, e.g. ``/a/star.hsp, /b/maximask.hsp``

    Returns
    -------
    list
        Map paths, stripped of surrounding whitespace, empty entries dropped

    """
    return [path.strip() for path in paths_str.split(",") if path.strip()]


def _covering_pixels(coverage, ra, dec):
    """Covering Pixels.

    The map's coverage pixels touched by these positions, padded by their
    neighbours and intersected with what the map actually holds.

    Parameters
    ----------
    coverage : healsparse.HealSparseCoverage
        Coverage index of the map
    ra : numpy.ndarray
        Right ascension in degrees
    dec : numpy.ndarray
        Declination in degrees

    Returns
    -------
    numpy.ndarray
        Coverage pixel indices to load, possibly empty

    """
    import hpgeom

    nside_coverage = coverage.nside_coverage
    touched = np.unique(
        hpgeom.angle_to_pixel(nside_coverage, ra, dec, nest=True)
    )
    padded = np.unique(
        np.concatenate(
            [touched, hpgeom.neighbors(nside_coverage, touched).ravel()]
        )
    )
    # neighbors() returns -1 for a non-existent neighbour.
    padded = padded[padded >= 0]

    return padded[coverage.coverage_mask[padded]]


def query_map_coverage(path, ra, dec):
    """Query Map With Coverage.

    Read only the part of a healsparse map these positions need, and return
    both its value at each position and whether each position is inside the
    map's coverage.

    Parameters
    ----------
    path : str
        Path to the healsparse map
    ra : numpy.ndarray
        Right ascension in degrees
    dec : numpy.ndarray
        Declination in degrees

    Returns
    -------
    tuple
        ``(values, in_coverage)`` — the map value at each position (the map's
        sentinel outside coverage) and a boolean array, both of length
        ``len(ra)``

    """
    import healsparse
    import hpgeom

    ra = np.asarray(ra)
    dec = np.asarray(dec)

    coverage = healsparse.HealSparseCoverage.read(path)
    nside_coverage = coverage.nside_coverage

    in_coverage = coverage.coverage_mask[
        hpgeom.angle_to_pixel(nside_coverage, ra, dec, nest=True)
    ]

    pixels = _covering_pixels(coverage, ra, dec)

    if pixels.size == 0:
        # healsparse raises when no requested pixel is in the coverage map, so
        # the empty case is answered without asking it: load one arbitrary
        # coverage pixel purely to learn the dtype and sentinel, and return
        # that sentinel everywhere. Same answer, one small read.
        covered = np.flatnonzero(coverage.coverage_mask)
        if covered.size == 0:
            raise ValueError(f"healsparse map {path} has empty coverage")
        mask_map = healsparse.HealSparseMap.read(
            path, pixels=[int(covered[0])]
        )
        values = np.full(ra.size, mask_map.sentinel, dtype=mask_map.dtype)
        return values, in_coverage

    mask_map = healsparse.HealSparseMap.read(
        path, pixels=[int(pixel) for pixel in pixels]
    )
    values = np.asarray(mask_map.get_values_pos(ra, dec, lonlat=True))

    return values, in_coverage


def query_map(path, ra, dec):
    """Query Map.

    Return a healsparse map's value at each world position, reading only the
    coverage pixels those positions touch.

    Parameters
    ----------
    path : str
        Path to the healsparse map
    ra : numpy.ndarray
        Right ascension in degrees
    dec : numpy.ndarray
        Declination in degrees

    Returns
    -------
    numpy.ndarray
        Map value at each position; positions outside the map's coverage carry
        the map's sentinel value

    """
    values, _ = query_map_coverage(path, ra, dec)

    return values


def flag_positions(paths, ra, dec, bits=None, w_log=None):
    """Flag Positions.

    Combine one or more healsparse masks into a single per-object integer flag.

    Each map contributes at each position:

    * boolean map: ``1`` where the map is ``True``, ``0`` elsewhere;
    * integer map: the map value, optionally restricted to ``bits``
      (``value & bits``); ``0`` where the value is zero or the position is
      outside coverage.

    Contributions are combined with a bitwise OR, so the returned flag is zero
    for a clean object and carries the union of the bits that fired otherwise.

    Parameters
    ----------
    paths : list
        Paths to the healsparse maps to query
    ra : numpy.ndarray
        Right ascension in degrees
    dec : numpy.ndarray
        Declination in degrees
    bits : int, optional
        Bit mask applied to integer maps; default ``None`` means any nonzero
        value flags
    w_log : logging.Logger, optional
        Logging instance

    Returns
    -------
    numpy.ndarray
        Integer flag per object, ``0`` for a clean object

    """
    ra = np.asarray(ra)
    dec = np.asarray(dec)
    flag = np.zeros(ra.size, dtype=np.int64)

    if ra.size == 0:
        return flag

    for path in paths:
        values, in_coverage = query_map_coverage(path, ra, dec)

        if values.dtype == bool:
            contribution = values.astype(np.int64)
        else:
            integer = values.astype(np.int64)
            # Off-coverage: the sentinel, negative by healsparse convention.
            # Zeroed rather than OR-ed in, see this module's docstring.
            integer = np.where(integer < 0, 0, integer)
            if bits is not None:
                integer &= bits
            contribution = integer

        flag |= contribution

        n_off = int(np.count_nonzero(~in_coverage))
        if w_log is not None:
            w_log.info(
                f"Mask query {path}: "
                f"{int(np.count_nonzero(contribution))}/{ra.size} objects "
                f"flagged, {n_off} outside coverage"
            )
            if n_off == ra.size:
                # Every object reads the sentinel, so the column is all-zero
                # and looks exactly like "nothing is masked here". Say it.
                w_log.warning(
                    f"Mask query {path}: NO object is inside this map's"
                    + " coverage — the resulting MASK_EXT contribution is"
                    + " zero everywhere because the map does not reach these"
                    + " positions, not because they are clean."
                )

    return flag
