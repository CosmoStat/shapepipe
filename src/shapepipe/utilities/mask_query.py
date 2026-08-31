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
* ``mask_query`` writes a single integer ``FLAG_EXT`` column onto the exposure
  SExtractor catalogue, combining the configured maps into "clean (0) or
  flagged (nonzero)" so that ``setools`` — whose expression language has no
  bitwise operators — can cut on ``FLAG_EXT == 0``.

Coverage
--------
``healsparse.HealSparseMap.get_values_pos`` returns a map's *sentinel* for
positions outside its coverage: ``False`` for boolean maps, and typically
``-1`` for integer maps. ``make_cat`` passes that sentinel through verbatim,
which is the documented off-map flag for the final catalogue.

``flag_positions`` instead treats off-coverage as **not flagged**, for both map
kinds. This makes the integer case agree with the boolean case (whose sentinel
is literally ``False``) rather than diverge from it, and it keeps a map whose
coverage does not reach an exposure from silently rejecting every star on it.
Off-coverage counts are logged so the situation is visible rather than silent.

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


def query_map(path, ra, dec):
    """Query Map.

    Read a healsparse map and return its value at each world position.

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
    import healsparse

    mask_map = healsparse.HealSparseMap.read(path)

    return np.asarray(mask_map.get_values_pos(ra, dec, lonlat=True))


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

    for path in paths:
        values = query_map(path, ra, dec)

        if values.dtype == bool:
            contribution = values.astype(np.int64)
            n_off = 0
        else:
            integer = values.astype(np.int64)
            # Off-coverage: the sentinel, negative by healsparse convention.
            # Zeroed rather than OR-ed in, see this module's docstring.
            off_coverage = integer < 0
            n_off = int(np.count_nonzero(off_coverage))
            integer = np.where(off_coverage, 0, integer)
            if bits is not None:
                integer &= bits
            contribution = integer

        flag |= contribution

        if w_log is not None:
            w_log.info(
                f"Mask query {path}: "
                f"{int(np.count_nonzero(contribution))}/{ra.size} objects "
                f"flagged, {n_off} outside coverage"
            )

    return flag
