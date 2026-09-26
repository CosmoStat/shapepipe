"""COVERAGE_MAP_BUILDER

Stamp per-CCD sky footprints into a HealSparse coverage (nexp) map.

Each polygon is one CCD footprint. Since the CCDs of a single exposure do not
overlap, stamping value 1 per CCD polygon and accumulating makes the map count,
per sky pixel, the number of exposures with a valid PSF model covering that
pixel.

Arrays in, map out. The only caller is the ``coverage_map`` rule
(``workflow/scripts/coverage_map.py``), which reads the corners out of the
per-exposure ``exp_footprint.json`` records; the nine-column ``exp_ra_dec.txt``
the pre-Snakemake chain appended to is gone, and so is the CLI that parsed it.
What is left here is survey geometry — the RA-seam guard, the pole guard, the
nside validation and the median smoothing — none of which depends on how the
corners arrived.

Author: Mike Hudson, Martin Kilbinger <martin.kilbinger@cea.fr>

"""

import healsparse as hsp
import hpgeom as hpg
import numpy as np

# Declination beyond which a polygon is considered too close to a pole for
# HealSparse's planar polygon fill to be reliable. CCD footprints are ~10
# arcmin, so this is a guard against pathological headers, not a common path.
_DEC_POLE_LIMIT = 89.0


def unwrap_ra(ra):
    """Unwrap RA.

    Shift RA values onto a common branch when a polygon straddles the RA=0
    seam, so the corners describe the small footprint rather than its ~360°
    complement. Values above 180 deg are shifted by -360 deg when the raw
    spread exceeds 180 deg; otherwise the values are returned unchanged.

    Parameters
    ----------
    ra : array_like
        polygon RA corners, in degrees

    Returns
    -------
    numpy.ndarray
        RA corners on a common branch, in degrees

    """
    ra = np.asarray(ra, dtype=float)
    if ra.max() - ra.min() > 180.0:
        return np.where(ra > 180.0, ra - 360.0, ra)
    return ra


def check_nside(nside_coverage, nside):
    """Raise unless both HealSparse resolutions are positive powers of two.

    Checked before any polygon is stamped: healsparse itself complains only
    once a map is being made, by which time the caller has already paid for
    reading every footprint record on the products root.

    Parameters
    ----------
    nside_coverage : int
        HealSparse coverage nside
    nside : int
        HealSparse map nside

    Raises
    ------
    ValueError
        if either value is not a positive power of two

    """
    for name, n in (("nside_coverage", nside_coverage), ("nside", nside)):
        if n <= 0 or n & (n - 1):
            raise ValueError(f"{name} must be a power of 2, got {n}")


def build_map(
    ccd_ids, ra, dec, nside_coverage, nside, verbose=False
):
    """Stamp one polygon per CCD footprint into a HealSparse nexp map.

    Since the CCDs of a single exposure do not overlap, accumulating value 1
    per CCD polygon makes the map count, per sky pixel, the number of
    exposures with a valid PSF model covering it.

    Parameters
    ----------
    ccd_ids : array_like
        length-N CCD IDs, used only in the skipped-polygon warnings
    ra : array_like
        ``(N, 4)`` polygon RA corners, in degrees
    dec : array_like
        ``(N, 4)`` polygon Dec corners, in degrees
    nside_coverage : int
        HealSparse coverage nside
    nside : int
        HealSparse map nside
    verbose : bool, optional
        print progress

    Returns
    -------
    healsparse.HealSparseMap
        the nexp map

    """
    check_nside(nside_coverage, nside)

    ra = np.atleast_2d(np.asarray(ra, dtype=float))
    dec = np.atleast_2d(np.asarray(dec, dtype=float))

    if verbose:
        print(
            f"Creating HealSparse map (nside_coverage={nside_coverage},"
            f" nside={nside})"
        )
    m = hsp.HealSparseMap.make_empty(nside_coverage, nside, np.uint16)

    if verbose:
        print("Adding polygons to map")

    n_added = 0
    n_skipped = 0
    for i in range(len(ccd_ids)):
        dec_i = dec[i]

        # Pole guard: HealSparse's planar polygon fill degrades near the
        # poles. CCD footprints never reach here, so warn and skip.
        if np.any(np.abs(dec_i) >= _DEC_POLE_LIMIT):
            print(
                f"Warning: skipping CCD {ccd_ids[i]} with |dec| >= "
                f"{_DEC_POLE_LIMIT} (too close to a pole)"
            )
            n_skipped += 1
            continue

        # RA-wrap guard: put corners on a common branch across the seam.
        ra_i = unwrap_ra(ra[i])

        m += hsp.Polygon(ra=list(ra_i), dec=list(dec_i), value=1)
        n_added += 1

        if verbose and i % 1000 == 0:
            print(f"{i:6d} / {len(ccd_ids):6d}")

    print(f"Added {n_added} polygons to map")
    if n_skipped > 0:
        print(f"Skipped {n_skipped} polygons near a pole")

    return m


def median_filter(hsp_map, n_iterations=1):
    """Smooth an nexp map: each pixel becomes its neighbourhood median.

    The neighbourhood is a pixel and its eight HEALPix neighbours; out-of-map
    neighbours read as the map's sentinel and so pull an edge pixel down,
    which is the intended behaviour for a coverage map (a lone pixel is noise,
    not depth). The ``coverage_map`` rule does NOT smooth — the map it writes
    is the raw exposure count, which is what sp_validation's mask application
    wants — so this is the offline step, for looking at a map rather than
    applying one.

    Parameters
    ----------
    hsp_map : healsparse.HealSparseMap
        input nexp map
    n_iterations : int, optional
        number of smoothing passes

    Returns
    -------
    healsparse.HealSparseMap
        filtered map

    """
    nside = hsp_map.nside_sparse

    for _ in range(n_iterations):
        new_hsp = hsp_map.copy()
        pixs = hsp_map.valid_pixels
        n = hpg.neighbors(nside, pixs)
        n = np.where(n >= 0, n, 0)
        pix2 = pixs.reshape(len(pixs), 1)
        n = np.hstack([n, pix2])
        mn = hsp_map[n]
        new = np.median(mn, axis=1).astype(np.uint16)
        new_hsp[pixs] = new
        hsp_map = new_hsp

    return hsp_map
