#!/usr/bin/env python3
"""Print one ngmix chunk's closed object-ID range as bash exports.

Run inside the chunk's shell, from the tile's own sexcat, because the range is
only knowable at execution time (PRD D4)::

    eval "$(ngmix_range.py --run-dir $SP_RUN --chunk 3 --n-chunks 8)"
    # -> export NGMIX_ID_MIN=751; export NGMIX_ID_MAX=1125

SExtractor's NUMBER column (ngmix's obj_id) runs 1..N contiguous, so covering
[1, N] processes every object exactly once. The first N-1 chunks get equal
shares; the last takes the remainder, with a CLOSED upper bound at n_obj —
never ID_OBJ_MAX = -1, which ngmix reads as unbounded and would double-count.
Each tile is sized from its OWN count (the bash monolith used the campaign
average and left the last chunk open).
"""

import argparse
from pathlib import Path


def id_ranges(n_obj: int, n_chunks: int) -> list[tuple[int, int]]:
    base, rem = divmod(n_obj, n_chunks)
    ranges, lo = [], 1
    for k in range(1, n_chunks + 1):
        hi = lo + base + (rem if k == n_chunks else 0) - 1
        ranges.append((lo, hi))
        lo = hi + 1
    return ranges


def object_count(run_dir: Path) -> int:
    """NAXIS2 of the last HDU of this tile's sexcat (get_number_objects.py)."""
    from astropy.io import fits

    cats = sorted((run_dir / "output" / "run_sp_tile_Sx").glob(
        "sextractor_runner/output/sexcat*.fits"))
    if not cats:
        raise SystemExit(f"[ngmix_range] FATAL: no sexcat under {run_dir}")
    with fits.open(cats[0]) as hdul:
        return int(hdul[-1].header["NAXIS2"])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--chunk", required=True, type=int)
    p.add_argument("--n-chunks", required=True, type=int)
    a = p.parse_args()
    lo, hi = id_ranges(object_count(a.run_dir), a.n_chunks)[a.chunk - 1]
    print(f"export NGMIX_ID_MIN={lo}; export NGMIX_ID_MAX={hi}")


if __name__ == "__main__":
    main()
