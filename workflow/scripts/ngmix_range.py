#!/usr/bin/env python3
"""Print one ngmix chunk's closed object-ID range as bash exports.

Run inside the chunk's shell, from the tile's own sexcat, because the range is
only knowable at execution time (PRD D4)::

    eval "$(ngmix_range.py --run-dir $SP_RUN --chunk 3 --n-chunks 8)"
    # -> export NGMIX_ID_MIN=751; export NGMIX_ID_MAX=1125

SExtractor's NUMBER column (ngmix's obj_id) runs 1..N contiguous, so covering
[1, N] processes every object exactly once. The bounds are CLOSED, never
ID_OBJ_MAX = -1 — ngmix treats ``id_obj_max <= 0`` as unbounded
(``ngmix_package/ngmix.py:804-806``), so an open-ended last chunk silently
re-measures the whole tile instead of its share; rule tile_ngmix carries the
straggler that taught this.

CHUNKS ARE BALANCED BY EPOCH-WEIGHTED COST, NOT BY OBJECT COUNT. ngmix fits an
object jointly across every exposure it lands on, so its cost scales with
(object, epoch) pairs, not objects: the campaign measured 0.2714 CPU-s per
pair with an intercept consistent with zero, plus ~0.05 CPU-s per-object of
setup. The old equal-count split therefore produced chunks whose cost varied
1.6x WITHIN a single tile, tracking their epoch counts at R^2 0.91-1.000. Those
chunks are siblings in the fused tile_shape group job, which ends when its
SLOWEST member does, so the spread was pure wall clock: 7.7 hours over the
34-tile smk-g4 campaign, and on the worst tile (196.307) the slowest chunk ran
26.3 min past the median and took the group to 97.4% of its wall limit.

Re-splitting all 34 of that campaign's sexcats: the old boundaries leave the
slowest chunk at 1.131x-1.627x its tile's median predicted cost (worst
200.302), the new ones at 1.0000x-1.0002x, and the sum over tiles of
slowest-chunk cost falls 132,875 -> 113,710 predicted CPU-s, 14.4%.

Moving the boundaries is scientifically free: ngmix seeds its RNG from object
POSITION (SEED_FROM_POSITION), so which chunk an object falls in cannot change
its measurement.

DETERMINISM HERE IS CORRECTNESS, NOT TIDINESS. The n_chunks chunks are separate
processes that each run this script independently and must derive the IDENTICAL
range set; if two of them disagree, objects are silently measured twice or
silently dropped and nothing downstream notices — merge_sep_cats concatenates
whatever it is given. Hence integer arithmetic end to end (weights scale to
milli-epochs so no float ever decides a boundary), and NO FALLBACK: if the
EPOCH extensions cannot be read this script exits non-zero. Such a fallback
would be the worst available behaviour precisely because it would apply only to
the chunks that hit the failure, shredding the tile's coverage instead of
failing it.
"""

import argparse
from pathlib import Path

# Weights are integers in milli-epochs (one epoch = 1000) so every boundary is
# decided by exact integer comparison, identically on all n_chunks processes.
MILLI_EPOCH = 1000

# The per-object setup cost, expressed in epochs so one integer weight carries
# both terms of the measured cost law: 0.05 CPU-s of setup / 0.2714 CPU-s per
# (object, epoch) = 0.184 epoch-equivalents. It is what keeps the zero-epoch
# objects (13 of 35,298 on tile 186.307) from weighing nothing — they still
# cost ~6% of a typical object, and a chunk handed thousands of them for free
# would be a straggler of a new kind. Only the RATIO of the two costs moves a
# boundary, which makes this robust: refitting on 186.307's eight measured
# chunk CPU times with the setup term held fixed gives 0.2619 per GEOMETRIC
# epoch (see object_epochs for why that is below 0.2714), i.e. ALPHA 0.191,
# and every boundary on that tile then moves by at most one object. Zeroing
# ALPHA entirely moves them by at most 44.
ALPHA_MILLI_EPOCHS = 184


def id_ranges(epochs, n_chunks: int) -> list[tuple[int, int]]:
    """Split object IDs ``1..len(epochs)`` into ``n_chunks`` closed ranges.

    ``epochs[i]`` is object ``i + 1``'s geometric epoch count. The ranges are
    CONTIGUOUS — ``NGMIX_ID_MIN``/``NGMIX_ID_MAX`` is an interval, not a set —
    and tile ``[1, n_obj]`` exactly, so every object is measured once.

    The objective is the SLOWEST chunk, not the average one, because the
    group job waits for it. So this minimises the maximum chunk weight
    exactly: binary-search the smallest feasible capacity, then fill left to
    right under it. Ties in that maximum break toward the earlier chunks,
    which is arbitrary but fixed, and fixed is the property that matters.

    With ``n_obj < n_chunks`` the first ``n_obj`` chunks take one object each
    and the remainder are EMPTY, written ``(n_obj + 1, n_obj)``: ``lo = hi+1``
    is the canonical empty closed interval and preserves the chain
    ``ranges[k][0] == ranges[k-1][1] + 1``. Deliberately not ``(1, 0)``, which
    is what the old equal-count split emitted here — ``ID_OBJ_MAX = 0`` is
    ngmix's unbounded sentinel, so each empty chunk would have re-measured the
    ENTIRE tile.
    """
    if n_chunks < 1:
        raise ValueError(f"n_chunks must be >= 1, got {n_chunks}")
    n_obj = len(epochs)
    if n_obj < 1:
        raise ValueError("cannot split a catalogue of zero objects")

    if n_obj < n_chunks:
        return ([(k, k) for k in range(1, n_obj + 1)]
                + [(n_obj + 1, n_obj)] * (n_chunks - n_obj))

    weights = [MILLI_EPOCH * int(e) + ALPHA_MILLI_EPOCHS for e in epochs]
    cap = _min_feasible_capacity(weights, n_chunks)

    ranges: list[tuple[int, int]] = []
    start, load = 0, 0
    for i, w in enumerate(weights):
        # Close before object i when the current chunk holds something and
        # either it cannot take i under `cap`, or the objects still to come
        # (n_obj - i) no longer outnumber the chunks still to open — the second
        # clause is what guarantees exactly n_chunks non-empty ranges when a
        # heavy head would otherwise leave the tail with nothing to hold.
        if start < i and (load + w > cap
                          or n_obj - i < n_chunks - len(ranges)):
            ranges.append((start + 1, i))
            start, load = i, 0
        load += w
    ranges.append((start + 1, n_obj))
    return ranges


def _chunks_needed(weights: list[int], cap: int) -> int:
    """Chunks a left-to-right fill uses when none may exceed ``cap``."""
    used, load = 1, 0
    for w in weights:
        if load + w > cap:
            used, load = used + 1, w
        else:
            load += w
    return used


def _min_feasible_capacity(weights: list[int], n_chunks: int) -> int:
    """Smallest ``cap`` that fits ``weights`` into ``n_chunks`` chunks.

    ``_chunks_needed`` is monotone non-increasing in ``cap``, so bisection on
    the integers ``[max(weights), sum(weights)]`` lands on the exact optimum
    without a float ever entering the comparison.
    """
    lo, hi = max(weights), sum(weights)
    while lo < hi:
        mid = (lo + hi) // 2
        if _chunks_needed(weights, mid) <= n_chunks:
            hi = mid
        else:
            lo = mid + 1
    return lo


def object_epochs(run_dir: Path):
    """Per-object geometric epoch count, from this tile's own sexcat.

    tile_detect's SExtractor post-process (``MAKE_POST_PROCESS`` in
    ``config_tile_Sx.ini``) writes one ``EPOCH_<k>`` extension per exposure
    overlapping the tile — so the extension COUNT is tile-specific and is
    discovered by name, never assumed — each with ``n_obj`` rows in NUMBER
    order and ``CCD_N < 0`` where the object misses that exposure. Summing
    ``CCD_N >= 0`` across them reproduces the final catalogue's ``N_EPOCH``
    column exactly — checked row by row against 186.307's
    ``run_sp_Mc/.../final_cat-186-307.fits``, all 35,298 of them, 7 extensions,
    116,727 pairs, mean 3.31. The post-process is upstream of the whole
    tile_shape group, so the extensions always exist by the time ngmix runs;
    their absence is a broken tile, not a case to accommodate.

    N_EPOCH is a GEOMETRIC count and mildly over-states the work, because ngmix
    drops epochs it cannot fit. The same catalogue's NGMIX_N_EPOCH is lower for
    2,203 of the 35,298 objects and never higher: 113,947 pairs against
    116,727, 2.4%. That gap is the whole of the ~3.4% by which this cost model
    over-predicts the tile's measured chunk times — refit on NGMIX_N_EPOCH the
    law is 0.2681 CPU-s per pair at R^2 0.997, against 0.2619 at R^2 0.973 on
    the geometric count. The better number is unavailable before ngmix runs and
    is not worth wanting anyway: splitting on it moves boundaries by up to 177
    objects and improves the true slowest chunk by 1.07%.

    Only the EPOCH extensions are read. ``LDAC_OBJECTS`` carries a 10 kB
    VIGNET per row (376 MB on 186.307) and pulling it in would cost more than
    the straggle this split exists to remove; the EPOCH tables are ~530 kB
    each, and its NAXIS2 comes from the header alone as the row-count check.
    """
    import numpy as np
    from astropy.io import fits

    cats = sorted((run_dir / "output" / "run_sp_tile_Sx").glob(
        "sextractor_runner/output/sexcat*.fits"))
    if not cats:
        raise SystemExit(f"[ngmix_range] FATAL: no sexcat under {run_dir}")
    with fits.open(cats[0], memmap=True) as hdul:
        epoch_hdus = [h for h in hdul if h.name.startswith("EPOCH")]
        if not epoch_hdus:
            raise SystemExit(
                f"[ngmix_range] FATAL: no EPOCH extensions in {cats[0]}; "
                "tile_detect's SExtractor post-process did not run. "
                "Refusing to fall back to an equal-object split: it would "
                "apply to only the chunks that saw this failure, and the "
                "tile's objects would be double-measured or dropped."
            )
        if "LDAC_OBJECTS" not in hdul:
            raise SystemExit(
                f"[ngmix_range] FATAL: no LDAC_OBJECTS in {cats[0]}"
            )
        n_obj = int(hdul["LDAC_OBJECTS"].header["NAXIS2"])
        expected = np.arange(1, n_obj + 1, dtype=np.int64)
        counts = np.zeros(n_obj, dtype=np.int64)
        for hdu in epoch_hdus:
            data = hdu.data
            # NUMBER is asserted, not assumed: the whole scheme is an ID
            # INTERVAL, so a permuted or gappy NUMBER column would make the
            # weights describe different objects than the bounds select.
            if len(data) != n_obj or not np.array_equal(
                np.asarray(data["NUMBER"], dtype=np.int64), expected
            ):
                raise SystemExit(
                    f"[ngmix_range] FATAL: {cats[0]}[{hdu.name}] is not "
                    f"{n_obj} rows of NUMBER = 1..{n_obj} in order"
                )
            counts += np.asarray(data["CCD_N"]) >= 0
    return counts


def main() -> None:
    """Emit this chunk's range as the two bash exports the config reads."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--chunk", required=True, type=int)
    p.add_argument("--n-chunks", required=True, type=int)
    a = p.parse_args()
    # Explicit, because a negative --chunk would index from the end and hand
    # this process some other chunk's range without any error.
    if not 1 <= a.chunk <= a.n_chunks:
        raise SystemExit(
            f"[ngmix_range] FATAL: --chunk {a.chunk} outside 1..{a.n_chunks}"
        )
    lo, hi = id_ranges(object_epochs(a.run_dir), a.n_chunks)[a.chunk - 1]
    print(f"export NGMIX_ID_MIN={lo}; export NGMIX_ID_MAX={hi}")


if __name__ == "__main__":
    main()
