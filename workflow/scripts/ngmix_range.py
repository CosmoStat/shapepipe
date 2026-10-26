#!/usr/bin/env python3
"""Materialise the ngmix chunk partition once per tile, then serve it.

``tile_vignets`` computes the partition and writes it to the ``tile_shape``
group's node-local scratch. Each chunk shell reads its own row::

    # tile_vignets, once (see tile.smk's TILE_NGMIX_RANGES)
    ngmix_range.py --run-dir $SP_RUN --n-chunks 8 --write $SP_LOCAL/ngmix_ranges.json

    # each chunk, at its own start
    eval "$(ngmix_range.py --read $SP_LOCAL/ngmix_ranges.json --chunk 3)"
    # -> export NGMIX_ROW_MIN=751; export NGMIX_ROW_MAX=1125

The file lives on ``$SP_LOCAL`` for the group job's lifetime and is not a rule
output; see ``workflow/rules/tile.smk`` for group scratch ownership. Ranges are
computed at execution time because Snakemake params evaluate before the sexcat
exists.

@sc [label:coupling] ngmix-range-row-partition
Ranges are closed, 1-based sexcat row positions, not NUMBER values, matching
``ngmix_package.ngmix.chunk_rows``. They must cover [1, N] exactly once.
Never use a nonpositive upper bound: ngmix treats it as unbounded and can
re-measure the whole tile.

Chunks balance epoch-weighted cost because ngmix fits each object jointly
across its exposures and the group waits for its slowest chunk. Measurements
on a 34-tile nibi run give 0.2714 CPU-s per (object, epoch) pair plus ~0.05
CPU-s of per-object setup. Equal-count chunks vary in cost by 1.6x within a
tile; epoch-weighted splitting reduces the sum of predicted slowest-chunk
costs from 132,875 to 113,710 CPU-s (14.4%).

Moving the boundaries does not change measurements: ngmix's RNG is seeded per object
from its sky position, so which chunk an object falls in cannot change its
measurement (see ``ngmix_package.ngmix.position_seed``).

@sc [label:coupling] ngmix-range-single-writer
All sibling chunks must read the partition written once by ``tile_vignets``.
Read mode must fail rather than recompute missing ranges: disagreeing chunk
boundaries silently duplicate or drop objects, and ``merge_sep_cats`` merely
concatenates the results.

Integer weights make partitioning reproducible across attempts. Missing EPOCH
extensions fail write mode. See the Snakefile's ``NGMIX_RANGE_HASH`` for
invalidation across resumes after splitter edits, and ``bin/sp`` for launch
code snapshots.
"""

import argparse
import json
from pathlib import Path

# Integer milli-epoch weights make boundary comparisons exact and reproducible.
MILLI_EPOCH = 1000

# Setup cost in epoch-equivalents: 0.05 / 0.2714 = 0.184. This gives
# zero-epoch objects nonzero weight (~6% of a typical object). On tile 186.307,
# refitting with geometric epochs gives 0.191 and moves boundaries by at most
# one object; zeroing the setup term moves them by at most 44.
ALPHA_MILLI_EPOCHS = 184


def row_ranges(epochs, n_chunks: int) -> list[tuple[int, int]]:
    """Split catalogue rows ``1..len(epochs)`` into ``n_chunks`` closed ranges.

    ``epochs[i]`` is row ``i + 1``'s geometric epoch count. The ranges are
    contiguous — ``NGMIX_ROW_MIN``/``NGMIX_ROW_MAX`` is an interval, not a set —
    and tile ``[1, n_obj]`` exactly, so every object is measured once.

    The objective is the slowest chunk, not the average one, because the
    group job waits for it. So this minimises the maximum chunk weight
    exactly: binary-search the smallest feasible capacity, then fill left to
    right under it. Ties in that maximum break toward the earlier chunks,
    which is arbitrary but fixed, and fixed is the property that matters.

    With ``n_obj < n_chunks`` the first ``n_obj`` chunks take one object each
    and the remainder are empty, written ``(n_obj + 1, n_obj)``. This preserves
    ``ranges[k][0] == ranges[k-1][1] + 1`` and keeps the upper bound positive
    as required by the module's row-partition contract.
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
    discovered by name, never assumed — each with ``n_obj`` rows in
    ``LDAC_OBJECTS`` row order and ``CCD_N < 0`` where the object misses that
    exposure. Summing
    ``CCD_N >= 0`` across them reproduces the sexcat's geometric
    ``N_EPOCH`` — verified on tile 186.307: 35,298 rows, 7 extensions,
    116,727 pairs, mean 3.31. The post-process is upstream of the whole
    tile_shape group, so the extensions always exist by the time ngmix runs;
    their absence is a broken tile, not a case to accommodate.

    This geometric count mildly overstates the work, because ngmix
    fits only epochs with a validated PSF and
    drops epochs it cannot fit. The same catalogue's NGMIX_N_EPOCH is lower for
    2,203 of the 35,298 objects and never higher: 113,947 pairs against
    116,727, a 2.4% difference. A fit using NGMIX_N_EPOCH gives 0.2681 CPU-s
    per pair at R^2 0.997, against 0.2619 at R^2 0.973 using geometric counts.
    Fitted epoch counts are unavailable before ngmix; on this tile they improve
    the predicted slowest chunk by only 1.07%.

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
        counts = np.zeros(n_obj, dtype=np.int64)
        number = None
        for hdu in epoch_hdus:
            data = hdu.data
            # Row alignment is asserted, not assumed: the weights are indexed
            # by row, so every EPOCH extension must carry the same NUMBER
            # sequence, row for row, as the first one.
            this_number = np.asarray(data["NUMBER"], dtype=np.int64)
            if number is None:
                number = this_number
            if len(data) != n_obj or not np.array_equal(this_number, number):
                raise SystemExit(
                    f"[ngmix_range] FATAL: {cats[0]}[{hdu.name}] is not "
                    f"{n_obj} rows aligned with {epoch_hdus[0].name}"
                )
            counts += np.asarray(data["CCD_N"]) >= 0
    return counts


def partition(epochs, n_chunks: int) -> dict:
    """The whole partition, as the dict that gets serialised to JSON.

    Minimal on purpose: the ranges themselves, plus the two numbers a reader
    needs to check that the file is the one it expects (``n_obj`` so a mismatch
    against the catalogue is visible, ``n_chunks`` so a chunk index is validated
    against what was actually written, not against what the caller believes).
    ``chunk`` is 1-based, matching SP_NGMIX_CHUNK and the run-directory suffix.
    """
    ranges = row_ranges(epochs, n_chunks)
    return {
        "n_obj": len(epochs),
        "n_chunks": n_chunks,
        "chunks": [
            {"chunk": k, "row_min": lo, "row_max": hi}
            for k, (lo, hi) in enumerate(ranges, start=1)
        ],
    }


def chunk_range(doc: dict, chunk: int) -> tuple[int, int]:
    """Chunk ``chunk``'s closed range out of a ``partition()`` document.

    Validated rather than indexed: a negative ``chunk`` would index from the end
    and hand this process some other chunk's range without any error.
    """
    n_chunks = doc["n_chunks"]
    if not 1 <= chunk <= n_chunks:
        raise SystemExit(
            f"[ngmix_range] FATAL: --chunk {chunk} outside 1..{n_chunks}"
        )
    row = doc["chunks"][chunk - 1]
    if row["chunk"] != chunk:
        raise SystemExit(
            f"[ngmix_range] FATAL: ranges file row {chunk - 1} is chunk "
            f"{row['chunk']}, not {chunk}"
        )
    return int(row["row_min"]), int(row["row_max"])


def read_ranges(path: Path) -> dict:
    """Load the ranges file, or die saying what should have written it."""
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        raise SystemExit(
            f"[ngmix_range] FATAL: no ranges file at {path}.\n"
            "  tile_vignets writes it once per tile, on the node-local scratch "
            "this group job shares.\n"
            "  Its absence means tile_vignets did NOT run in this group job "
            "(same cause as a missing vignette store).\n"
            "  FIX: rm the tile's tile_vignets.json and resume.\n"
            "  NOT recomputed here on purpose: a per-chunk fallback would let "
            "chunks disagree about the partition."
        )


def main() -> None:
    """Write the whole partition, or emit one chunk's range as bash exports."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", type=Path, metavar="JSON",
                   help="compute the whole partition and write it here")
    p.add_argument("--read", type=Path, metavar="JSON",
                   help="look one chunk's range up in a file --write made")
    p.add_argument("--run-dir", type=Path, help="--write: the tile's run root")
    p.add_argument("--n-chunks", type=int, help="--write: chunks to split into")
    p.add_argument("--chunk", type=int, help="--read: 1-based chunk index")
    a = p.parse_args()

    if bool(a.write) == bool(a.read):
        raise SystemExit(
            "[ngmix_range] FATAL: pass exactly one of --write / --read"
        )

    if a.write:
        if a.run_dir is None or a.n_chunks is None:
            raise SystemExit(
                "[ngmix_range] FATAL: --write needs --run-dir and --n-chunks"
            )
        doc = partition(object_epochs(a.run_dir), a.n_chunks)
        # Temp name plus rename, the same all-or-nothing publish tile_local()
        # uses for the WCS store: a reader must never see a half-written file.
        tmp = a.write.with_name(f".{a.write.name}.tmp")
        tmp.write_text(json.dumps(doc))
        tmp.replace(a.write)
        return

    if a.chunk is None:
        raise SystemExit("[ngmix_range] FATAL: --read needs --chunk")
    lo, hi = chunk_range(read_ranges(a.read), a.chunk)
    print(f"export NGMIX_ROW_MIN={lo}; export NGMIX_ROW_MAX={hi}")


if __name__ == "__main__":
    main()
