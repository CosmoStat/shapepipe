#!/usr/bin/env python3
"""Archive PR #922 epoch-cut diagnostics without retaining module logs.

Chunk completeness checks read the single end-of-loop summary before temp()
removes the chunk directory. The make-cat check gathers those structured chunk
records, and its manifest is published beside final_cat. Backfills use the same
schema and parser. Missing or ambiguous records are errors, never zero counts.
These are diagnostics only: no cuts or catalogue values are changed.
"""

import json
import re
from pathlib import Path

FIELDS = (
    "considered",
    "masked_fraction",
    "central_veto",
    "failed",
    "objects_emptied",
)


def validate_counts(counts):
    """Validate an exact, non-negative integer count record."""
    if not isinstance(counts, dict) or set(counts) != set(FIELDS):
        raise ValueError(f"epoch cuts must contain exactly {FIELDS}")
    if any(type(v) is not int or v < 0 for v in counts.values()):
        raise ValueError("epoch cuts must be non-negative integers")
    if (
        sum(counts[k] for k in ("masked_fraction", "central_veto", "failed"))
        > counts["considered"]
    ):
        raise ValueError("epoch cuts exceed epochs considered")
    if counts["objects_emptied"] > counts["considered"]:
        raise ValueError("objects emptied exceeds epochs considered")
    return dict(counts)


def parse_epoch_cuts(text):
    """Parse one complete end-of-loop summary; reject retry ambiguity."""
    lines = [
        line.split("epoch cuts:", 1)[1].strip()
        for line in text.splitlines()
        if "epoch cuts:" in line
    ]
    if len(lines) != 1:
        raise ValueError(
            f"expected one epoch cuts summary, found {len(lines)}"
        )
    counts = {}
    for token in lines[0].split():
        match = re.fullmatch(r"([a-z_]+)=(\d+)", token)
        if not match or match[1] in counts:
            raise ValueError(f"invalid epoch cuts token: {token!r}")
        counts[match[1]] = int(match[2])
    return validate_counts(counts)


def read_chunk_logs(run_root):
    """Read ngmix worker logs only, not copies or ShapePipe's aggregate log."""
    paths = sorted(
        (Path(run_root) / "ngmix_runner" / "logs").glob("process-*.log")
    )
    if len(paths) != 1:
        raise ValueError(
            f"expected one ngmix worker log under {run_root}, "
            f"found {len(paths)}"
        )
    return parse_epoch_cuts(paths[0].read_text())


def aggregate_chunks(chunks, n_chunks):
    """Validate the complete 1..N chunk set and sum each diagnostic once."""
    if type(n_chunks) is not int or n_chunks < 1:
        raise ValueError("n_chunks must be a positive integer")
    expected = {str(k) for k in range(1, n_chunks + 1)}
    if set(chunks) != expected:
        raise ValueError(
            f"expected chunks 1..{n_chunks}, found {sorted(chunks)}"
        )
    chunks = {
        str(k): validate_counts(chunks[str(k)]) for k in range(1, n_chunks + 1)
    }
    return {
        "schema_version": 1,
        "n_chunks": n_chunks,
        "chunks": chunks,
        "totals": {
            key: sum(c[key] for c in chunks.values()) for key in FIELDS
        },
    }


def gather_manifests(tile_dir, n_chunks):
    """Gather successful chunk records after the merge DAG edge."""
    chunks = {}
    for k in range(1, n_chunks + 1):
        path = Path(tile_dir) / "manifests" / f"tile_ngmix_{k}.json"
        manifest = json.loads(path.read_text())
        if (
            manifest.get("stage") != "tile_ngmix"
            or manifest.get("unit") != Path(tile_dir).name
            or manifest.get("status") != "complete"
        ):
            raise ValueError(f"not a successful chunk manifest: {path}")
        chunks[str(k)] = manifest.get("epoch_cuts")
    return aggregate_chunks(chunks, n_chunks)


def require_durable_manifest(final_cat, tile, n_chunks=None):
    """Validate the persisted catalogue and its completeness/count record."""
    final_cat = Path(final_cat)
    if not final_cat.is_file() or final_cat.stat().st_size == 0:
        raise ValueError(f"no persisted final_cat: {final_cat}")
    path = final_cat.parent / "tile_make_cat.json"
    manifest = json.loads(path.read_text())
    if (
        manifest.get("stage") != "tile_make_cat"
        or manifest.get("unit") != tile
        or manifest.get("status") != "complete"
    ):
        raise ValueError(
            f"not a successful make-cat manifest for {tile}: {path}"
        )
    record = manifest.get("epoch_cuts") or {}
    expected = aggregate_chunks(
        record.get("chunks", {}),
        n_chunks if n_chunks is not None else record.get("n_chunks"),
    )
    if record != expected:
        raise ValueError(f"invalid epoch cuts aggregate: {path}")
    return manifest
