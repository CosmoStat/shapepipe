#!/usr/bin/env python3
"""Build the run index (``run_index.sqlite``) that drives the compute DAG.

The index is *parse-time data*, never a rule input: the Snakefile loads it once
at parse time into plain dicts, so extending a run's tile list changes which
jobs exist without touching the mtime chain of completed work.

It records, per run:

    tiles(tile_id, ra_dir, n_exp, status)
    exposures(exp_id)                          -- deduplicated union over tiles
    tile_exposures(tile_id, exp_id)            -- the tile->exposure edges

The tile->exposure edges are *data-derived*: they are read from each tile's
``find_exposures`` output (``exp_numbers-<IDra>-<IDdec>.txt``), which
``find_exposures_runner`` produces by parsing the tile FITS ``HISTORY`` header.
So this is a **plain script**, not a DAG node: it runs between the two static
invocations — after ``snakemake prepare_tiles`` (Git+Uz+Fe, keep-going) and
before ``snakemake prepare_exposures`` parses the now-existing index. It
iterates the *declared* tile list and checks each tile's Fe output at its
deterministic path (no globbing — O(tile-list) existence checks, respecting the
no-``ls``-at-scale ban). A bad tile costs only that tile: it is recorded in
``missing.json`` and the index is built over the rest. The script exits nonzero
only if the missing *fraction* exceeds ``--missing-threshold`` (default 0.0),
so a keep-going download storm that lost a few tiles does not cost the run.
There is NO ``--allow-missing`` flag (that Snakemake flag does not exist).

Exposure IDs are stored with their trailing single-char suffix stripped
(``2243881p`` -> ``2243881``); that ``exp_base`` is the dedup key and the
exposure-rule wildcard, matching the ``exp/<prefix>/<expbase>/`` forest layout.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path


def read_exposure_list(exp_numbers_file: Path) -> list[str]:
    """Return the exp_base IDs listed in one tile's find_exposures output.

    Each line is an exposure name like ``2243881p``; we strip a trailing
    alphabetic suffix so the base ID is the dedup key.
    """
    ids = []
    for line in exp_numbers_file.read_text().splitlines():
        name = line.strip()
        if not name:
            continue
        ids.append(name[:-1] if name[-1].isalpha() else name)
    return ids


def build(tile_ids: list[str], run_dir: Path, db_path: Path,
          missing_threshold: float = 0.0) -> dict:
    """Build the index over ``tile_ids``; return a summary dict.

    For each tile, check its ``exp_numbers-<IDra>-<IDdec>.txt`` at the tile's
    deterministic ``find_exposures`` run dir (RUN_DATETIME=False, no glob). A
    tile whose exposure list is missing is recorded in ``missing.json`` and the
    index is built over the rest (a bad tile costs that tile, not the run). The
    build is fatal only if the missing fraction exceeds ``missing_threshold``.
    """
    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path)
    con.executescript(
        """
        DROP TABLE IF EXISTS tiles;
        DROP TABLE IF EXISTS exposures;
        DROP TABLE IF EXISTS tile_exposures;
        CREATE TABLE tiles(tile_id TEXT PRIMARY KEY, ra_dir TEXT, n_exp INTEGER);
        CREATE TABLE exposures(exp_id TEXT PRIMARY KEY);
        CREATE TABLE tile_exposures(
            tile_id TEXT, exp_id TEXT,
            PRIMARY KEY (tile_id, exp_id));
        """
    )

    missing = []
    all_exposures: set[str] = set()
    for tile_id in tile_ids:
        ra_dir = tile_id.split(".")[0]
        idra, iddec = tile_id.split(".")
        # Flat forest, deterministic run dir (RUN_DATETIME=False).
        exp_file = (run_dir / "tiles" / tile_id / "output" / "run_sp_tile_Fe" /
                    "find_exposures_runner" / "output" /
                    f"exp_numbers-{idra}-{iddec}.txt")
        if not exp_file.exists():
            missing.append(tile_id)
            continue
        exp_ids = read_exposure_list(exp_file)
        con.execute("INSERT INTO tiles VALUES (?,?,?)",
                    (tile_id, ra_dir, len(exp_ids)))
        for exp_id in exp_ids:
            all_exposures.add(exp_id)
            con.execute("INSERT OR IGNORE INTO tile_exposures VALUES (?,?)",
                        (tile_id, exp_id))

    con.executemany("INSERT OR IGNORE INTO exposures VALUES (?)",
                    [(e,) for e in sorted(all_exposures)])
    con.commit()
    con.close()

    (db_path.parent / "missing.json").write_text(json.dumps(missing, indent=2))
    frac = len(missing) / len(tile_ids) if tile_ids else 0.0
    if frac > missing_threshold:
        raise SystemExit(
            f"Missing exposure lists for {len(missing)}/{len(tile_ids)} tile(s) "
            f"(fraction {frac:.3f} > threshold {missing_threshold}): {missing}. "
            f"Re-run prepare_tiles for them, or raise --missing-threshold.")
    return {"n_tiles": len(tile_ids) - len(missing),
            "n_exposures": len(all_exposures),
            "n_missing": len(missing)}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tile-list", required=True, type=Path,
                   help="file of tile IDs, one per line")
    p.add_argument("--run-dir", required=True, type=Path,
                   help="$SP_RUN: root of the tiles/ work-dir forest")
    p.add_argument("--db", required=True, type=Path,
                   help="output run_index.sqlite path")
    p.add_argument("--missing-threshold", type=float, default=0.0,
                   help="fatal if the missing-tile fraction exceeds this "
                        "(default 0.0: any missing tile is fatal)")
    args = p.parse_args()

    tile_ids = [ln.strip() for ln in args.tile_list.read_text().splitlines()
                if ln.strip()]
    summary = build(tile_ids, args.run_dir, args.db, args.missing_threshold)
    print(f"run_index: {summary['n_tiles']} tiles, "
          f"{summary['n_exposures']} exposures, "
          f"{summary['n_missing']} missing -> {args.db}", file=sys.stderr)


if __name__ == "__main__":
    main()
