#!/usr/bin/env python3
"""Build the run index (``run_index.sqlite``) that drives the compute DAG.

The index is *parse-time data*, never a rule input: the Snakefile loads it once
at parse time into plain dicts, so extending a run's tile list changes which
jobs exist without touching the mtime chain of completed work.

It ACCUMULATES ACROSS INVOCATIONS, but is AUTHORITATIVE FOR THE CURRENT TILE
LIST (D1). Precisely:

  * a tile in the current list that has find_exposures output has its ``tiles``
    row replaced and its ``tile_exposures`` edges DELETED AND REBUILT — the Fe
    output is the truth, so a tile whose exposure list shrank or changed must
    not keep stale edges to exposures it no longer reads;
  * a tile NOT in the current list is left completely untouched, which is what
    makes the index span the campaign and lets a later ``clean_exposure`` see
    every consuming tile;
  * ``exposures`` rows only ever accumulate (``INSERT OR IGNORE``). An exposure
    no tile references any more is a harmless orphan: nothing reads that table
    except by joining through ``tile_exposures``.

It records:

    tiles(tile_id, ra_dir, n_exp, status)
    exposures(exp_id)                          -- deduplicated union over tiles
    tile_exposures(tile_id, exp_id)            -- the tile->exposure edges

The tile->exposure edges are *data-derived*: they are read from each tile's
``find_exposures`` output (``exp_numbers-<IDra>-<IDdec>.txt``), which
``find_exposures_runner`` produces by parsing the tile FITS ``HISTORY`` header.
So the build is not a DAG node: ``build()`` is imported and called at PARSE TIME
by the Snakefile of the COMPUTE invocation ONLY — ``SP_PHASE == "compute"``,
which ``bin/sp`` sets — after the PREPARE invocation has produced the tiles'
find_exposures output. No other parse builds anything: a prepare parse or a
passthrough invocation (``sp --unlock``, ``sp --dag``) just loads whatever is
already on disk. (There is no ``sp index`` verb; the CLI below stays for
hand-inspection.) ``build`` below documents the missing-tile policy and the
write ordering.

Exposure IDs are stored with their trailing single-char suffix stripped
(``2243881p`` -> ``2243881``); that ``exp_base`` is the dedup key and the
exposure-rule wildcard, matching the sharded ``exp/<prefix>/<expbase>/`` store.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path


def read_exposure_list(exp_numbers_file: Path) -> list[tuple[str, str]]:
    """Return ``(exp_id, name)`` pairs from one tile's find_exposures output.

    Each line is an exposure *name* like ``2243881p``; the bare base ID (suffix
    stripped) is the dedup key everywhere in the DAG, but the original name is
    kept in the index — the fabricated per-unit ``exp_numbers`` list must carry
    it verbatim (``get_images`` matches ``<name>.fits.fz`` in the store; the
    bare ID matches nothing).
    """
    pairs = []
    for line in exp_numbers_file.read_text().splitlines():
        name = line.strip()
        if not name:
            continue
        pairs.append((name[:-1] if name[-1].isalpha() else name, name))
    return pairs


def exp_list_path(run_dir: Path, tile_id: str) -> Path:
    """This tile's find_exposures output, at its deterministic path.

    Sharded store (D2), fixed run dir (RUN_DATETIME=False) — an existence check,
    never a glob (no ``ls`` at scale).
    """
    idra, iddec = tile_id.split(".")
    return (run_dir / "tiles" / tile_id[:2] / tile_id / "output" /
            "run_sp_tile_Fe" / "find_exposures_runner" / "output" /
            f"exp_numbers-{idra}-{iddec}.txt")


def build(tile_ids: list[str], run_dir: Path, db_path: Path,
          missing_threshold: float | None = 0.0) -> dict:
    """Build the index over ``tile_ids``; return a summary dict.

    For each tile, check its ``exp_numbers-<IDra>-<IDdec>.txt`` at the tile's
    deterministic ``find_exposures`` run dir (RUN_DATETIME=False, no glob). A
    tile whose exposure list is missing is recorded in ``missing.json`` and the
    index is built over the rest (a bad tile costs that tile, not the run). The
    build is fatal only if the missing fraction exceeds ``missing_threshold``.

    ORDER MATTERS: the threshold is evaluated FIRST, from the missing set, and
    the database + ``missing.json`` are written only if it passes. A build that
    aborts must leave no trace — an aborted parse that had already mutated
    durable state was the bug this ordering fixes.

    The write is idempotent, which is what makes it acceptable that the compute
    parse runs it even under ``-n``: re-running over an unchanged tree produces
    an identical database.
    """
    missing = [t for t in tile_ids if not exp_list_path(run_dir, t).exists()]
    frac = len(missing) / len(tile_ids) if tile_ids else 0.0
    if missing_threshold is not None and frac > missing_threshold:
        raise SystemExit(
            f"Missing exposure lists for {len(missing)}/{len(tile_ids)} tile(s) "
            f"(fraction {frac:.3f} > threshold {missing_threshold}): {missing}. "
            f"Re-run prepare_tiles for them, or raise --missing-threshold.")

    db_path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(db_path, timeout=60)
    # No DROP: the index accumulates across invocations (D1).
    con.executescript(
        """
        CREATE TABLE IF NOT EXISTS tiles(
            tile_id TEXT PRIMARY KEY, ra_dir TEXT, n_exp INTEGER);
        CREATE TABLE IF NOT EXISTS exposures(
            exp_id TEXT PRIMARY KEY, name TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS tile_exposures(
            tile_id TEXT, exp_id TEXT,
            PRIMARY KEY (tile_id, exp_id));
        """
    )

    missing_set = set(missing)
    all_exposures: set[tuple[str, str]] = set()
    for tile_id in tile_ids:
        if tile_id in missing_set:
            continue
        ra_dir = tile_id.split(".")[0]
        exp_pairs = read_exposure_list(exp_list_path(run_dir, tile_id))
        con.execute("INSERT OR REPLACE INTO tiles VALUES (?,?,?)",
                    (tile_id, ra_dir, len(exp_pairs)))
        # Replace this tile's edge set wholesale. INSERT OR IGNORE alone only
        # ever added, so a tile whose exposure list shrank kept edges to
        # exposures it no longer reads — and those stale edges would block it in
        # the report and pin those exposures against cleanup.
        con.execute("DELETE FROM tile_exposures WHERE tile_id = ?", (tile_id,))
        con.executemany("INSERT INTO tile_exposures VALUES (?,?)",
                        [(tile_id, exp_id) for exp_id, _ in exp_pairs])
        all_exposures.update(exp_pairs)

    # Exposures accumulate: OR IGNORE, never REPLACE (the name never changes,
    # and orphans left by a shrunken tile are harmless).
    con.executemany("INSERT OR IGNORE INTO exposures VALUES (?,?)",
                    sorted(all_exposures))
    con.commit()
    con.close()

    (db_path.parent / "missing.json").write_text(json.dumps(missing, indent=2))
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
