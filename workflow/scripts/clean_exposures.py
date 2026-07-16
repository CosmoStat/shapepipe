#!/usr/bin/env python3
"""Index-driven exposure-store cleanup (`sp clean-exposures`).

temp() is deliberately NOT used on the shared exposure directories (findings
4/5/6/19): they overlap tiles by construction, so temp() deletion would cascade
destructive reruns across spatial neighbours whenever a tile is appended. Instead
the store is persistent and space is reclaimed explicitly here.

For each exposure in run_index.sqlite, delete its bulky intermediates
(run_sp_exp_Sp / run_sp_exp_Ma / run_sp_exp_SxSePsfPi) ONLY when every tile that
consumes it (tile_exposures) has its final_cat-<ID>.fits present across the whole
output root. An exposure still feeding an unfinished tile is left untouched.

Caveat: if a later batch appends a tile that shares a deleted exposure,
Snakemake will regenerate that exposure — cleanup trades peak storage for
recompute. Run it between campaign phases, not mid-DAG. --dry-run by default.
"""

import argparse
import shutil
import sqlite3
import sys
from pathlib import Path

PREFIXES = ("run_sp_exp_Sp", "run_sp_exp_Ma", "run_sp_exp_SxSePsfPi")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--index", required=True, type=Path)
    p.add_argument("--apply", action="store_true",
                   help="actually delete (default: dry-run, list only)")
    args = p.parse_args()

    con = sqlite3.connect(args.index)
    edges = con.execute("SELECT exp_id, tile_id FROM tile_exposures").fetchall()
    con.close()

    consumers: dict[str, list[str]] = {}
    for exp_id, tile_id in edges:
        consumers.setdefault(exp_id, []).append(tile_id)

    def tile_done(tile):
        return (args.run_dir / "tiles" / tile / f"final_cat-{tile}.fits").exists()

    freed, held = 0, 0
    for exp_id, tiles in consumers.items():
        if not all(tile_done(t) for t in tiles):
            held += 1
            continue
        for prefix in PREFIXES:
            d = args.run_dir / "exp" / exp_id / "output" / prefix
            if d.is_dir():
                freed += 1
                print(f"{'DELETE' if args.apply else 'would delete'} {d}")
                if args.apply:
                    shutil.rmtree(d)
    print(f"[clean_exposures] {freed} dirs {'freed' if args.apply else 'reclaimable'}, "
          f"{held} exposures held (consuming tile unfinished)", file=sys.stderr)


if __name__ == "__main__":
    main()
