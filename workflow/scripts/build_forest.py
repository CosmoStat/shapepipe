#!/usr/bin/env python3
"""Build a tile's exposure symlink forest (the SP_EXP view).

A plain script (not a run: block) so the tile chain stays group-compatible.
Reads the tile's exposures from run_index.sqlite and symlinks each exposure's
``exp/<base>/output`` into ``<forest>/<base>/output`` by exact name (no glob), so
the tile modules' $SP_EXP globs resolve to exactly this tile's exposures'
products. The forest is a convenience view; the DAG edge to the exposures is
declared in the rule's input (tile.smk), not here.
"""

import argparse
import sqlite3
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tile", required=True)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--index", required=True, type=Path)
    p.add_argument("--forest", required=True, type=Path)
    args = p.parse_args()

    con = sqlite3.connect(args.index)
    exps = [r[0] for r in con.execute(
        "SELECT exp_id FROM tile_exposures WHERE tile_id=?", (args.tile,))]
    con.close()

    args.forest.mkdir(parents=True, exist_ok=True)
    for e in exps:
        src = args.run_dir / "exp" / e / "output"
        dst = args.forest / e / "output"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src)
    print(f"[build_forest] {args.tile}: {len(exps)} exposures -> {args.forest}")


if __name__ == "__main__":
    main()
