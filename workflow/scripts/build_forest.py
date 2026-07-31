#!/usr/bin/env python3
"""Build a tile's exposure symlink forest (the SP_EXP view).

A plain script (not a run: block) so the tile chain stays group-compatible.
Reads the tile's exposures from run_index.sqlite and symlinks each exposure's
``exp/<prefix>/<base>/output`` into ``<forest>/<prefix>/<base>/output`` by exact
name (no glob). The 2-char ``<prefix>`` shard level is NOT cosmetic: ShapePipe's
``exp_utils.get_exp_output_files`` hardwires the sharded v2.0 layout into its
$SP_EXP glob (``<SP_EXP>/<prefix>/<base>/output/run_sp_*/...``), so a flat
forest makes every tile gather stage fail "No split_exp_runner output found".
(The exposure STORE is sharded the same way, for the filesystem's sake.)
The forest is a convenience view; the DAG edge to the exposures is declared in
the rule's input (tile.smk), not here.
"""

import argparse
import shutil
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
        src = args.run_dir / "exp" / e[:2] / e / "output"
        dst = args.forest / e[:2] / e / "output"   # sharded: the module glob's shape
        dst.parent.mkdir(parents=True, exist_ok=True)
        # A symlink (the normal case) is unlinked; a REAL directory left behind
        # by a hand-run or an older layout must be removed as a tree — unlink()
        # raises IsADirectoryError on it and would kill the job.
        if dst.is_symlink() or dst.exists():
            if dst.is_dir() and not dst.is_symlink():
                shutil.rmtree(dst)
            else:
                dst.unlink()
        dst.symlink_to(src)
    print(f"[build_forest] {args.tile}: {len(exps)} exposures -> {args.forest}")


if __name__ == "__main__":
    main()
