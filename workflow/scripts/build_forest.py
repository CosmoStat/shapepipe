#!/usr/bin/env python3
"""Build a tile's exposure symlink forest (the SP_EXP view).

A plain script rather than a ``run:`` block keeps the tile chain group-compatible.
It reads the tile's exposures from run_index.sqlite and creates links by exact
name, without globs. Links for current exposures are replaced; entries absent
from the current index query are not pruned.
The exposure dependencies belong to the rule's input in
``workflow/rules/tile.smk``, not to this convenience view.

@sc [label:coupling] exposure-forest-sharded-layout
Link ``exp/<prefix>/<base>/output`` at ``<forest>/<prefix>/<base>/output``, with
the first two base-ID characters as the shard. ``exp_utils.get_exp_output_files``
expects this layout in its ``$SP_EXP`` glob; a flat forest makes gather stages
fail to find exposure outputs.
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

    con = sqlite3.connect(args.index, timeout=60)
    exps = [r[0] for r in con.execute(
        "SELECT exp_id FROM tile_exposures WHERE tile_id=?", (args.tile,))]
    con.close()

    args.forest.mkdir(parents=True, exist_ok=True)
    for e in exps:
        src = args.run_dir / "exp" / e[:2] / e / "output"
        dst = args.forest / e[:2] / e / "output"   # sharded: the module glob's shape
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Unlink existing links; remove real directories as trees because
        # unlink() cannot delete them.
        if dst.is_symlink() or dst.exists():
            if dst.is_dir() and not dst.is_symlink():
                shutil.rmtree(dst)
            else:
                dst.unlink()
        dst.symlink_to(src)
    print(f"[build_forest] {args.tile}: {len(exps)} exposures -> {args.forest}")


if __name__ == "__main__":
    main()
