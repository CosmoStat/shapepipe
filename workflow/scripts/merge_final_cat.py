#!/usr/bin/env python3
"""Collect the campaign's per-tile final catalogues into one HDF5 file.

Run as the shell of the campaign-level ``final_cat_merge`` rule, never by hand.

``<products_dir>/final_cat_<run>.hdf5`` holds
one dataset per tile, carrying the columns named by the input type's
``final_cat.param`` (``workflow/config/cfis/`` for data,
``workflow/config/cfis_image_sims/`` for image sims), plus an ``n_tiles`` attribute on the
file root. sp_validation opens that file as its ``galaxy_cat_path``
(``sp_validation/catalog.py``); see ``spval_group`` for the group-path interface.

(sp_validation's own ``merge_catalogues`` is a different layer entirely: it
works over already-calibrated ``shape_catalog_comprehensive_*.fits``. It does
not do this merge, and this does not do that one.)

Column extraction uses ``scripts/python/create_final_cat.py``'s
``read_param_file``, ``read_data`` and ``copy_data``; those functions define
column selection and ordering. Discovery uses the workflow's products tree:
``tiles/<2-char prefix>/<ID>/final_cat-<ID>.hdf5``, one dataset per column.

``create_final_cat.py`` is loaded by path relative to this file, so it follows
``bin/sp``'s launch code snapshot. It is a script, not an installed module.

Reconciliation uses ``hdf5_reconcile.py``, shared with the star-catalogue merge;
see that module for source tracking and update semantics. This script does not
use ``create_final_cat.py``'s append-only ``process()``.

``build_index.campaign_tiles`` derives membership from ``--tile-list`` and
``--index-db`` rather than passing ~20k paths in a shell argument (limited to
128 KiB by ``MAX_ARG_STRLEN``). See the Snakefile's ``final_cat_merge`` rule
for dependencies and the membership fingerprint.

@sc [label:coupling] final-merge-campaign-membership
Merge only tiles in the tile-list/index intersection, not a glob of the shared
products root: unrelated tiles would evade the rule's membership fingerprint.
A selected tile with no catalogue is a hard error, not a skip.

@sc [decision:catalogue_assembly.failure_sentinels,label:selection] never-fit-rows-pass-through
Every row of every tile catalogue reaches the merged file, unchanged,
including objects ngmix never fit. Those carry `NGMIX_N_EPOCH == 0` with
sentinel values (ellipticities `-10`, `T == 0`), `NGMIX_MCAL_FLAGS` nonzero
(`LM_FUNC_NOTFINITE`) and `NGMIX_MCAL_TYPES_FAIL == 5`, so the consumer's
`NGMIX_MCAL_FLAGS == 0` cut rejects them. The merge neither fills these rows
nor drops them. Enforced by tests/unit/test_final_cat_merge_invariants.py.
"""

import argparse
import importlib.util
import sys
from pathlib import Path

# Same directory; the rule invokes this file by path, so it is sys.path[0].
import build_index
import hdf5_reconcile

# <repo>/scripts/python/create_final_cat.py, from <repo>/workflow/scripts/this.
CFC_PATH = (Path(__file__).resolve().parents[2]
            / "scripts" / "python" / "create_final_cat.py")


def spval_group(campaign: str) -> str:
    """The hdf5 group the campaign's per-tile datasets live under.

    @sc [label:schema] final-merge-spval-group-path
    Keep ``patches/<campaign>`` as the group path required by sp_validation's
    reader. CosmoStat/sp_validation#340 tracks the reader migration; ``patches``
    is a file-schema key, not a workflow unit.
    """
    return f"patches/{campaign}"


def load_create_final_cat():
    """The hdf5 layout's definition, loaded by path (see the module docstring)."""
    if not CFC_PATH.exists():
        sys.exit(f"merge_final_cat: {CFC_PATH} is not there — the launch code "
                 f"snapshot must carry scripts/python/ (see bin/sp).")
    spec = importlib.util.spec_from_file_location("create_final_cat", CFC_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def catalogues(products_dir: Path, tile_list: Path, index_db: Path) -> list:
    """``(tile ID, path)`` for the campaign's tiles, in ID order.

    Not a glob over the products root: see the module docstring on why the set
    is the campaign's and not the filesystem's.
    """
    out, missing = [], []
    for tile in sorted(build_index.campaign_tiles(tile_list, index_db)):
        path = (products_dir / "tiles" / tile[:2] / tile
                / f"final_cat-{tile}.hdf5")
        if path.exists():
            out.append((tile, path))
        else:
            missing.append(tile)
    if missing:
        sys.exit(f"merge_final_cat: {len(missing)} campaign tile(s) have no "
                 f"final catalogue: {' '.join(missing[:5])}"
                 f"{' ...' if len(missing) > 5 else ''}")
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--products-dir", required=True, type=Path,
                   help="the persistent root; per-tile catalogues are found "
                        "beneath it")
    p.add_argument("--tile-list", required=True, type=Path,
                   help="the campaign's tile list (config tile_list)")
    p.add_argument("--index-db", required=True, type=Path,
                   help="the campaign's run index (config outputs.index_db)")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--campaign", required=True,
                   help="names the campaign's group in the output file")
    p.add_argument("--param-file", required=True, type=Path,
                   help="the input type's final_cat.param — the column list")
    p.add_argument("--snapshot-json", type=Path, default=None,
                   help="sp run's code snapshot (bin/sp's "
                        "$STATE_DIR/code/snapshot.json); absent outside sp run")
    args = p.parse_args()

    cfc = load_create_final_cat()
    param_list = cfc.read_param_file(str(args.param_file), verbose=False)
    if not param_list:
        sys.exit(f"merge_final_cat: no columns read from {args.param_file}")
    # read_data/copy_data read their knobs out of this dict, exactly as
    # create_final_cat.py's own main() builds it.
    params = {"param_list": param_list, "verbose": False}

    tiles = catalogues(args.products_dir, args.tile_list, args.index_db)
    if not tiles:
        # An empty hdf5 would satisfy every downstream existence check and
        # produce an empty shear catalogue.
        sys.exit(f"merge_final_cat: no tile in {args.tile_list} is indexed in "
                 f"{args.index_db}, so there is nothing to merge")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    group_path = spval_group(args.campaign)
    digest = hdf5_reconcile.schema_digest(param_list)

    def read_tile(tile, path):
        """One tile's requested columns, via create_final_cat.py's own reader."""
        extracted, dtype = cfc.read_data(str(path), params)
        return cfc.copy_data(params["param_list"], extracted, dtype)

    todo = hdf5_reconcile.plan(args.output, group_path, tiles, digest)
    if todo.empty():
        print(f"[merge_final_cat] unchanged: {args.output} "
              f"({len(tiles)} tile(s))")
        return
    hdf5_reconcile.apply(args.output, group_path, todo, tiles, read_tile,
                         digest, "n_tiles",
                         hdf5_reconcile.code_provenance(args.snapshot_json))
    print(f"[merge_final_cat] {todo.describe()} -> {args.output} "
          f"({len(tiles)} tile(s), {len(param_list)} column(s), "
          f"group {group_path})")


if __name__ == "__main__":
    main()
