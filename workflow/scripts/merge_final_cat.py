#!/usr/bin/env python3
"""Collect the campaign's per-tile final catalogues into ONE hdf5 file.

Run as the shell of the campaign-level ``final_cat_merge`` rule, never by hand.

WHAT IT PRODUCES, AND FOR WHOM. ``<products_dir>/final_cat_<campaign>.hdf5``:
one dataset per tile, carrying the columns named by
``workflow/config/cfis/final_cat.param``, plus an ``n_tiles`` attribute on the
file root. sp_validation opens that file as its ``galaxy_cat_path``
(``sp_validation/catalog.py``), so its SCHEMA is an interface and not a choice —
see ``SPVAL_GROUP`` below for the one legacy literal in it.

(sp_validation's own ``merge_catalogues`` is a different layer entirely: it
works over already-calibrated ``shape_catalog_comprehensive_*.fits``. It does
not do this merge, and this does not do that one.)

WHAT IT REUSES, AND WHAT IT DOES NOT. The column extraction is
``create_final_cat.py``'s — ``read_param_file`` for the parameter list,
``read_data`` and ``copy_data`` for pulling those columns out of one catalogue
with their FITS dtypes — so the column grammar keeps exactly one definition.
Those three are REPRODUCIBLE FUNCTIONS, and this PR is what made them so: the
parameter list comes back ordered rather than through a set, ``copy_data``
allocates the requested columns alone rather than leaving every other column of
the source as uninitialised memory, and a missing column raises with its own
name instead of falling out of a bare ``except:`` as an UnboundLocalError. The
fixes are upstream, in that script, because a hand-run of it deserves them as
much as this rule does.
Its ``process()`` is NOT used and neither is any of its discovery: that function
walks a directory tree the workflow does not have and never will, and it groups
by a unit ShapePipe v2 no longer has. This script walks the workflow's own
products tree instead (``tiles/<2-char prefix>/<ID>/final_cat-<ID>.fits``) and
writes the hdf5 itself.

WHERE ``create_final_cat.py`` IS FOUND. Beside this workflow, at
``<repo>/scripts/python/create_final_cat.py`` — resolved relative to THIS file,
so it follows the launch code snapshot (``bin/sp``) exactly as
``workflow/scripts/*`` does, and a campaign never reads a mid-run edit. It is
loaded by path rather than imported: it is a script, not an installed module,
and the container's ``shapepipe`` install does not carry it.

IT RECONCILES, IT NEITHER REBUILDS NOR BLINDLY APPENDS, and the machinery for
that is ``hdf5_reconcile.py``, shared with the star side so the campaign's two
products cannot disagree about what an output owes its inputs. That module
carries the argument in full: an append reads the appended tiles, a source that
changed is re-read, a tile that left the campaign is deleted, a column-set
change refreshes everything, and a no-op leaves the file untouched.
``create_final_cat.py``'s own ``process()`` implements only the append-only half
— it skips a tile already in the file, whatever the file on disk now says —
which is right for a hand-driven update and wrong for a DAG output. (Its ``-s``
single-ID mode implements ``check`` and ``remove``; ``add`` is accepted by the
argument validator and then falls through to the ordinary walk, so it is not a
way to add one tile by hand.)

WHICH TILES — AND WHY THE JOB DERIVES THE SET RATHER THAN BEING TOLD IT. The set
is the CAMPAIGN's: every tile both declared in ``tile_list`` and present in the
index, which is exactly the Snakefile's TILES_READY, rebuilt here from the same
two files the Snakefile started from (``--tile-list`` and ``--index-db``, read
through ``build_index.campaign_tiles`` so there is one definition and not two
that can drift). It is derived rather than passed because at DR6 scale the set
is ~20k paths and a shell command reaches ``execve`` as a SINGLE argv entry
capped at 128 KiB by ``MAX_ARG_STRLEN``; the rule's ``input`` is the DAG edge
and its ``params`` carries a fingerprint of that same list, which is the rerun
trigger.

THE TWO SETS ARE THE SAME SET, which is the point of deriving it this way rather
than globbing ``<products_dir>/tiles``: a products root shared with an earlier,
larger tile list would hand the job tiles the fingerprint never saw and no rerun
trigger would notice. A tile in the derived set whose catalogue is missing is a
hard error here, not a skip — under the DAG it cannot happen, since every one of
them is a declared input of this job.
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

    ``patches/`` is a LEGACY KEY IN sp_validation's FILE SCHEMA, kept verbatim
    only so its reader works unchanged (CosmoStat/sp_validation#340 tracks
    removing it); it names nothing in this workflow, which has campaigns and
    tiles and no other unit. This is the one place the literal appears —
    everything else here says campaign.
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
                / f"final_cat-{tile}.fits")
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
                   help="workflow/config/cfis/final_cat.param — the column list")
    p.add_argument("--hdu", type=int, default=1)
    args = p.parse_args()

    cfc = load_create_final_cat()
    param_list = cfc.read_param_file(str(args.param_file), verbose=False)
    if not param_list:
        sys.exit(f"merge_final_cat: no columns read from {args.param_file}")
    # read_data/copy_data read their knobs out of this dict, exactly as
    # create_final_cat.py's own main() builds it.
    params = {"hdu_num": args.hdu, "param_list": param_list, "verbose": False}

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
                         digest, "n_tiles")
    print(f"[merge_final_cat] {todo.describe()} -> {args.output} "
          f"({len(tiles)} tile(s), {len(param_list)} column(s), "
          f"group {group_path})")


if __name__ == "__main__":
    main()
