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

IT REBUILDS THE WHOLE FILE, IT DOES NOT APPEND. ``create_final_cat.py``'s own
``process()`` skips tiles already in the file, which is right for a hand-driven
incremental update (``-s add`` / ``-s remove`` are that tool's job). A DAG rule
wants the opposite: the output must be a pure function of the input set, so that
a no-op rerun is byte-stable and a changed set is visibly a different file.
Appending would make the result depend on the order campaigns were run in, and
would silently keep a tile whose catalogue was later rebuilt. The cost is
reading every tile's catalogue on every run of the rule — real work at DR6 scale
(~20k tiles), which is why this is not a localrule.

BYTE-STABLE ON A NO-OP RERUN: written to a tmp path, compared, moved only if it
differs (the pattern ``persist_exp.py`` and ``clean_exposure.py`` use). Tiles
are visited in sorted ID order so the file is a function of the input set alone.
An unconditional rewrite would move the output's mtime every invocation.

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
import filecmp
import importlib.util
import sys
from pathlib import Path

import h5py

# Same directory; the rule invokes this file by path, so it is sys.path[0].
import build_index

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

    # tmp-then-cmp-then-mv; the tmp never outlives this process.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_name(args.output.name + ".tmp")
    try:
        tmp.unlink(missing_ok=True)          # h5py "a" would reopen a stale one
        with h5py.File(tmp, "w") as hdf5_file:
            group = hdf5_file.create_group(spval_group(args.campaign))
            for tile, path in tiles:
                extracted, dtype = cfc.read_data(str(path), params)
                data = cfc.copy_data(params["param_list"], extracted, dtype)
                group.create_dataset(tile, data=data, dtype=data.dtype)
            # The same attribute create_final_cat.py's print_list() writes, and
            # what sp_validation reads to know how many tiles it is holding.
            hdf5_file.attrs["n_tiles"] = len(tiles)

        if args.output.exists() and filecmp.cmp(tmp, args.output, shallow=False):
            print(f"[merge_final_cat] unchanged: {args.output}")
        else:
            tmp.replace(args.output)         # atomic: same filesystem
            print(f"[merge_final_cat] {len(tiles)} tile(s), "
                  f"{len(param_list)} column(s) -> {args.output} "
                  f"(group {spval_group(args.campaign)})")
    finally:
        tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
