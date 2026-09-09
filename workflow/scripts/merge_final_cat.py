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

IT RECONCILES, IT NEITHER REBUILDS NOR BLINDLY APPENDS. The output must be a
function of the input set — that is what makes the rule's fingerprint mean
something — but reading every tile's catalogue to add one tile is ~800 GB of IO
at DR6 scale for ~35 MB of new data. So the file is brought INTO AGREEMENT with
the campaign instead:

  * a campaign tile with no dataset is read and added;
  * a dataset whose tile is no longer in the campaign is deleted;
  * a dataset whose source catalogue has CHANGED is re-read. Each one records
    its source's size and mtime as attributes, and a mismatch is what "changed"
    means. This is the only reason a finished tile is ever read twice, and it is
    the reason the file cannot drift from its inputs the way an append-only
    tool does;
  * a dataset that agrees with its source is left alone, unread.

An append therefore reads exactly the appended tiles. ``create_final_cat.py``'s
own ``process()`` implements the append-only half of this — it skips a tile
already in the file, whatever the file on disk now says — which is right for a
hand-driven update and wrong for a DAG output; ``-s add`` / ``-s remove``
remain that tool's way to do this by hand.

WHAT IS AND IS NOT A FUNCTION OF THE INPUT SET. The file's CONTENT is: the same
tiles with the same catalogues give the same datasets, the same columns and the
same n_tiles, whether they arrived at once or one campaign at a time. Its BYTE
LAYOUT is not, because hdf5 lays out a group in the order things were added.
That is the trade for not re-reading the campaign, and it is why the no-op case
below compares actions rather than bytes.
UNTOUCHED ON A NO-OP RERUN, which is stronger than byte-stable and cheaper to
establish. Reconciling is planned before anything is written: if the plan is
empty the file is not opened for writing at all, so its mtime cannot move — and
mtime is a rerun trigger, so an unconditional rewrite would make every
invocation look like a change. When the plan is NOT empty the existing file is
copied to a tmp path, changed there and moved into place, so a crash mid-merge
leaves the old catalogue intact rather than a half-written one. The copy is a
fraction of the reading it replaces.

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
import shutil
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


class Plan:
    """What reconciling this campaign into this file requires: three tile lists.

    ``add`` and ``refresh`` are both "read the catalogue and write the dataset";
    they are separate only so the log can say which happened, because a refresh
    means a finished tile's catalogue moved under us and that is worth seeing.
    """

    def __init__(self, add, refresh, remove):
        self.add, self.refresh, self.remove = add, refresh, remove

    def empty(self):
        return not (self.add or self.refresh or self.remove)

    def describe(self):
        return (f"{len(self.add)} added, {len(self.refresh)} refreshed, "
                f"{len(self.remove)} removed")


def stamp(path: Path) -> tuple:
    """The source catalogue's identity, as recorded on its dataset.

    Size and mtime, not a checksum: the file is ~35 MB and the question is
    "did this change since we read it", which mtime answers for a pipeline
    that writes a catalogue once. A campaign that rewrites a final_cat in
    place with identical size and mtime would defeat it, and nothing does.
    """
    st = path.stat()
    return st.st_size, st.st_mtime_ns


def reconcile_plan(output: Path, group_path: str, tiles: list) -> Plan:
    """Compare the file on disk with the campaign, WITHOUT writing anything.

    Opened read-only, so a no-op invocation cannot move the output's mtime.
    """
    if not output.exists():
        return Plan([t for t, _ in tiles], [], [])

    want = {tile: path for tile, path in tiles}
    add, refresh = [], []
    with h5py.File(output, "r") as f:
        have = dict(f[group_path].items()) if group_path in f else {}
        present = set(have)
        for tile, path in tiles:
            if tile not in present:
                add.append(tile)
                continue
            attrs = have[tile].attrs
            if (int(attrs.get("src_bytes", -1)),
                    int(attrs.get("src_mtime_ns", -1))) != stamp(path):
                refresh.append(tile)
    return Plan(add, refresh, sorted(present - set(want)))


def apply_plan(output: Path, group_path: str, plan: Plan, tiles: list,
               cfc, params: dict) -> None:
    """Carry the plan out on a COPY, then move it into place.

    The copy is what makes a crash mid-merge leave the old catalogue intact,
    and it costs a fraction of the reading it replaces — an append that copies
    a 1 GB file to add one 35 MB tile still beats re-reading the campaign.
    """
    paths = dict(tiles)
    tmp = output.with_name(output.name + ".tmp")
    try:
        tmp.unlink(missing_ok=True)
        if output.exists():
            shutil.copy2(output, tmp)
        with h5py.File(tmp, "a") as f:
            group = f[group_path] if group_path in f else f.create_group(group_path)
            for tile in plan.remove:
                del group[tile]
            for tile in plan.refresh:
                del group[tile]
            for tile in plan.add + plan.refresh:
                path = paths[tile]
                extracted, dtype = cfc.read_data(str(path), params)
                data = cfc.copy_data(params["param_list"], extracted, dtype)
                dset = group.create_dataset(tile, data=data, dtype=data.dtype)
                # The dataset's own record of what it was read from; this is
                # what makes a later invocation able to leave it alone.
                dset.attrs["src_bytes"], dset.attrs["src_mtime_ns"] = stamp(path)
            # The same attribute create_final_cat.py's print_list() writes, and
            # what sp_validation reads to know how many tiles it is holding.
            f.attrs["n_tiles"] = len(group)
        tmp.replace(output)              # atomic: same filesystem
    finally:
        tmp.unlink(missing_ok=True)


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
    plan = reconcile_plan(args.output, group_path, tiles)
    if not plan.empty():
        apply_plan(args.output, group_path, plan, tiles, cfc, params)
        print(f"[merge_final_cat] {plan.describe()} -> {args.output} "
              f"({len(tiles)} tile(s), {len(param_list)} column(s), "
              f"group {group_path})")
    else:
        print(f"[merge_final_cat] unchanged: {args.output} "
              f"({len(tiles)} tile(s))")


if __name__ == "__main__":
    main()
