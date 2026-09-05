#!/usr/bin/env python3
"""Reclaim one finished tile's scratch store and leave a tombstone (PRD #848 D5).

Run as the shell of the in-DAG ``clean_tile`` rule, never by hand: the rule's
``input:`` is the tile's ``final_cat`` on the PERSISTENT root, so by the time
this executes the tile has published its final catalogue.
A tile's scratch store has no reader outside that tile — tiles read exposures,
nothing reads another tile's store — so unlike the exposure case there is no
consumer set to close over and no eligibility test to make. Writer, then
cleaner, and the DAG edge is the whole ordering argument.

What it deletes: the tile's whole ``<run_dir>/tiles/<shard>/<tile>/`` directory.
Measured on 186.307 (smk-g4, a finished 34-tile-campaign tile): 1,279,231,196
bytes across 137 inodes (71 regular files, 9 symlinks, 57 directories) —
``output/run_sp_tile_Sx`` 745 MB, ``run_sp_tile_Uz`` 382 MB, ``run_sp_tile_Mc``
46 MB, ``run_sp_tile_Ms`` 39 MB, ``run_sp_tile_Mh_exp`` 11 MB, and a non-``output/``
remainder of 62 inodes totalling 15,726 bytes. Deleting only ``output/`` is NOT
enough: 62 x 23,114 DR6 tiles is 1.43M inodes against a 1M quota, so the inode
bound binds on its own and the directory has to go as a whole.

The four retained paths
------------------------
The retained paths are used by other mechanisms after the tile is complete.
after this tile is done. Nothing is kept for tidiness, and the three that
already exist are required: ``require_survivors`` checks all of them before
anything is deleted, because a silent drift in one of these paths would not show
up as a broken clean — it would show up much later as a campaign that cannot
resume. The fourth, ``cleaned.json``, is this job's own output and is the one
thing here that cannot be pre-checked; it is written first instead (see ORDER
below).

  1. ``cleaned.json`` — this script's own tombstone, and the tile's surviving
     record. It absorbs every manifest verbatim so ``sp report`` can still
     report a reclaimed tile as ``cleaned`` rather than as "not run"; see
     ``run_report.absorb_tombstones``.

  2. ``manifests/tile_vignets.json`` — CLEAN_EXPOSURE'S CURRENCY, not ours.
     Exposure reclamation keys campaign-wide eligibility on exactly this path:
     ``clean_targets()`` tests it for out-of-scope consumers and ``rule
     clean_exposure`` declares it as an input for in-scope ones. Deleting it
     would silently strand every exposure shared with an out-of-scope tile
     (~7-10 tiles read each exposure), and would fail an in-scope
     ``clean_exposure`` job outright on a missing input if the tile finished
     mid-invocation. It is also not stale in any dangerous sense: it attests
     that ``tile_vignets`` succeeded, which stays true forever, and its product
     — the vignette store — was node-local and died with its job, so it never
     lived on scratch to be contradicted.

  3. ``manifests/tile_find_exposures.json`` — the PREPARE invocation's target.
     ``rule prepare_all_tiles`` declares this manifest for EVERY tile in the
     list, and the tile list accumulates across the campaign, so a cleaned tile
     is still demanded by every later ``sp run``. Delete it and the whole
     ``tile_prep`` group (get_images -> uncompress -> find_exposures) reruns per
     cleaned tile per invocation, which does not merely cost jobs: uncompress
     writes 382 MB back into the store this job just emptied.

  4. ``output/run_sp_tile_Fe/find_exposures_runner/output/exp_numbers-*.txt`` —
     the campaign's tile->exposure edge data, and the reason this is a survivor
     rather than a manifest is that NOTHING re-derives it. ``build_index.build()``
     runs at the parse of EVERY compute invocation, over the whole declared tile
     list, and checks this exact path (``build_index.exp_list_path``); a tile
     whose file is gone counts as missing, and the default
     ``SP_MISSING_THRESHOLD`` of 0.0 makes ONE missing tile a fatal parse. So
     deleting it does not degrade the campaign, it stops it: the first `sp run`
     after the first reclaimed tile cannot build a DAG at all.
     Only the one file is kept — ``run_sp_tile_Fe``'s own ``logs/``, ``tmp/``
     and the runner's process log are ordinary residue and go.

The rest of ``run_sp_tile_Fe``'s parents come along because a file cannot
outlive its directories: 10 inodes per cleaned tile in total (the tile dir,
``cleaned.json``, ``manifests/`` + 2 manifests, and the 4-deep Fe path + its
file), i.e. ~231k inodes at DR6 against the 1M scratch quota — versus 3.2M if
nothing were reclaimed and 1.4M if only ``output/`` were.

What is lost: the per-tile audit trail for both tools that
read it — ``sp_tilecost.py`` and ``sp_costmodel.py``. They attribute the fused
``tile_shape`` group job's cost per tile by reading the tile's SExtractor
catalogue (NAXIS2 of the sexcat = the object count, the cost model's independent
variable, ~380 MB and unkeepable; sp_costmodel also reads its ``EPOCH_k``
extensions for the geometric epoch count) and the eight
``tile_ngmix_<k>.benchmark.tsv`` files. The sexcat is gone for good; the
benchmark ROWS are absorbed into the tombstone under ``benchmarks``, because
they are two lines each and they are the measured-memory feed D4 sizes
``mem_mb`` from — exposure.smk moves ``exp_psf``'s benchmark outside
``manifests/`` for exactly this reason. They are absorbed rather than left on
disk (8 more inodes per tile = 185k at DR6 buys nothing a JSON blob does not),
so ``sp_tilecost.py`` will not find them at their old paths: the record
survives, the tool's current reader does not. FOLLOW-UP, deliberately not done
here: teach ``sp_tilecost.py`` to fall back to ``cleaned.json`` for a tile whose
benchmark TSVs are gone. Until it does, per-chunk cost attribution stops at the
first reclaimed tile even though the numbers are still on disk.

Deletion is symlink-safe (``clean_exposure`` gives the general reason). A
finished tile holds nine symlinks in two classes, and
the second is the one that matters:

  * ``exp_forest/<shard>/<exp>/output`` — 7 links into the EXPOSURE stores,
    each shared with 7-10 other tiles. Rebuildable, but only by re-running those
    chains from VOS.
  * ``output/run_sp_tile_Git/get_images_runner/output/CFIS_{image,weight}-*``
    — 2 links into ``/project/def-mjhudson/unions-wl/tiles``, the staged survey
    imaging: 621 GB across 2,536 files, on the BACKED-UP, GROUP-SHARED
    ``/project``, and not this campaign's to lose. get_images RETRIEVE=symlink
    is what puts them there, so every tile in every campaign carries a pair.

Both classes are handed WHOLESALE to ``shutil.rmtree`` — ``exp_forest/`` as a
top-level entry, ``run_sp_tile_Git/`` as one inside ``output/``, which ``prune``
descends only because the Fe survivor lives there. So the safety rests on
rmtree's own semantics (it unlinks a symlinked entry rather than recursing
through it), NOT on ``prune``'s ``is_symlink()`` test, which fires only for a
link that is itself a direct entry of a level prune walks. Both were verified on
the fixture: every link unlinked, no target followed. Whoever edits ``prune``
next should know that the worst case is not a scratch store they could rebuild —
it is a rmtree walking into half a terabyte of shared, backed-up survey data.

Logs are deleted, not absorbed, as in ``clean_exposure``. Here the
duplication is exact: on a finished tile every ``logs/<stage>.json`` is
BYTE-IDENTICAL to the ``manifests/<stage>.json`` beside it (verified across all
16 stage records of 186.307), because a tile with a failed stage has no
final_cat and so is never cleaned.

Order: write the tombstone first, then delete — ``clean_exposure``'s docstring
the crash window.

There is no ``consumers`` field and no consumer-set staleness to detect, because
a tile has no consumers. What reruns this job is the ordinary machinery
(``script_hash``, and ``mtime`` if final_cat is rewritten). A rerun over an
already-pruned tree deletes nothing and, because absorption is ADDITIVE, does
not blank the record either — see ``previous_record``.
"""

import argparse
import csv
import json
import shutil
import time
from pathlib import Path


def survivor_paths(tile_dir: Path, tile: str) -> dict:
    """``{what it is: path}`` for the three PRE-EXISTING survivors.

    Three, not four: ``cleaned.json`` is this job's own output and does not
    exist yet when this is called.

    The Fe path is spelled out rather than globbed and MUST agree with
    ``build_index.exp_list_path`` — that function is what re-reads it at every
    compute parse. It is not imported, because this script runs inside the
    container as a job shell and stays stdlib-only and import-free like
    ``clean_exposure.py``; ``require_survivors`` below turns a drift between the
    two into a loud failure on the first tile instead of a fatal parse later.
    """
    idra, iddec = tile.split(".")
    return {
        "clean_exposure eligibility (clean_targets / rule clean_exposure input)":
            tile_dir / "manifests" / "tile_vignets.json",
        "prepare_all_tiles target (rule prepare_all_tiles input)":
            tile_dir / "manifests" / "tile_find_exposures.json",
        "index build input (build_index.exp_list_path)":
            tile_dir / "output" / "run_sp_tile_Fe" / "find_exposures_runner"
            / "output" / f"exp_numbers-{idra}-{iddec}.txt",
    }


def require_survivors(survivors: dict, tile: str) -> None:
    """Abort before deleting anything if the survivor contract is not met.

    Fatal, not a warning. Each of these is read by a mechanism OUTSIDE this
    tile, and each failure mode is silent at deletion time and loud much later:
    a missing vignets manifest strands shared exposures, a missing Fe manifest
    reruns the prepare chain, a missing exposure list makes the next compute
    parse exit on the missing-tile threshold. Failing here costs one red job in
    a keep-going run (clean_tile is a leaf, so it poisons no cone) and leaves
    the store intact for inspection.
    """
    missing = {k: p for k, p in survivors.items() if not p.exists()}
    if missing:
        lines = [f"[clean_tile] {tile}: refusing to reclaim — "
                 f"{len(missing)} survivor(s) not on disk:"]
        lines += [f"    {p}\n        owned by: {k}" for k, p in missing.items()]
        lines.append("  Nothing was deleted. Either this tile's store is not in "
                     "the state a finished tile should be in, or one of these "
                     "paths has moved and clean_tile.py has not followed it.")
        raise SystemExit("\n".join(lines))


def previous_record(tombstone: Path) -> tuple:
    """``(manifests, benchmarks)`` from an existing tombstone, or two empties.

    THE ABSORPTION IS ADDITIVE, and this is why. A second clean of an
    already-cleaned tile finds only the two surviving manifests on disk, so a
    fresh absorption would overwrite a complete record with a two-entry one and
    the tile's history would be gone — silently, and for good. Found by running
    the fixture clean twice; not hypothetical, since ``script_hash`` reruns this
    job on any edit to this file. So the previous record is the BASE and what is
    on disk is laid over it: a rebuilt stage's fresh manifest still wins, and a
    reclaimed one keeps the only copy that exists.
    """
    if not tombstone.exists():
        return {}, {}
    try:
        old = json.loads(tombstone.read_text())
    except (OSError, json.JSONDecodeError):
        return {}, {}          # unreadable: start over rather than refuse
    return (dict(old.get("manifests") or {}),
            dict(old.get("benchmarks") or {}))


def absorb_manifests(mdir: Path) -> dict:
    """Every ``manifests/*.json`` verbatim, keyed by file stem.

    By GLOB, so it takes whatever is there; ``run_report.py`` re-keys on each
    body's own ``stage`` field. Same shape as ``clean_exposure``'s absorption,
    including its tolerance for a legacy ``<stage>.failed.json``.
    """
    out = {}
    if mdir.is_dir():
        for f in sorted(mdir.glob("*.json")):
            try:
                out[f.stem] = json.loads(f.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                out[f.stem] = {"unreadable": str(exc)}
    return out


def absorb_benchmarks(mdir: Path) -> dict:
    """The ngmix chunks' benchmark rows, keyed by file stem.

    Two lines each (header + values), so this is a few hundred bytes for the
    campaign's only per-chunk record of runtime, RSS and mean load — the feed D4
    sizes ``tile_ngmix``'s ``mem_mb`` from. Numbers are kept as strings: this is
    an archive of what the TSV said, not a re-measurement.
    """
    out = {}
    if mdir.is_dir():
        for f in sorted(mdir.glob("*.benchmark.tsv")):
            try:
                rows = list(csv.DictReader(f.read_text().splitlines(),
                                           delimiter="\t"))
            except OSError as exc:
                out[f.name] = {"unreadable": str(exc)}
                continue
            if rows:
                out[f.name] = dict(rows[0])
    return out


def prune(root: Path, keep: set, removed: list) -> None:
    """Delete everything under ``root`` except ``keep`` and the dirs leading to it.

    A whitelist walk rather than an rmtree-with-exceptions, so adding a survivor
    is one line and can never be half-implemented: a path is kept iff it is a
    survivor, recursed into iff it is a real directory on the way to one, and
    deleted otherwise.

    ``is_symlink()`` is tested BEFORE ``is_dir()`` and before the recursion —
    see the module docstring on the two symlink classes.
    """
    ancestors = {p for k in keep for p in k.parents}
    for entry in sorted(root.iterdir()):
        if entry in keep:
            continue
        if entry in ancestors and not entry.is_symlink():
            prune(entry, keep, removed)
            continue
        if entry.is_symlink():
            entry.unlink()
        elif entry.is_dir():
            shutil.rmtree(entry)          # unlinks nested symlinks, never follows
        else:
            entry.unlink()
        removed.append(str(entry))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tile-dir", required=True, type=Path)
    p.add_argument("--tile", required=True)
    p.add_argument("--tombstone", required=True, type=Path)
    args = p.parse_args()

    survivors = survivor_paths(args.tile_dir, args.tile)
    require_survivors(survivors, args.tile)

    mdir = args.tile_dir / "manifests"
    manifests, benchmarks = previous_record(args.tombstone)
    manifests.update(absorb_manifests(mdir))
    benchmarks.update(absorb_benchmarks(mdir))

    # Tombstone first, complete — then delete (see the module docstring).
    args.tombstone.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.tombstone.with_suffix(".json.tmp")
    tmp.write_text(json.dumps({
        "tile": args.tile,
        "cleaned_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "kept": sorted(str(p) for p in survivors.values()),
        "manifests": manifests,
        "benchmarks": benchmarks,
    }, indent=2, sort_keys=True) + "\n")
    tmp.replace(args.tombstone)   # atomic: no half-written tombstone, ever

    removed: list = []
    prune(args.tile_dir, set(survivors.values()) | {args.tombstone}, removed)
    print(f"[clean_tile] {args.tile}: removed {len(removed)} path(s); kept "
          f"{len(survivors)} survivor(s) + the tombstone")


if __name__ == "__main__":
    main()
