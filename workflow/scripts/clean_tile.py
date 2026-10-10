#!/usr/bin/env python3
"""Prune a finished tile's scratch store and preserve its tombstone and survivors.

Run through the in-DAG ``clean_tile`` rule, never by hand. The rule orders
cleanup after publication of ``final_cat`` on the persistent root; see
``workflow/rules/tile.smk``. Tiles do not read each other's scratch stores,
so tile cleanup needs no consumer-set eligibility check.

Cleanup removes everything under ``<run_dir>/tiles/<shard>/<tile>/`` except
four retained files and their parents. A finished tile measured on a nibi run
held 1.28 GB across 137 inodes, including 62 inodes outside ``output/``.
Deleting only ``output/`` would leave 1.43 M inodes across 23,114 DR6 tiles,
exceeding the 1 M quota; preserving just the survivors needs about 231k.

Retained paths
--------------
The tombstone ``cleaned.json`` absorbs manifests and benchmark rows for later
reporting; see ``run_report.absorb_tombstones``. The other three files must
already exist before cleanup; ``survivor_paths`` defines them:

  * ``manifests/tile_vignets.json`` attests completion for exposure-cleanup
    eligibility; see ``clean_targets`` in the Snakefile and ``clean_exposure``
    in ``workflow/rules/exposure.smk``.
  * ``manifests/tile_find_exposures.json`` satisfies ``prepare_all_tiles`` in
    the Snakefile, preventing preparation from recreating reclaimed products.
  * ``output/run_sp_tile_Fe/find_exposures_runner/output/exp_numbers-*.txt``
    supplies the index's tile-to-exposure edges; see ``build_index.exp_list_path``.
    Removing it trips the default zero-missing-tile threshold on the next build.

Their parent directories survive too: ten inodes per cleaned tile in total.

Cleanup deletes the SExtractor catalogue and benchmark TSV files.
Benchmark rows survive under ``benchmarks`` in the tombstone instead of costing
eight extra inodes per tile (185k at DR6 scale). Object counts and ``EPOCH_k``
extensions from the deleted catalogue do not survive; cost readers needing
those inputs cannot reconstruct them from the tombstone.

Symlink safety is enforced by ``prune``; links into exposure stores and staged
survey imaging must be unlinked without touching their targets. Logs are deleted
rather than absorbed; successful verdicts are already in the manifests (see
``completeness.py``).

The tombstone is published before deletion; see ``clean_exposure.py`` for the
crash-window rationale. A rerun preserves absorbed records through
``previous_record`` even when their source files have been reclaimed.
"""

import argparse
import csv
import json
import shutil
import time
from pathlib import Path


def survivor_paths(tile_dir: Path, tile: str) -> dict:
    """Return ``{purpose: path}`` for the three pre-existing survivors.

    ``cleaned.json`` is excluded because this job creates it.

    @sc [label:coupling] clean-tile-index-survivor-path
    The explicit Fe path must match ``build_index.exp_list_path`` so compute
    parses can read a cleaned tile's exposure edges. It is duplicated to keep
    this job stdlib-only; ``require_survivors`` detects missing paths before
    deletion rather than leaving the next parse unable to build its index.
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
    """Abort before deletion if any required survivor is missing.

    @sc [label:hazard] clean-tile-required-survivors
    All three survivor paths must exist before cleanup. Each supports a later
    campaign operation (see the module docstring); failing here leaves the
    store intact instead of making that operation fail after reclamation.
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

    @sc [label:custody] clean-tile-additive-tombstone
    Use the readable previous record as the base and overlay files still on
    disk. A rerun over a pruned tile finds only two manifests; fresh absorption
    alone would destroy the other stage records and benchmarks. A rebuilt
    stage's current record supersedes its archived one.
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

    Globbing includes ``<stage>.failed.json`` records if present;
    ``run_report.py`` re-keys on each body's own ``stage`` field.
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
    """The ngmix chunks' benchmark rows, keyed by file name.

    Archive the first data row from each TSV, retaining numeric values as
    strings. This preserves runtime, RSS and mean load without extra inodes.
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

    @sc [label:hazard] clean-tile-no-follow-deletion
    Never recurse through symlinks: direct entries are unlinked before the
    directory test, and ``shutil.rmtree`` unlinks nested links. Both are needed
    because the deleted trees contain links to shared exposure stores and
    backed-up survey imaging outside this tile's ownership.
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
