#!/usr/bin/env python3
"""Persist epoch-cut counts and reclaim finished tiles of an IDLE campaign.

Run --dry-run first. No DAG or Snakemake resume is needed. Products must be
outside scratch's run root. Both historical .fits and .hdf5 catalogues are
supported; a missing, empty or ambiguous persisted catalogue is never cleaned.
The parser and schema are shared with the live completeness checks, and
clean_tile owns all keep rules, the dry-run traversal and actual reclamation.
Do not run while jobs or another cleaner are writing this campaign.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import clean_tile
from completeness import write_if_changed
from epoch_cuts import (
    aggregate_chunks,
    read_chunk_logs,
    require_durable_manifest,
    validate_counts,
)


def resolve_config(config_path):
    """Resolve frozen base defaults plus the run's overlay."""
    import yaml
    from run_config import apply_defaults, merge

    base = config_path.parent / "workflow" / "config.yaml"
    config = apply_defaults(
        merge(
            yaml.safe_load(base.read_text()),
            yaml.safe_load(config_path.read_text()),
        )
    )
    return (
        Path(config["outputs"]["run_dir"]),
        Path(config["outputs"]["products_dir"]),
        int(config["ngmix_chunks"]),
    )


def catalogue_path(products_dir, tile):
    """Return one nonempty persisted catalogue, or None if none exists."""
    root = products_dir / "tiles" / tile[:2] / tile
    paths = [root / f"final_cat-{tile}.{ext}" for ext in ("fits", "hdf5")]
    paths = [p for p in paths if p.is_file() and p.stat().st_size > 0]
    if len(paths) > 1:
        raise ValueError(f"ambiguous persisted final_cat for {tile}: {paths}")
    return paths[0] if paths else None


def backfill_manifest(tile_dir, final_cat, n_chunks):
    """Build the durable manifest, or validate and reuse an existing one."""
    durable = final_cat.parent / "tile_make_cat.json"
    if durable.exists():
        return require_durable_manifest(final_cat, tile_dir.name, n_chunks)
    path = tile_dir / "manifests" / "tile_make_cat.json"
    manifest = json.loads(path.read_text())
    if (
        manifest.get("stage") != "tile_make_cat"
        or manifest.get("unit") != tile_dir.name
        or manifest.get("status") != "complete"
    ):
        raise ValueError(f"not a successful make-cat manifest: {path}")
    chunks = {}
    for k in range(1, n_chunks + 1):
        chunk_path = tile_dir / "manifests" / f"tile_ngmix_{k}.json"
        chunk = json.loads(chunk_path.read_text())
        if (
            chunk.get("stage") != "tile_ngmix"
            or chunk.get("unit") != tile_dir.name
            or chunk.get("status") != "complete"
        ):
            raise ValueError(f"not a successful chunk manifest: {chunk_path}")
        if "epoch_cuts" in chunk:
            chunks[str(k)] = validate_counts(chunk["epoch_cuts"])
        else:
            name = f"run_sp_tile_ngmix_Ng{k}u"
            root = tile_dir / "logs" / "modules" / name
            if not root.exists():
                root = tile_dir / "output" / name
            chunks[str(k)] = read_chunk_logs(root)
    manifest["epoch_cuts"] = aggregate_chunks(chunks, n_chunks)
    return manifest


def process_tile(tile_dir, products_dir, n_chunks, dry_run):
    """Validate one tile, report reclaimable bytes and optionally reclaim."""
    tile = tile_dir.name
    final_cat = catalogue_path(products_dir, tile)
    if final_cat is None:
        print(f"[backfill] {tile}: SKIP no persisted final_cat")
        return "skipped", 0, None
    survivors = clean_tile.survivor_paths(tile_dir, tile)
    clean_tile.require_survivors(survivors, tile)
    manifest = backfill_manifest(tile_dir, final_cat, n_chunks)
    tombstone = tile_dir / "cleaned.json"
    freed = clean_tile.prune(
        tile_dir, set(survivors.values()) | {tombstone}, [], dry_run=True
    )
    record = manifest["epoch_cuts"]
    action = "WOULD FREE" if dry_run else "FREE"
    totals_text = json.dumps(record["totals"], sort_keys=True)
    print(
        f"[backfill] {tile}: {action} {freed} bytes; "
        f"{n_chunks} chunks; totals={totals_text}"
    )
    if not dry_run:
        # Publish atomically BEFORE the cleaner can touch scratch. Never run a
        # Snakemake science stage merely to retrofit this small durable record.
        write_if_changed(
            final_cat.parent / "tile_make_cat.json",
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        )
        clean_tile.reclaim(tile_dir, tile, tombstone, final_cat)
    return "ready", freed, record["totals"]


def main(argv=None):
    """Run the backfill CLI; refusals leave the affected tile intact."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--products-dir", type=Path)
    parser.add_argument("--n-chunks", type=int)
    parser.add_argument(
        "--config",
        type=Path,
        help="run overlay beside frozen workflow/config.yaml",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.config:
        run, products, chunks = resolve_config(args.config)
        if run.resolve() != args.run_dir.resolve():
            parser.error("--run-dir disagrees with the run config")
        if (
            args.products_dir
            and args.products_dir.resolve() != products.resolve()
        ):
            parser.error("--products-dir disagrees with the run config")
        if args.n_chunks is not None and args.n_chunks != chunks:
            parser.error("--n-chunks disagrees with the run config")
        args.products_dir, args.n_chunks = products, chunks
    if args.products_dir is None or args.n_chunks is None or args.n_chunks < 1:
        parser.error(
            "provide --config, or both --products-dir and positive --n-chunks"
        )
    run = args.run_dir.resolve()
    products = args.products_dir.resolve()
    if products == run or run in products.parents:
        parser.error("products must be outside the scratch run root")
    if not (run / "tiles").is_dir():
        parser.error("run has no tiles directory")
    totals = {}
    counts = dict(ready=0, skipped=0, refused=0)
    freed = 0
    for tile_dir in sorted((run / "tiles").glob("*/*")):
        if not re.fullmatch(r"\d+\.\d+", tile_dir.name):
            continue
        try:
            if (
                tile_dir.is_symlink()
                or tile_dir.parent.is_symlink()
                or tile_dir.parent.name != tile_dir.name[:2]
            ):
                raise ValueError(f"unexpected tile/shard path: {tile_dir}")
            status, size, record = process_tile(
                tile_dir, products, args.n_chunks, args.dry_run
            )
            counts[status] += 1
            freed += size
            if record:
                for key, value in record.items():
                    totals[key] = totals.get(key, 0) + value
        except (OSError, ValueError, SystemExit) as exc:
            counts["refused"] += 1
            print(
                f"[backfill] {tile_dir.name}: REFUSED {exc}", file=sys.stderr
            )
    action = "DRY RUN" if args.dry_run else "DONE"
    print(
        f"[backfill] {action}: {counts['ready']} ready, "
        f"{counts['skipped']} skipped, {counts['refused']} refused; "
        f"{freed} bytes ({freed / 2**30:.3f} GiB) "
        f"{'would be freed' if args.dry_run else 'reclaimed'}"
    )
    print(f"[backfill] campaign totals={json.dumps(totals, sort_keys=True)}")
    return 1 if counts["refused"] else 0


if __name__ == "__main__":
    sys.exit(main())
