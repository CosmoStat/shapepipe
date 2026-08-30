#!/usr/bin/env python3
"""Reclaim ONE exposure's store and leave a tombstone (PRD #848 D5, S5).

Run as the shell of the in-DAG ``clean_exposure`` rule, never by hand: the rule's
``input:`` is every consuming tile's ``tile_vignets`` manifest, so by the time
this executes, every campaign tile that reads this exposure has already extracted
its postage stamps. Writer, then readers, then cleaner — DAG-ordered, race-free.

What it deletes: the exposure's whole ``output/`` tree (the bulk store —
run_sp_exp_Gie/Sp/Ma/SxSePsfPi), its ``manifests/`` and its ``logs/``, and its
star-catalogue link farms (``star_cat_exp``, plus the legacy ``star_cat_tiles``). The farms are
reclaimed for consistency, not for bytes: ``exp_star_cat``'s manifest is deleted
here like every other, so the exposure's chain must read as unbuilt, and 40
symlinks left behind are a farm no rule now owns. The catalogue itself lives in
the run-independent cache, so rebuilding the farm costs a relink and no query.

Deletion is SYMLINK-SAFE: a target that is itself a symlink is ``unlink``ed, not
``rmtree``d. Legacy unit dirs carry ``star_cat_exp`` as a link into the old
shared pool, and an rmtree would recurse through it and delete the shared cache
for every other exposure in the campaign.

Deleting the manifests is deliberate and load-bearing, not tidiness:

  * the manifests are the exposure rules' DECLARED outputs. If they survived, a
    tile appended later would find the exposure chain "up to date" and run
    tile_vignets against products that are no longer on disk. With them gone the
    DAG sees the chain as unbuilt and regenerates it — the accepted cost of a
    late append (D5), expressed as ordinary Snakemake bookkeeping rather than as
    a special case.
  * Snakemake only demands a missing intermediate when something downstream of it
    needs to run, so tiles already finished are NOT rerun by their exposures'
    manifests vanishing.

``logs/`` goes with them, and for the same reason rather than for bytes. Each
log holds the completeness verdict of one stage, written on every run and kept by
snakemake through failures; a log left behind would attest "complete" for a store
that is no longer there, contradicting the unbuilt chain the DAG must now see.
Its content for a successful stage is byte-identical to the manifest beside it,
so absorbing the logs into the tombstone would duplicate what the manifests
already carry — they are deleted, not copied.

Nothing is lost to the report: every ``manifests/*.json`` is copied verbatim into
the tombstone under ``manifests``, and ``run_report.py`` reads a cleaned
exposure's record out of the tombstone — it reports the unit as ``cleaned``,
warn counts and shortfalls intact, instead of "not run".

The absorption is by GLOB, so it takes whatever is in ``manifests/``, keyed by
file stem; the report re-keys on each manifest's own ``stage`` field. A legacy
``<stage>.failed.json`` from the pre-``log:`` convention is therefore carried
through unremarkably — it should never be there (an exposure with a failed stage
has no complete vignets consumer and so is not eligible for cleaning), but it
costs nothing to be right about.

Order matters, and it is the reverse of the obvious one: the tombstone is
written FIRST, complete, and only then is anything deleted. A crash between the
two leaves a tombstone beside a store that is still there — the next invocation
treats the exposure as cleaned and only the disk is lost. Deleting first would
put the crash window where the manifests are already gone and the record that
replaces them was never written, and the report would be blind to that exposure
forever.

The exp_psf benchmark tsv lives beside ``manifests/`` and ``logs/``, not inside
either, so it survives this job — it is the measured-memory feed for resource sizing (D4).

The tombstone records the consumer set it was cleaned against. The rule carries
that same set as a ``params`` value, so when the index grows a new consumer the
tombstone goes stale under the default ``params`` rerun-trigger and the clean job
is rescheduled after the new tile's vignets — the exposure is cleaned once per
consumer set, not once per campaign.
"""

import argparse
import json
import shutil
import time
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exp-dir", required=True, type=Path)
    p.add_argument("--exp", required=True)
    p.add_argument("--tombstone", required=True, type=Path)
    p.add_argument("--consumers", default="",
                   help="comma-separated tile ids this exposure was cleaned against")
    args = p.parse_args()

    consumers = [t for t in args.consumers.split(",") if t]

    # Absorb the manifests before they go: the tombstone becomes the exposure's
    # surviving record.
    manifests = {}
    mdir = args.exp_dir / "manifests"
    if mdir.is_dir():
        for f in sorted(mdir.glob("*.json")):
            try:
                manifests[f.stem] = json.loads(f.read_text())
            except (OSError, json.JSONDecodeError) as exc:
                manifests[f.stem] = {"unreadable": str(exc)}

    # is_symlink() first, and OR'd with exists(): exists() follows the link, so a
    # dangling legacy star_cat_exp would otherwise be skipped and survive.
    candidates = (args.exp_dir / "output", mdir, args.exp_dir / "logs",
                  args.exp_dir / "star_cat_exp", args.exp_dir / "star_cat_tiles")
    targets = [t for t in candidates if t.is_symlink() or t.exists()]

    # Tombstone first, complete — then delete (see the module docstring).
    args.tombstone.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.tombstone.with_suffix(".json.tmp")
    tmp.write_text(json.dumps({
        "exp": args.exp,
        "cleaned_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "consumers": consumers,
        "removed": [str(t) for t in targets],
        "manifests": manifests,
    }, indent=2) + "\n")
    tmp.replace(args.tombstone)   # atomic: no half-written tombstone, ever

    removed = []
    for target in targets:
        # NEVER rmtree a symlink (see the module docstring).
        if target.is_symlink():
            target.unlink()
        else:
            shutil.rmtree(target)
        removed.append(str(target))
    print(f"[clean_exposure] {args.exp}: removed {len(removed)} tree(s) after "
          f"{len(consumers)} consuming tile(s)")


if __name__ == "__main__":
    main()
