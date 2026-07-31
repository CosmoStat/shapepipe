#!/usr/bin/env python3
"""Reclaim ONE exposure's store and leave a tombstone (PRD #848 D5, S5).

Run as the shell of the in-DAG ``clean_exposure`` rule, never by hand: the rule's
``input:`` is every consuming tile's ``tile_vignets`` manifest, so by the time
this executes, every campaign tile that reads this exposure has already extracted
its postage stamps. Writer, then readers, then cleaner — DAG-ordered, race-free.

What it deletes: the exposure's whole ``output/`` tree (the bulk store —
run_sp_exp_Gie/Sp/Ma/SxSePsfPi) AND its ``manifests/``. Deleting the manifests is
deliberate and load-bearing, not tidiness:

  * the manifests are the exposure rules' DECLARED outputs. If they survived, a
    tile appended later would find the exposure chain "up to date" and run
    tile_vignets against products that are no longer on disk. With them gone the
    DAG sees the chain as unbuilt and regenerates it — the accepted cost of a
    late append (D5), expressed as ordinary Snakemake bookkeeping rather than as
    a special case.
  * Snakemake only demands a missing intermediate when something downstream of it
    needs to run, so tiles already finished are NOT rerun by their exposures'
    manifests vanishing.

Nothing is lost to the report: each manifest's content is copied verbatim into
the tombstone under ``manifests`` before deletion, so `sp report` can still say
what this exposure produced and why any runner was short.

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

    removed = []
    for target in (args.exp_dir / "output", mdir):
        if target.exists():
            shutil.rmtree(target)
            removed.append(str(target))

    args.tombstone.parent.mkdir(parents=True, exist_ok=True)
    args.tombstone.write_text(json.dumps({
        "exp": args.exp,
        "cleaned_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "consumers": consumers,
        "removed": removed,
        "manifests": manifests,
    }, indent=2) + "\n")
    print(f"[clean_exposure] {args.exp}: removed {len(removed)} tree(s) after "
          f"{len(consumers)} consuming tile(s)")


if __name__ == "__main__":
    main()
