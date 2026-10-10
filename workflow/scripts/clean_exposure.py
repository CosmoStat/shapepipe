#!/usr/bin/env python3
"""Reclaim one exposure's products and stage records, leaving a tombstone.

Run through the in-DAG ``clean_exposure`` rule, never by hand. Its dependencies
order cleanup after persistence and all consuming tiles' vignette extraction;
see ``workflow/rules/exposure.smk`` and the Snakefile's ``clean_targets``.
The tombstone records that consumer set; the rule's ``params`` trigger cleanup
again when the set changes.

This script deletes ``output/``, ``manifests/`` and ``logs/``, not the exposure
directory itself. The ``exp_psf`` benchmark beside those directories survives
for resource sizing.

@sc [label:operations] clean-exposure-remove-stage-records
Remove the stage manifests with the products: surviving declared outputs would
make a later tile treat the exposure as built while its vignette inputs are
absent. Missing manifests let Snakemake regenerate the chain when demanded;
finished tiles do not rerun solely because an intermediate is missing.
Remove logs too so their verdicts do not contradict the reclaimed store.
Successful log verdicts duplicate manifests; see ``completeness.py``.

Current ``manifests/*.json`` records are absorbed under ``manifests`` in the
tombstone, keyed by file stem, including ``<stage>.failed.json`` if present.
``run_report.py`` reads these records to report the unit as ``cleaned`` with
its counts and shortfalls. This script does not merge an existing tombstone.

@sc [label:custody] clean-exposure-record-before-delete
Publish the complete tombstone atomically before deleting any tree. A crash
between publication and deletion leaves the record and unreclaimed disk space;
deleting first would risk losing manifests before their record is preserved.

@sc [label:hazard] clean-exposure-no-follow-deletion
Unlink symlink targets, including dangling links, rather than passing them to
``rmtree``. Nested links are unlinked by ``rmtree`` itself so deletion does not
follow them into shared stores.
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

    # is_symlink() first, and OR'd with exists(): exists() follows the link, so
    # a dangling link would otherwise be skipped and survive.
    candidates = (args.exp_dir / "output", mdir, args.exp_dir / "logs")
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
        # See the module's no-follow deletion contract.
        if target.is_symlink():
            target.unlink()
        else:
            shutil.rmtree(target)
        removed.append(str(target))
    print(f"[clean_exposure] {args.exp}: removed {len(removed)} tree(s) after "
          f"{len(consumers)} consuming tile(s)")


if __name__ == "__main__":
    main()
