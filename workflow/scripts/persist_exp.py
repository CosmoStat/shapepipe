#!/usr/bin/env python3
"""Copy ONE exposure's keepable PSF products off scratch, and record what went.

Run as the shell of the in-DAG ``exp_persist`` rule, never by hand.

WHY A COPY AND NOT AN EXEMPTION FROM CLEANUP. The obvious alternative — teach
``clean_exposure`` to spare these files — does not work, because reclamation is
not what threatens them. The exposure store lives on ``run_dir``, which is
/scratch: a 60-day purge takes everything there whether or not this workflow
ever cleaned it. ``products_dir`` is /project, backed up and not purged. So the
only way a per-exposure product outlives its campaign is to LEAVE THE
FILESYSTEM, and that is a copy. Reclamation ordering then falls out for free:
``clean_exposure`` takes this rule's manifest as an input, so the store is never
deleted before its keepers have been written elsewhere.

WHY A SEPARATE RULE AND NOT A ``cp`` APPENDED TO ``exp_psf``. The list of what
to keep is a decision that will be revisited — rho statistics want one file
today, a residual study may want three tomorrow — and ``exp_psf`` is four hours
per exposure. The list rides on this rule's ``params``, so editing it makes
snakemake rerun THIS rule (seconds of cp) and leaves the PSF chain alone. Folded
into ``exp_psf``, the same edit would re-derive every PSF model in the campaign.

WHAT IT SEARCHES. ``<exp-dir>/output/run_sp_exp_SxSePsfPi/*/output/`` — the four
module output dirs of the PSF config (sextractor, setools, psfex, psfex_interp)
— RECURSIVELY. The recursion is not laziness: setools does not write flat, it
writes into ``mask/``, ``rand_split/``, ``new_cat/``, ``plot/`` and ``stat/``
beneath its own output dir, so a caller who wrote ``star_split_ratio_80-*.fits``
meaning "the training star sample" would match nothing under a non-recursive
glob. Patterns are therefore plain FILE names and the layout is ours to know,
not the config author's.

ZERO MATCHES FOR ONE PATTERN IS A WARNING, NOT A FAILURE. setools rejects sparse
CCDs (~0.2% attrition, tolerated by exp_psf's own count floor), so per-CCD
counts are not fixed, and a pattern naming an optional diagnostic may legitimately
find nothing. ZERO FILES IN TOTAL IS A FAILURE: it means the store was not what
we think it is, and writing a green manifest over that would let
``clean_exposure`` delete an exposure whose products were never saved.

The destination is FLAT — one ``psf/`` dir per exposure, no module subtree —
because the module a file came from is already in its name and the consumer
(rho/tau statistics) globs the directory. A name collision between two modules
is therefore a hard error rather than a silent overwrite; nothing in the current
config can produce one, and if a future one can we want to hear about it.

The manifest is the rule's ONLY declared output, and it lives on the persistent
root beside the copies (``<products_dir>/exp/<shard>/<exp>/manifests/``), NOT in
the exposure's scratch ``manifests/`` dir which ``clean_exposure`` deletes
wholesale. It is deliberately NOT a ``directory()`` output: what was copied, and
how big each file was, is provenance we want written down, and a directory
output attests only that some directory exists.

It carries no timestamp and is written tmp-then-``cmp``-then-``mv`` (the pattern
``exp_star_cat`` uses), so a rerun that copies the same files leaves the mtime
alone — mtime is a rerun trigger, and an unconditional rewrite would make every
downstream ``clean_exposure`` look out of date once per invocation.
"""

import argparse
import filecmp
import json
import shutil
import sys
from pathlib import Path

# The PSF chain's run dir (RUN_NAME in config_exp_psfex.ini). Hardcoded rather
# than passed: this rule persists the PSF stage's products and nothing else, and
# a knob here would be a knob for "persist some other stage", which is a
# different rule.
RUN_NAME = "run_sp_exp_SxSePsfPi"


def collect(exp_dir: Path, patterns: list) -> tuple:
    """Matched files per pattern, in a stable order, plus the empty patterns."""
    root = exp_dir / "output" / RUN_NAME
    found, empty = {}, []
    for pat in patterns:
        # One glob per module output dir, recursive beneath it (see the module
        # docstring on setools' subdirectories). sorted() over the union keeps
        # the manifest byte-stable across filesystem readdir order.
        hits = sorted({p for mod in sorted(root.glob("*/output"))
                       for p in mod.rglob(pat) if p.is_file()})
        if hits:
            found[pat] = hits
        else:
            empty.append(pat)
    return found, empty


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exp-dir", required=True, type=Path,
                   help="the exposure's scratch store")
    p.add_argument("--exp", required=True)
    p.add_argument("--dest", required=True, type=Path,
                   help="<products_dir>/exp/<shard>/<exp>/psf")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--pattern", action="append", default=[],
                   help="repeatable; a plain file-name glob")
    args = p.parse_args()

    if not args.pattern:
        sys.exit("persist_exp: no --pattern given (config persist_exp is empty)")

    found, empty = collect(args.exp_dir, args.pattern)
    if not found:
        sys.exit(f"persist_exp: {args.exp}: no file matched any of "
                 f"{args.pattern} under {args.exp_dir}/output/{RUN_NAME}")

    args.dest.mkdir(parents=True, exist_ok=True)
    seen, files = {}, []
    for pat, hits in found.items():
        for src in hits:
            dst = args.dest / src.name
            if src.name in seen:
                sys.exit(f"persist_exp: {args.exp}: two source files are both "
                         f"named {src.name} ({seen[src.name]} and {src}); the "
                         f"destination is flat, so this would silently overwrite")
            seen[src.name] = src
            # Skip a byte-identical copy: it is not just an I/O saving, it keeps
            # the destination's mtimes still for anything downstream that reads
            # them.
            if not (dst.exists() and filecmp.cmp(src, dst, shallow=False)):
                tmp = dst.with_name(dst.name + ".tmp")
                shutil.copy2(src, tmp)
                tmp.replace(dst)          # atomic: no half-copied product
            files.append({"name": src.name, "pattern": pat,
                          "src": str(src), "bytes": dst.stat().st_size})

    body = {
        "stage": "exp_persist", "level": "exp", "unit": args.exp,
        "status": "complete",
        "dest": str(args.dest),
        "patterns": list(args.pattern),
        # The warning the docstring argues for: named patterns that matched
        # nothing. Present as a key even when empty, so a reader never has to
        # wonder whether an old manifest predates the field.
        "patterns_unmatched": empty,
        "n_files": len(files),
        "bytes": sum(f["bytes"] for f in files),
        "files": sorted(files, key=lambda f: f["name"]),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.manifest.with_name(args.manifest.name + ".tmp")
    tmp.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
    if args.manifest.exists() and filecmp.cmp(tmp, args.manifest, shallow=False):
        tmp.unlink()                      # unchanged: leave the mtime alone
    else:
        tmp.replace(args.manifest)

    warn = f" ({len(empty)} pattern(s) matched nothing: {empty})" if empty else ""
    print(f"[persist_exp] {args.exp}: {len(files)} file(s), "
          f"{body['bytes'] / 1e6:.1f} MB -> {args.dest}{warn}")


if __name__ == "__main__":
    main()
