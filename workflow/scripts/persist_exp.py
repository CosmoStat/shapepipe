#!/usr/bin/env python3
"""Pack ONE exposure's keepable PSF products into a tar off scratch, and record what went.

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

The manifest lists every member (name, pattern, source path, bytes), so a reader
knows what the tar holds without opening it.

ONE UNCOMPRESSED TAR PER EXPOSURE, ``<dest>/<exp>.tar``, NOT LOOSE COPIES.
Inodes, not bytes, are what bind on /project: the group quota is ~1 M files,
and loose per-CCD copies are ~200 per exposure with all candidates on — ~25k for
a 64-tile campaign, ~2 M at DR6 scale, against ~7 GB of bytes. A tar collapses
that to one inode per exposure and costs nothing to read: FITS members go
``tarfile.open(t).extractfile(m).read()`` -> ``fits.open(io.BytesIO(...))``,
which is why a tar rather than a multi-HDU FITS bundle (the keep list mixes
FITS, ``.psf`` and ``.txt``; a FITS container could not hold the last two).
Uncompressed because FITS barely compresses and a plain tar is seekable.

Members are FLAT — file name only, no module subtree — because the module a
file came from is already in its name and the consumer globs member names. A
name collision between two modules is therefore a hard error rather than a
silent overwrite; nothing in the current config can produce one, and if a
future one can we want to hear about it.

The tar is written DETERMINISTICALLY (ownership zeroed, members in sorted
order, source mtimes kept), tmp-then-``cmp``-then-``mv``: a rerun over an
unchanged store produces a byte-identical tar and leaves the existing one's
mtime alone.

The manifest is the rule's ONLY declared output, and it lives on the persistent
root beside the tar (``<products_dir>/exp/<shard>/<exp>/manifests/``, beside the tar's ``psf/``), NOT in
the exposure's scratch ``manifests/`` dir which ``clean_exposure`` deletes
wholesale. It is deliberately NOT a ``directory()`` output: what was copied, and
how big each file was, is provenance we want written down, and a directory
output attests only that some directory exists.

It carries no timestamp and is written tmp-then-``cmp``-then-``mv`` (the pattern
``clean_exposure`` uses), so a rerun that packs the same files leaves the mtime
alone — mtime is a rerun trigger, and an unconditional rewrite would make every
downstream ``clean_exposure`` look out of date once per invocation.
"""

import argparse
import filecmp
import json
import sys
import tarfile
from pathlib import Path

# The PSF chain's run dir (RUN_NAME in config_exp_psfex.ini). Hardcoded rather
# than passed: this rule persists the PSF stage's products and nothing else, and
# a knob here would be a knob for "persist some other stage", which is a
# different rule.
RUN_NAME = "run_sp_exp_SxSePsfPi"

# --- the product catalogue (CosmoStat/shapepipe#844) ------------------------
# THE SINGLE SOURCE OF TRUTH for what an exposure can keep. `persist_exp:` in
# config.yaml names PRODUCTS, not globs: `psf_model`, not `*.psf`. The glob is
# an implementation detail of the module that writes the file, and a keep list
# written in globs is a keep list nobody can read — the argument that produced
# #844 and the 2026-09-08 call's request to keep the PSF model, which had to be
# spelled `*.psf` to be said at all.
#
# Each entry is (glob, per-exposure size, what keeping it buys). Sizes are for
# 40 CCDs, measured on smk-m2 (127 exposures, 64 tiles); "?" means not yet
# measured. `persist_exp.py --list-products` renders this table, and
# config.yaml's block is that rendering rather than a second copy of it.
#
# ORDER IS THE ORDER OF THE CHAIN — sextractor, setools, psfex, psfex_interp —
# so the table reads as the pipeline runs.
PRODUCTS = {
    "star_selection": (
        "star_selection-*.fits", 24_500_000,
        "setools' PRE-SPLIT selection. The only file that answers which stars "
        "the selection cuts rejected and why; the split samples have already "
        "lost the rejects."),
    "star_train": (
        "star_split_ratio_80-*.fits", 19_900_000,
        "the 80% TRAINING sample, the stars PSFEx actually fitted. Rows "
        "duplicate star_selection."),
    "star_test": (
        "star_split_ratio_20-*.fits", 7_100_000,
        "the 20% VALIDATION sample — the positions psf_validation's rows "
        "correspond to. Rows duplicate star_selection."),
    "star_stats": (
        "star_stat-*.txt", None,
        "setools' per-CCD STAT block: star counts, stars/deg^2, FWHM mode and "
        "cuts. The selection's summary without its catalogue."),
    "psf_model": (
        "*.psf", 2_800_000,
        "the PSFEx model itself. Keeping it means the PSF can be "
        "re-interpolated at ANY position later without rebuilding the exposure "
        "chain — the single most capability-adding entry here."),
    "psfex_cat": (
        "psfex_cat-*.cat", None,
        "PSFEx's own output catalogue (FITS_LDAC): the per-star FLAGS_PSF and "
        "CHI2_PSF, i.e. WHICH stars outlier rejection clipped. Not recoverable "
        "from anything else — the .psf header keeps only the LOADED/ACCEPTED "
        "counts."),
    "psf_validation": (
        "validation_psf-*.fits", 2_000_000,
        "the psfex_interp validation catalogue, one per CCD: the input to the "
        "rho/tau statistics, and to the star_cat_merge rule that stacks them "
        "into the campaign's full_starcat."),
}

# PSFEx residual/check images and its XML diagnostics are deliberately absent:
# the committed default.psfex sets CHECKIMAGE_TYPE NONE and WRITE_XML N, so
# nothing is emitted to match. They are a config change first, a catalogue
# entry second.

# A raw glob is still accepted, as an escape hatch for a file the catalogue does
# not name yet. The test is syntactic and deliberately cheap: a product name is
# a bare identifier, so anything carrying a glob metacharacter or a dot is a
# glob. That makes `*.psf`, `star_stat-*.txt` and `default.psfex` globs, and
# `psf_model` a name, with no ambiguity a user could stumble into.
_GLOBBY = set("*?[]. ")


def is_glob(entry: str) -> bool:
    """True when this keep-list entry is a raw glob rather than a product name."""
    return any(ch in _GLOBBY for ch in entry)


def resolve(entry: str) -> str:
    """The file-name glob for one keep-list entry, name or raw glob."""
    if is_glob(entry):
        return entry
    try:
        return PRODUCTS[entry][0]
    except KeyError:
        raise KeyError(
            f"unknown persist_exp product {entry!r}; the products are "
            f"{', '.join(PRODUCTS)} (or write a raw glob such as '*.psf')"
        ) from None


def product_of(entry: str) -> str:
    """The NAME to record for an entry — the entry itself for a raw glob."""
    return entry


def render_products() -> str:
    """The catalogue as a table, for --list-products and for config.yaml."""
    width = max(len(n) for n in PRODUCTS)
    lines = [f"{'product'.ljust(width)}  {'glob'.ljust(26)}  size/exposure",
             f"{'-' * width}  {'-' * 26}  -------------"]
    for name, (glob, size, why) in PRODUCTS.items():
        size_s = "unmeasured" if size is None else f"{size / 1e6:.1f} MB"
        lines.append(f"{name.ljust(width)}  {glob.ljust(26)}  {size_s}")
        for i, chunk in enumerate(_wrap(why, 66)):
            lines.append(f"{' ' * width}      {chunk}")
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list:
    out, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    return out


def collect(exp_dir: Path, patterns: list) -> tuple:
    """Matched files per ENTRY, in a stable order, plus the entries that matched
    nothing. Entries are product names or raw globs; resolve() takes either."""
    root = exp_dir / "output" / RUN_NAME
    found, empty = {}, []
    for entry in patterns:
        pat = resolve(entry)
        # One glob per module output dir, recursive beneath it (see the module
        # docstring on setools' subdirectories). sorted() over the union keeps
        # the manifest byte-stable across filesystem readdir order.
        hits = sorted({p for mod in sorted(root.glob("*/output"))
                       for p in mod.rglob(pat) if p.is_file()})
        if hits:
            found[entry] = hits
        else:
            empty.append(entry)
    return found, empty


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exp-dir", type=Path,
                   help="the exposure's scratch store")
    p.add_argument("--exp")
    p.add_argument("--dest", type=Path,
                   help="<products_dir>/exp/<shard>/<exp>/psf; the tar is "
                        "<dest>/<exp>.tar")
    p.add_argument("--manifest", type=Path)
    p.add_argument("--pattern", action="append", default=[],
                   help="repeatable; a product name (see --list-products) or a "
                        "raw file-name glob")
    p.add_argument("--list-products", action="store_true",
                   help="print the product catalogue and exit")
    args = p.parse_args()

    # --list-products is a QUERY, not a run: it answers "what can I keep?" and
    # needs no exposure, so the run arguments are optional at the parser and
    # required here instead.
    if args.list_products:
        print(render_products())
        return
    missing = [f"--{n.replace('_', '-')}" for n in
               ("exp_dir", "exp", "dest", "manifest")
               if getattr(args, n) is None]
    if missing:
        p.error(f"the following arguments are required: {', '.join(missing)}")

    if not args.pattern:
        sys.exit("persist_exp: no --pattern given (config persist_exp is empty)")

    for entry in args.pattern:               # loud, and before any work
        try:
            resolve(entry)
        except KeyError as exc:
            sys.exit(f"persist_exp: {exc.args[0]}")

    found, empty = collect(args.exp_dir, args.pattern)
    if not found:
        sys.exit(f"persist_exp: {args.exp}: no file matched any of "
                 f"{args.pattern} under {args.exp_dir}/output/{RUN_NAME}")

    args.dest.mkdir(parents=True, exist_ok=True)
    tar_path = args.dest / f"{args.exp}.tar"
    # A file matched by TWO patterns is one file, not a collision. Keep lists
    # overlap on purpose — `validation_psf-*.fits` alongside `*.fits` is a
    # perfectly ordinary way to say "the validation catalogues, and everything
    # else FITS while we are here" — and treating the second match as a name
    # clash failed every exposure in the campaign. What must still be fatal is
    # two DIFFERENT paths landing on one flat member name, which would silently
    # overwrite; that is a same-name/different-source test, and the first
    # pattern to match a file is the one recorded for it.
    seen, files = {}, []
    for pat, hits in found.items():
        for src in hits:
            if src.name in seen:
                if seen[src.name][0] == src:
                    continue          # same file, a second matching pattern
                sys.exit(f"persist_exp: {args.exp}: two source files are both "
                         f"named {src.name} ({seen[src.name][0]} and {src}); tar "
                         f"members are flat, so this would silently overwrite")
            seen[src.name] = (src, pat)
            files.append({"name": src.name, "product": pat,
                          "pattern": resolve(pat),
                          "src": str(src), "bytes": src.stat().st_size})
    files.sort(key=lambda f: f["name"])

    def anonymous(ti: tarfile.TarInfo) -> tarfile.TarInfo:
        # Ownership is the one thing that would differ between two writes of
        # the same files from different accounts/nodes; drop it. mtime stays:
        # it is the product's, and it is stable while the store is.
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        return ti

    # tmp-then-cmp-then-mv, and the tmp NEVER outlives a failure: an orphaned
    # .tmp on /project is an inode nothing revisits — the leak this whole tar
    # design exists to avoid, one per failed attempt at DR6 scale.
    tmp = tar_path.with_name(tar_path.name + ".tmp")
    try:
        with tarfile.open(tmp, "w", format=tarfile.PAX_FORMAT) as tf:
            for f in files:
                tf.add(seen[f["name"]][0], arcname=f["name"], filter=anonymous)
        if tar_path.exists() and filecmp.cmp(tmp, tar_path, shallow=False):
            tmp.unlink()                  # unchanged: leave the mtime alone
        else:
            tmp.replace(tar_path)         # atomic: no half-written archive
    finally:
        tmp.unlink(missing_ok=True)

    body = {
        "stage": "exp_persist", "level": "exp", "unit": args.exp,
        "status": "complete",
        "tar": str(tar_path),
        "products": list(args.pattern),
        "patterns": [resolve(e) for e in args.pattern],
        # The warning the docstring argues for: named patterns that matched
        # nothing. Present as a key even when empty, so a reader never has to
        # wonder whether an old manifest predates the field.
        "patterns_unmatched": empty,
        "n_files": len(files),
        "bytes": sum(f["bytes"] for f in files),
        "files": files,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.manifest.with_name(args.manifest.name + ".tmp")
    try:
        tmp.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
        if args.manifest.exists() and filecmp.cmp(tmp, args.manifest, shallow=False):
            tmp.unlink()                  # unchanged: leave the mtime alone
        else:
            tmp.replace(args.manifest)
    finally:
        tmp.unlink(missing_ok=True)

    warn = f" ({len(empty)} pattern(s) matched nothing: {empty})" if empty else ""
    print(f"[persist_exp] {args.exp}: {len(files)} file(s), "
          f"{body['bytes'] / 1e6:.1f} MB -> {tar_path}{warn}")


if __name__ == "__main__":
    main()
