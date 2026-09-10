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

RETENTION IS ADDITIVE, AND THAT IS A SAFETY PROPERTY. The keep list rides on
the rule's ``params``, so SHRINKING it reruns this script — and a naive rerun
would rewrite the tar without the products that were dropped, deleting them
from the backed-up filesystem because someone edited a config, with the scratch
store they came from usually long gone. An existing tar is therefore a FLOOR:
its members are carried into the new one whatever the current list says, and a
config change can only ever add. Removing a product is a deliberate act on
products_dir, not a config edit.

THE KEEP LIST IS WHAT THE CAMPAIGN KEEPS ON TOP OF THE MERGE'S INPUTS.
``psf_validation`` is packed unconditionally (see ALWAYS below); ``persist_exp:``
is purely optional retention, and an EMPTY one is a coherent instruction — the
tar then holds the star catalogue's inputs and nothing else.

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
# THE STAR CATALOGUE'S INPUTS ARE NOT A USER CHOICE. star_cat_merge stacks
# every CCD's psf_validation into the campaign's full_starcat, so exp_persist
# ALWAYS packs it, whatever `persist_exp:` says. Two reasons, and neither is
# about taste. It is the merged catalogue's PROVENANCE: a full_starcat with no
# per-exposure inputs beside it cannot be audited, re-cut or recomputed after a
# purge. And it is what keeps APPENDING TILES CHEAP: a tile added next month
# brings exposures whose validation catalogues must join the existing stack, and
# if the earlier ones are gone the merge either shrinks or rebuilds their chains
# from VOS. ~2 MB per exposure, so ~40 GB and ~40k inodes at DR6 scale, against
# a group quota of ~1 M inodes — the cost of being able to say where the number
# came from.
ALWAYS = "psf_validation"

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
                   help=f"repeatable; a product name (see --list-products) or "
                        f"a raw file-name glob. {ALWAYS} is packed whether or "
                        f"not it is named — star_cat_merge needs it")
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

    # The merge's input first and always, then whatever the campaign chose to
    # keep on top of it (see ALWAYS). Deduped, so naming it explicitly in
    # persist_exp: is harmless rather than a repeated pattern.
    entries = [ALWAYS] + [e for e in args.pattern if e != ALWAYS]

    for entry in entries:                    # loud, and before any work
        try:
            resolve(entry)
        except KeyError as exc:
            sys.exit(f"persist_exp: {exc.args[0]}")

    found, empty = collect(args.exp_dir, entries)
    # THE MERGE'S INPUTS ARE NOT ALLOWED TO BE MISSING, and this is a harder
    # rule than "something matched". An exposure whose psfex_interp failed but
    # whose PSFEx model landed has a non-empty match set under the default
    # retention list, so it used to get a green manifest — and clean_exposure
    # takes that manifest as its go-ahead and deletes the store, taking the
    # stars with it. There is no recovering them afterwards short of rebuilding
    # the chain from VOS, so a missing psf_validation fails the job here, while
    # the store is still on disk. Retention products that match nothing stay
    # warnings: they are optional by construction.
    if ALWAYS not in found:
        sys.exit(f"persist_exp: {args.exp}: nothing matched {ALWAYS} "
                 f"({resolve(ALWAYS)}) under {args.exp_dir}/output/{RUN_NAME}. "
                 f"That is the star catalogue's input and it is not optional — "
                 f"refusing to write a manifest that would let clean_exposure "
                 f"reclaim this store.")
    if not found:
        sys.exit(f"persist_exp: {args.exp}: no file matched any of "
                 f"{entries} under {args.exp_dir}/output/{RUN_NAME}")

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

    # --- RETENTION IS ADDITIVE: an existing tar is a FLOOR, never a draft ----
    # Shrinking `persist_exp:` used to rerun this rule (the list rides on
    # params, which is the whole point of the rule) and overwrite the tar with
    # a smaller one — deleting products from the BACKED-UP filesystem because
    # someone edited a config. The scratch store they came from is usually gone
    # by then, so nothing could put them back. Whatever is already in the tar
    # therefore stays in it: a config change can only ever ADD.
    #
    # Removing a product is consequently not a config edit. It is a deliberate
    # act on products_dir, and it should look like one.
    carried, prior_products = [], {}
    if tar_path.exists():
        prior = args.manifest
        if prior.exists():
            try:
                prior_products = {f["name"]: f.get("product", "?")
                                  for f in json.loads(prior.read_text())["files"]}
            except (OSError, ValueError, KeyError):
                pass                      # a damaged manifest loses only labels
        try:
            old_read = tarfile.open(tar_path)
        except tarfile.TarError as exc:
            sys.exit(f"persist_exp: {args.exp}: cannot read the existing "
                     f"{tar_path}: {exc}. Refusing to write a new one — the "
                     f"old tar is left exactly as it is, and it may still hold "
                     f"products nothing else has. Move it aside deliberately "
                     f"if you have decided it is lost.")
        with old_read as tf:
            for ti in tf.getmembers():
                if ti.name in seen or not ti.isfile():
                    continue              # a live source supersedes it
                carried.append(ti.name)
                files.append({"name": ti.name,
                              "product": prior_products.get(ti.name, "?"),
                              "pattern": None, "src": None, "bytes": ti.size})

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
        # One pass in sorted member order, taking each member from whichever
        # side has it: a live source on disk, or the existing tar. Members are
        # copied across with their own TarInfo, so a carried member is
        # byte-for-byte what it was and a rerun that changes nothing still
        # produces an identical archive.
        with tarfile.open(tmp, "w", format=tarfile.PAX_FORMAT) as tf:
            # Already proven readable above, where the members were listed.
            old_tar = (tarfile.open(tar_path) if carried else None)
            try:
                for f in files:
                    if f["name"] in seen:
                        tf.add(seen[f["name"]][0], arcname=f["name"],
                               filter=anonymous)
                    else:
                        ti = anonymous(old_tar.getmember(f["name"]))
                        tf.addfile(ti, old_tar.extractfile(f["name"]))
            finally:
                if old_tar is not None:
                    old_tar.close()
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
        "products": entries,
        "patterns": [resolve(e) for e in entries],
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

    warn = (f" ({len(empty)} retention product(s) matched nothing: {empty})"
            if empty else "")
    if carried:
        warn += f" ({len(carried)} member(s) carried from the existing tar)"
    print(f"[persist_exp] {args.exp}: {len(files)} file(s), "
          f"{body['bytes'] / 1e6:.1f} MB -> {tar_path}{warn}")


if __name__ == "__main__":
    main()
