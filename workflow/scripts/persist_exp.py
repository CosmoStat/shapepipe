#!/usr/bin/env python3
"""Pack one exposure's retained PSF products into a persistent tar and manifest.

Run through the in-DAG ``exp_persist`` rule, never by hand.
The rule copies products from ``run_dir`` to ``products_dir`` so they survive
scratch's 60-day purge; see ``workflow/rules/exposure.smk`` for the dependency
that orders persistence before cleanup. Keeping retention in a separate rule
lets a keep-list edit repack products without rerunning the four-hour PSF chain.

Product names resolve to file-name globs, searched recursively beneath
``<exp-dir>/output/run_sp_exp_SxSePsf/*/output/``. Recursion includes setools'
``rand_split/``, ``new_cat/``, ``plot/`` and ``stat/`` subdirectories.
An optional product with no matches produces a warning; missing
``psf_validation`` fails before publication, while the scratch store remains.
An empty keep list still packs those mandatory star-catalogue inputs.

One uncompressed tar per exposure limits inode use on /project: loose copies
of all candidates cost about 200 files per exposure, or 2 M at DR6 scale,
against a group quota of about 1 M. A tar holds FITS, ``.psf`` and ``.txt``
products together and lets readers open FITS members through ``io.BytesIO``.
FITS compresses poorly, and an uncompressed archive is seekable.
Flat member names support consumer globs; distinct sources sharing a name fail
rather than overwrite. Ownership is zeroed, members are sorted, and source
mtimes are kept for deterministic archives.

The manifest records each member's product, pattern, source, size and sha256.
It is the rule's sole declared output, on the persistent root at
``<products_dir>/exp/<shard>/<exp>/manifests/``, beside the tar's ``psf/`` directory.
Neither archive nor manifest is replaced when its bytes are unchanged, keeping
mtimes stable for downstream rerun checks. The manifest carries no timestamp.

@sc [label:safety] persist-exp-additive-and-always-validation
`persist_exp:` only ever adds: existing tar members absent from the live match
set are carried forward, and live sources replace matching members.
`psf_validation` is packed on every run whether or not the list names it.
A subtractive keep-list rewrite would delete persistent products whose scratch
sources may be gone; omitting validation would deprive `star_cat_merge` of inputs.
Removing retained products requires a deliberate act on `products_dir`, not a
config edit. Enforced by tests/unit/test_persist_exp_props.py.
"""

import argparse
import filecmp
import hashlib
import json
import sys
import tarfile
from pathlib import Path

# The shared PSF run name from config_exp_psfex.ini and config_exp_mccd.ini.
# This script persists only PSF-stage products; test_workflow_run_names.py
# checks agreement with the configs.
RUN_NAME = "run_sp_exp_SxSePsf"

# --- product catalogue ----------------------------------------------------
# Entries are (glob, per-exposure size, retention benefit), ordered by the
# chain: sextractor, setools, psfex, psfex_interp. Sizes cover 40 CCDs, measured
# on a nibi run of 127 exposures and 64 tiles; None means unmeasured.
# --list-products renders this source table for users and config.yaml.
#
# Validation inputs support audit, re-cutting and incremental star-cat merges
# after scratch is purged (~2 MB per exposure). See the module retention
# contract for their mandatory inclusion.
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
        "setools' per-CCD STAT block: star counts, FWHM mode and cuts. The "
        "selection's summary without its catalogue."),
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

# Raw globs cover products absent from the catalogue. Glob metacharacters,
# dots and spaces distinguish them from bare product identifiers: `*.psf`
# and `default.psfex` are globs; `psf_model` is a product name.
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
    """Return the entry unchanged as the product label, including raw globs."""
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
    """Return stable file lists per entry and the entries with no matches.

    Entries are product names or raw globs; resolve() takes either.
    """
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

    # --list-products needs no exposure; require run arguments only below.
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
    # Require validation even when optional PSF products matched: publishing
    # without the merge's inputs would authorize cleanup and lose those stars.
    # See the module retention contract.
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
    # Overlapping patterns may match the same file; record its first pattern.
    # Distinct source paths with the same flat member name are a collision.
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

    # Carry unmatched members forward (see the module retention contract).
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

    # Compare before replacement; finally removes the tmp on ordinary failures
    # so failed attempts do not consume persistent inodes.
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

    # Each member's sha256, read back from the tar as published. The manifest
    # is the DAG edge star_cat_merge waits on: a refit that changes values but
    # no sizes must change it, and a rerun over the same bytes must not.
    with tarfile.open(tar_path) as tf:
        for f in files:
            f["sha256"] = hashlib.file_digest(
                tf.extractfile(f["name"]), "sha256").hexdigest()

    body = {
        "stage": "exp_persist", "level": "exp", "unit": args.exp,
        "status": "complete",
        "tar": str(tar_path),
        "products": entries,
        "patterns": [resolve(e) for e in entries],
        # Always include unmatched optional entries, even when the list is empty.
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
