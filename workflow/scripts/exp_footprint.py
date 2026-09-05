#!/usr/bin/env python3
"""Record ONE exposure's per-CCD sky footprint, for the CCDs that have a PSF model.

Run as the shell of the in-DAG ``exp_footprint`` rule, never by hand.

WHAT THIS REPLACES. The v1.x coverage chain answered "which CCDs have a valid
PSF, and where are they on the sky" by downloading ~25k exposure headers from
VOSpace (``header_downloader.py``), scraping a finished campaign's patch
directories for a ``missing_job_32_all.txt`` and SUBTRACTING it from all 40*N
candidates (``ccd_psf_handler.py``), then appending one row per CCD to a shared
``exp_ra_dec.txt``. Every one of those three moves exists because the pre-
Snakemake pipeline kept no per-unit record. The workflow keeps two, both already
on disk, so this rule is a local read of things it produced itself.

THE TWO INPUTS, AND WHY EACH IS THE RIGHT ONE.

1. WHICH CCDS HAVE A PSF -> ``exp_persist.json``'s ``files[].name``, the
   ``validation_psf-<exp>-<ccd>.fits`` members of the persisted tar. This is
   EXACT, not inferred: ``psfex_interp`` returns WITHOUT writing that file on
   NOT_ENOUGH_STARS, BAD_CHI2 and FILE_NOT_FOUND, and writes it on success, so
   the member list IS the valid-PSF set. It is also DURABLE (persistent root,
   survives ``clean_exposure``) and a DECLARED RULE OUTPUT, so the DAG orders
   this rule after it for free. The alternative considered and rejected was
   ``exp_psf.json``'s psfex_interp entry: that is a COUNT, not a list of names,
   and it lives inside the directory ``clean_exposure`` deletes wholesale.

2. WCS -> ``headers-<exp>.npy``, written by ``split_exp`` beside the per-CCD
   images. A length-N ``dtype=object`` array; element ``i`` is
   ``{"WCS": WCS(h), "header": h.tostring()}``, loaded with ``allow_pickle=True``
   exactly as ``merge_headers.py`` does.

THE LOAD-BEARING INVARIANT, PINNED BY tests/unit/test_exp_footprint.py:
array index ``i`` == the CCD index in ``image-<exp>-<i>.fits`` == ``<exp>-<i>``
== ``validation_psf-<exp>-<i>.fits``. ``split_exp`` writes both the image and the
array element from the same ``idx-1`` in one loop, so the alignment is
structural — but it is a cross-component contract (split_exp's numbering vs.
psfex_interp's filenames vs. this record's ids), and nothing else asserts it now
that the old ``summary.get_all_shdus`` test is gone. No index is ever re-derived
here: a CCD's id comes from its position in the array, never from parsing a name.

WHY THE SHAPE COMES FROM THE STORED HEADER AND NOT FROM THE WCS. astropy hands
back the DECOMPRESSED header for a tile-compressed HDU, so this npy carries true
``NAXIS1/2`` and no ``ZIMAGE`` — while the old VOSpace text headers carried the
binary-table ``NAXIS`` and needed ``ZNAXIS1/2``. ``_image_shape`` handles both,
and is imported rather than reimplemented for exactly that reason. ``WCS`` drops
the ``Z*`` keywords, so it cannot answer the question either way.

DATA AND MANIFEST IN ONE FILE. 40 rows of 8 floats is a few KB, and inodes are
what bind on /project (persist_exp.py argues the quota arithmetic). Per-exposure
JSON, never an appended shared text file: a single appended file is precisely
what cannot survive 20k parallel jobs, and it is what the v1.x chain used.
``ccds_no_psf`` is carried explicitly so the record is self-describing about
attrition — a reader can tell "this CCD is not in the map" from "this CCD was
never looked at".

Written byte-stable, no timestamp, tmp-then-``cmp``-then-``mv``, exactly as
``persist_exp.py`` does and for the same reason: mtime is a rerun trigger, and an
unconditional rewrite would make ``clean_exposure`` look out of date once per
invocation.
"""

import argparse
import filecmp
import json
import sys
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
from astropy.io.fits import Header

from shapepipe.utilities.ccd_footprint import (
    _ccd_corners,
    _image_shape,
)

# The split stage's run dir (RUN_NAME in config_exp_Sp.ini) and its module
# output dir. Hardcoded rather than passed, for persist_exp.py's reason: this
# rule reads the split stage's headers and nothing else, and a knob here would
# be a knob for "read some other stage".
HEADERS_DIR = "output/run_sp_exp_Sp/split_exp_runner/output"

# The name psfex_interp gives a CCD's PSF model at the validation positions.
# Only the PREFIX is ours to know; the rest of the name is the unit id and the
# CCD index, which is what makes the set a set of indices.
PSF_PREFIX = "validation_psf"


def valid_ccds(manifest, exp):
    """The CCD indices with a PSF model, read off an ``exp_persist`` manifest.

    Raises if the manifest was written by a keep list that never asked for the
    PSF files: an empty answer would then be indistinguishable from an exposure
    that genuinely lost every CCD, and the coverage map would silently lose a
    whole exposure. The Snakefile guards the same precondition at parse time;
    this is the per-unit half, and it is the one that sees the manifest.
    """
    body = json.loads(Path(manifest).read_text())
    patterns = body.get("patterns") or []
    # Does the keep list ask for these files AT ALL? Asked by matching a
    # hypothetical member name rather than by string-matching the pattern, so
    # any spelling that would have packed them counts.
    probe = f"{PSF_PREFIX}-{exp}-0.fits"
    if not any(fnmatch(probe, pat) for pat in patterns):
        sys.exit(
            f"exp_footprint: {exp}: {manifest} was written by a keep list "
            f"({patterns}) that packs no {PSF_PREFIX} files, so it names no "
            f"valid-PSF set. Add a pattern matching {probe} to `persist_exp:` "
            f"in config.yaml and rerun exp_persist.")

    ccds = set()
    for entry in body.get("files") or []:
        name = entry.get("name", "")
        if not name.startswith(f"{PSF_PREFIX}-{exp}-"):
            continue
        stem = name[len(f"{PSF_PREFIX}-{exp}-"):].removesuffix(".fits")
        if not stem.isdigit():
            sys.exit(f"exp_footprint: {exp}: cannot read a CCD index out of "
                     f"tar member {name!r}")
        ccds.add(int(stem))
    return ccds


def footprint(headers, exp, with_psf):
    """One record's ``ccds`` and ``ccds_no_psf``, in array order.

    ``headers`` is the ``headers-<exp>.npy`` array; index ``i`` IS the CCD
    index (see the module docstring). Corners are computed only for the CCDs in
    ``with_psf`` — a CCD with no PSF model contributes nothing to a map that
    counts exposures with a valid PSF, and computing its corners anyway would
    invite a later reader to use them.
    """
    ccds, no_psf = [], []
    for i, entry in enumerate(headers):
        ccd_id = f"{exp}-{i}"
        if i not in with_psf:
            no_psf.append(ccd_id)
            continue
        # The WCS is the pickled object split_exp built; the SHAPE must come
        # from the stored header text (module docstring).
        shape = _image_shape(Header.fromstring(entry["header"]))
        ra, dec = _ccd_corners(entry["WCS"], shape)
        # float(), because _ccd_corners hands back numpy scalars and this record
        # has to be byte-stable: a plain double's repr is, a numpy type's
        # serialisation is json's business rather than ours.
        ccds.append({"id": ccd_id,
                     "ra": [float(x) for x in ra],
                     "dec": [float(x) for x in dec]})
    return ccds, no_psf


def write_stable(path, body):
    """Write JSON tmp-then-``cmp``-then-``mv``; an unchanged body keeps its mtime."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        tmp.write_text(json.dumps(body, indent=2, sort_keys=True) + "\n")
        if path.exists() and filecmp.cmp(tmp, path, shallow=False):
            tmp.unlink()                  # unchanged: leave the mtime alone
        else:
            tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exp-dir", required=True, type=Path,
                   help="the exposure's scratch store")
    p.add_argument("--exp", required=True)
    p.add_argument("--persist-manifest", required=True, type=Path,
                   help="<products_dir>/exp/<shard>/<exp>/manifests/"
                        "exp_persist.json; names the valid-PSF CCDs")
    p.add_argument("--manifest", required=True, type=Path)
    args = p.parse_args()

    npy = args.exp_dir / HEADERS_DIR / f"headers-{args.exp}.npy"
    if not npy.exists():
        sys.exit(f"exp_footprint: {args.exp}: no WCS array at {npy} — the "
                 f"split stage's store is gone or never wrote one "
                 f"(config_exp_Sp.ini's OUTPUT_SUFFIX must include `image`).")
    headers = np.load(npy, allow_pickle=True)

    with_psf = valid_ccds(args.persist_manifest, args.exp)
    # A PSF file for a CCD the split never produced means the two stages
    # disagree about the focal plane — the one failure the id alignment cannot
    # absorb, so it is loud rather than silently dropped.
    stray = sorted(i for i in with_psf if i >= len(headers))
    if stray:
        sys.exit(f"exp_footprint: {args.exp}: {args.persist_manifest} names "
                 f"CCD(s) {stray} but {npy} holds only {len(headers)}")

    ccds, no_psf = footprint(headers, args.exp, with_psf)

    write_stable(args.manifest, {
        "stage": "exp_footprint", "level": "exp", "unit": args.exp,
        "status": "complete",
        "source_manifest": str(args.persist_manifest),
        "n_ccd_headers": len(headers),
        "n_valid_psf": len(ccds),
        "ccds": ccds,
        "ccds_no_psf": no_psf,
    })

    warn = f" ({len(no_psf)} without a PSF model)" if no_psf else ""
    print(f"[exp_footprint] {args.exp}: {len(ccds)}/{len(headers)} CCD "
          f"footprint(s){warn} -> {args.manifest}")


if __name__ == "__main__":
    main()
