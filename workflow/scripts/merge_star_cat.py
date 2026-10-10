#!/usr/bin/env python3
"""Collect the campaign's per-CCD PSF validation catalogues into one HDF5 file.

Run as the shell of the campaign-level ``star_cat_merge`` rule, never by hand.

``<products_dir>/full_starcat_<run>.hdf5``:
one dataset per exposure at ``exposures/<exp>``, holding that exposure's every
CCD's ``validation_psf-<exp>-<ccd>.fits`` rows stacked, with a ``CCD_NB`` column
recording which CCD each row came from. It is the input to the rho/tau
statistics. sp_validation reads it once CosmoStat/sp_validation#340 lands;
until then it opens the flat ``full_starcat-0000000.fits``.

Per-exposure datasets bound memory to one exposure and allow incremental
updates rather than restacking ~40 GB of members at DR6 scale to add ~2 MB.
See ``hdf5_reconcile.py`` for reconciliation semantics, shared with
``merge_final_cat.py``.

Required columns retain their stored dtypes; optional columns use ``OPTIONAL``'s
fixed dtypes. ``CCD_NB`` is an int32 parsed from the member name
(``validation_psf-<exp>-<ccd>.fits``).

Members are read directly from the persistent tar, without unpacking ~20k
archives x ~40 members into ~800k files. See ``persist_exp.py`` for archive
layout and retention policy. The counting pass reads FITS headers only.

@sc [label:schema] star-merge-optional-column-fill
Check MAG/SNR/ACCEPTED availability per member: a campaign can mix ordinary
and pix2wcs-converted catalogues. Zero-fill only members lacking the column,
so converted files remain readable without discarding ordinary files' values.

``build_index.campaign_exposures`` derives exposures from ``--tile-list`` and
``--index-db`` rather than passing ~20k paths in a shell argument (limited to
128 KiB by ``MAX_ARG_STRLEN``). See the Snakefile's ``star_cat_merge`` rule
for dependencies and the membership fingerprint.

@sc [label:coupling] star-merge-campaign-membership
Consider only exposures used by the tile-list/index intersection, then select
those with persistent manifests and matching validation members. Do not glob
the shared products root: unrelated exposures would evade the rule's fingerprint.

Manifests select archives; ``read_exposure`` selects matching member names from
each archive. Member digests are not verified here.
"""

import argparse
import json
import sys
import tarfile
from fnmatch import fnmatch
from pathlib import Path

import numpy as np
from astropy.io import fits

# Same directory; the rule invokes this file by path, so it is sys.path[0].
import build_index
import hdf5_reconcile
import persist_exp

# Use persist_exp's mandatory validation product and its pattern; see
# persist_exp.py for the packing requirement.
MEMBER_PRODUCT = persist_exp.ALWAYS
MEMBER_PATTERN = persist_exp.resolve(MEMBER_PRODUCT)

# Group holding the per-exposure datasets.
GROUP = "exposures"

# psfex_interp writes a SExtractor-style file: empty primary, header-carrying
# image extension, then the validation_psf table.
HDU = 2

# Required columns in consumer order; optional columns follow below.
COLUMNS = ("X", "Y", "RA", "DEC",
           "HSM_G1_PSF", "HSM_G2_PSF", "HSM_T_PSF",
           "HSM_M4_1_PSF", "HSM_M4_2_PSF", "HSM_RHO4_PSF",
           "HSM_G1_STAR", "HSM_G2_STAR", "HSM_T_STAR",
           "HSM_M4_1_STAR", "HSM_M4_2_STAR", "HSM_RHO4_STAR",
           "HSM_FLAG_PSF", "HSM_FLAG_STAR")
# @sc [label:schema] star-merge-optional-dtypes
# Pin optional dtypes even when every member lacks these columns; otherwise
# exposure datasets can have incompatible structured dtypes at concatenation.
OPTIONAL = {"MAG": np.float32, "SNR": np.float32, "ACCEPTED": np.int32}
CCD_COLUMN = "CCD_NB"
ALL_COLUMNS = COLUMNS + tuple(OPTIONAL) + (CCD_COLUMN,)


def ccd_number(member_name: str) -> int:
    """The CCD this member's rows belong to: ``validation_psf-<exp>-<ccd>.fits``.

    Always digits, which is why the column is an int; a member name that does
    not carry one is a tar we do not understand, and saying so beats writing a
    sentinel into the catalogue.
    """
    ccd = member_name.rsplit(".", 1)[0].rsplit("-", 1)[-1]
    if not ccd.isdigit():
        sys.exit(f"merge_star_cat: cannot read a CCD number out of member "
                 f"name {member_name!r}")
    return int(ccd)


def is_member(entry: dict) -> bool:
    """Is this manifest entry one of the members this merge reads?

    Accept the product label or matching filename. Filename matching supports
    manifests with no product field and keep lists expressed as raw globs.
    This check diagnoses labels; ``tars`` uses filenames for merge membership.
    """
    return (entry.get("product") == MEMBER_PRODUCT
            or fnmatch(entry["name"], MEMBER_PATTERN))


def manifests(products_dir: Path, tile_list: Path, index_db: Path) -> list:
    """``(exposure, manifest path)`` for the campaign's packed exposures.

    Not a glob over the products root: see the module docstring on why the set
    is the campaign's and not the filesystem's.
    """
    out = []
    for exp in build_index.campaign_exposures(tile_list, index_db):
        path = (products_dir / "exp" / exp[:2] / exp / "manifests"
                / "exp_persist.json")
        if path.exists():
            out.append((exp, path))
    return out


def tars(manifest_paths: list) -> tuple:
    """``[(exposure, tar path)]`` for the merge, and the exposures with nothing.

    Check each selected tar's existence before reading rows. The tar is the
    reconciliation source; see ``hdf5_reconcile.py`` for source stamps.
    """
    chosen, empty = [], []
    for exp, man_path in manifest_paths:
        man = json.loads(man_path.read_text())
        # Select by filename, as read_exposure does; a product label alone
        # does not guarantee a selectable member in the tar.
        if not any(fnmatch(f["name"], MEMBER_PATTERN) for f in man["files"]):
            if any(is_member(f) for f in man["files"]):
                # Report a product label inconsistent with the member name.
                print(f"[merge_star_cat] {exp}: manifest labels a "
                      f"{MEMBER_PRODUCT} member whose name does not match "
                      f"{MEMBER_PATTERN}; not merging it")
            empty.append(exp)
            continue
        tar_path = Path(man["tar"])
        if not tar_path.exists():
            sys.exit(f"merge_star_cat: {man_path} names a tar that is not "
                     f"there: {tar_path}")
        chosen.append((exp, tar_path))
    return chosen, empty


def read_exposure(exp: str, tar_path: Path) -> np.ndarray:
    """One exposure's every CCD, stacked, as a structured array.

    Count rows from FITS headers, allocate once, then fill member slices in
    sorted name order. Only one exposure is held in memory. Repacking any
    retention product changes the archive source stamp and can refresh this
    exposure even when its validation members are unchanged.
    """
    try:
        tf = tarfile.open(tar_path)
    except tarfile.TarError as exc:
        sys.exit(f"merge_star_cat: cannot read {tar_path}: {exc}. That tar is "
                 f"this exposure's only copy of its PSF products — do not "
                 f"delete it; re-pack the exposure if its scratch store is "
                 f"still there, and treat the exposure as lost if it is not.")
    with tf:
        names = sorted(n for n in tf.getnames()
                       if fnmatch(n, MEMBER_PATTERN))
        if not names:
            sys.exit(f"merge_star_cat: {tar_path} holds no {MEMBER_PATTERN}")

        # --- pass 1: row counts and dtypes, from headers alone --------------
        counts, dtypes, n_total = [], None, 0
        for name in names:
            with fits.open(tf.extractfile(name), memmap=False,
                           ignore_missing_simple=True) as hdul:
                hdu = hdul[HDU]
                counts.append(hdu.header["NAXIS2"])
                # @sc [label:schema] star-merge-unscaled-required-columns
                # Required validation_psf columns must be unscaled:
                # ColDefs.dtype ignores TSCAL/TZERO. Scaled columns require
                # dtypes from .data to avoid narrowing their decoded values.
                if dtypes is None:
                    dtypes = hdu.columns.dtype
            n_total += counts[-1]

        fields = [(c, dtypes[c]) for c in COLUMNS]
        # See OPTIONAL for the cross-exposure dtype contract.
        fields += list(OPTIONAL.items())
        fields += [(CCD_COLUMN, np.int32)]
        data = np.empty(n_total, dtype=np.dtype(fields))

        # --- pass 2: fill ---------------------------------------------------
        at = 0
        for name, n_rows in zip(names, counts):
            with fits.open(tf.extractfile(name), memmap=False,
                           ignore_missing_simple=True) as hdul:
                rows = hdul[HDU].data
                have = set(rows.dtype.names or ())
                sl = slice(at, at + n_rows)
                for col in COLUMNS:
                    data[col][sl] = rows[col]
                for col in OPTIONAL:
                    # See the module's per-member optional-column contract.
                    data[col][sl] = rows[col] if col in have else 0
                data[CCD_COLUMN][sl] = ccd_number(name)
            at += n_rows

    if at != n_total:
        raise ValueError(f"merge_star_cat: {tar_path}: pass 1 counted "
                         f"{n_total} rows, pass 2 filled {at}")
    return data


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--products-dir", required=True, type=Path,
                   help="the persistent root; exp_persist manifests and tars "
                        "are found beneath it")
    p.add_argument("--tile-list", required=True, type=Path,
                   help="the campaign's tile list (config tile_list)")
    p.add_argument("--index-db", required=True, type=Path,
                   help="the campaign's run index (config outputs.index_db)")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--campaign", required=True,
                   help="named in the log; the group name is fixed")
    p.add_argument("--snapshot-json", type=Path, default=None,
                   help="sp run's code snapshot (bin/sp's "
                        "$STATE_DIR/code/snapshot.json); absent outside sp run")
    args = p.parse_args()

    manifest_paths = manifests(args.products_dir, args.tile_list, args.index_db)
    chosen, empty = tars(manifest_paths)
    if not chosen:
        # Not a no-op: an empty star catalogue would pass every downstream
        # existence check and produce meaningless rho statistics.
        sys.exit(f"merge_star_cat: no {MEMBER_PRODUCT} member in any of "
                 f"{len(manifest_paths)} exp_persist manifest(s) for this "
                 f"campaign. persist_exp packs {MEMBER_PRODUCT} for every "
                 f"exposure, so this means the manifests are not what we think "
                 f"they are.")
    if empty:
        print(f"[merge_star_cat] {len(empty)} exposure(s) persisted no "
              f"{MEMBER_PRODUCT}: {', '.join(sorted(empty)[:5])}"
              f"{' ...' if len(empty) > 5 else ''}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    digest = hdf5_reconcile.schema_digest(ALL_COLUMNS)
    todo = hdf5_reconcile.plan(args.output, GROUP, chosen, digest)
    if todo.empty():
        print(f"[merge_star_cat] unchanged: {args.output} "
              f"({len(chosen)} exposure(s))")
        return
    hdf5_reconcile.apply(args.output, GROUP, todo, chosen, read_exposure,
                         digest, "n_exposures",
                         hdf5_reconcile.code_provenance(args.snapshot_json))
    print(f"[merge_star_cat] {todo.describe()} -> {args.output} "
          f"({len(chosen)} exposure(s), {len(ALL_COLUMNS)} column(s), "
          f"campaign {args.campaign})")


if __name__ == "__main__":
    main()
