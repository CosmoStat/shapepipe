#!/usr/bin/env python3
"""Collect the campaign's per-CCD PSF validation catalogues into ONE hdf5 file.

Run as the shell of the campaign-level ``star_cat_merge`` rule, never by hand.

WHAT IT PRODUCES, AND FOR WHOM. ``<products_dir>/full_starcat_<campaign>.hdf5``:
one dataset per exposure at ``exposures/<exp>``, holding that exposure's every
CCD's ``validation_psf-<exp>-<ccd>.fits`` rows stacked, with a ``CCD_NB`` column
recording which CCD each row came from. It is the input to the rho/tau
statistics. Historically this was ``combine_runs.bash psf`` plus a
``merge_starcat_runner`` pass producing one flat FITS table,
``full_starcat-0000000.fits``, and sp_validation still opens that name today;
its readers move to this hdf5 under CosmoStat/sp_validation#340, the same
migration that retires the ``patches/`` key on the galaxy side.

WHY HDF5, AND WHY ONE DATASET PER EXPOSURE. The campaign's two products should
behave the same way, and one flat table cannot: appending a tile meant
restacking every exposure the campaign had ever seen — ~40 GB of members at DR6
scale to add ~2 MB. Per-exposure datasets make the file RECONCILABLE
(hdf5_reconcile.py carries that argument, and merge_final_cat.py is the same
machinery on the tile side), so an append reads the appended exposures and
nothing else while the file still cannot drift from its inputs. Memory follows:
one exposure at a time, not one campaign.

NATIVE DTYPES. Columns are written as the validation_psf files store them —
float32 stays float32. The FITS writer this replaces widened every float column
to ``1D``, doubling both the file and the peak memory of the job that wrote it,
for no information.

CCD_NB IS AN INTEGER. It is parsed out of the member name
(``validation_psf-<exp>-<ccd>.fits``), where it is always digits, so a string
buys nothing — and an int column costs 4 bytes a row against the 8 a
two-character fixed-width string does.

IT READS THE TARS, IT DOES NOT UNPACK THEM. ``exp_persist`` packs each
exposure's keepers into one uncompressed tar on the persistent root
(``<products_dir>/exp/<shard>/<exp>/psf/<exp>.tar``) precisely because inodes,
not bytes, bind on /project. Unpacking ~20k tars x ~40 members to merge them
would materialise ~800k files on the filesystem that design exists to protect.
Members are read through the archive's own file object — seekable, the tar
being uncompressed by design — so the counting pass costs a header rather than
a member.

THE OPTIONAL COLUMNS ARE A PER-FILE QUESTION. A pix2wcs-converted catalogue has
no MAG/SNR/ACCEPTED where an ordinary one does, and a campaign can hold both.
Deciding once for the merge is wrong in both directions: it either fails on the
first converted file or silently zeroes the real values of every ordinary one.
Each file is asked for its own schema, and only the files that lack a column are
zero-filled.

WHICH EXPOSURES — AND WHY THE JOB DERIVES THE SET RATHER THAN BEING TOLD IT.
The set is the CAMPAIGN's: every exposure read by a tile that is both declared
in ``tile_list`` and present in the index, which is the Snakefile's TILES_READY
walked one edge further. This script rebuilds it from the same two files the
Snakefile started from (``--tile-list`` and ``--index-db``, both small, both on
the persistent root, both read through ``build_index.campaign_exposures`` so
there is one query and not two that can drift), and then takes the exposures
whose ``exp_persist`` manifest is on the persistent root.

It is derived rather than passed because at DR6 scale the set is ~20k paths and
a shell command reaches ``execve`` as a SINGLE argv entry capped at 128 KiB by
``MAX_ARG_STRLEN``. So the rule's ``input`` is the DAG EDGE — what must exist
before this runs — and its ``params`` carries a FINGERPRINT of the same set,
which is what makes the merge rerun when the set changes. A glob over
``<products_dir>/exp`` would NOT be the same set: it would sweep in exposures of
an earlier, larger tile list sharing the products root, stacking rows the
fingerprint never saw and no rerun trigger would notice.

THE MANIFEST, NOT THE TAR, IS WHAT IT READS FIRST: the manifest records what was
actually packed, member by member, with sizes and the product each came from, so
this script never guesses at tar contents.
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

# The members this merge consumes, named as the keep list names them and
# resolved through the same catalogue persist_exp packs by — so the glob has one
# definition and adding a product cannot leave the two disagreeing. They are
# always there to find: persist_exp packs this product for every exposure
# whatever `persist_exp:` says, and fails the pack rather than writing a
# manifest without it.
MEMBER_PRODUCT = persist_exp.ALWAYS
MEMBER_PATTERN = persist_exp.resolve(MEMBER_PRODUCT)

# The group holding the per-exposure datasets. Unlike the galaxy side's
# `patches/`, this name is ours and says what it holds.
GROUP = "exposures"

# The validation_psf table's HDU: what MergeStarCatPSFEX defaulted to and what
# psfex_interp writes — a SExtractor-style file, empty primary, header-carrying
# image extension, then the table.
HDU = 2

# The columns, in the order the FITS full_starcat carried them, which is the
# order every consumer has seen. The optional three are zero-filled per file.
COLUMNS = ("X", "Y", "RA", "DEC",
           "HSM_G1_PSF", "HSM_G2_PSF", "HSM_T_PSF",
           "HSM_G1_STAR", "HSM_G2_STAR", "HSM_T_STAR",
           "HSM_FLAG_PSF", "HSM_FLAG_STAR")
OPTIONAL = ("MAG", "SNR", "ACCEPTED")
CCD_COLUMN = "CCD_NB"
ALL_COLUMNS = COLUMNS + OPTIONAL + (CCD_COLUMN,)


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

    BY PRODUCT NAME, OR FAILING THAT BY FILE NAME. persist_exp records the
    product every member came from and always packs psf_validation, so the name
    is the answer for anything it writes today. The glob is the fallback, and it
    earns its place twice over: a tar packed before the product field existed
    has no label at all, and a keep list written as a raw glob
    (`validation_psf-*.fits` rather than `psf_validation`) labels its members
    with the glob. Neither should make the campaign's star catalogue silently
    empty.
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

    Every tar is checked for existence HERE, so a products root missing a file
    fails before a single row is read rather than an hour in. The tar is also
    the unit's SOURCE for reconciling: its size and mtime are what a later
    invocation compares against to decide whether this exposure changed.
    """
    chosen, empty = [], []
    for exp, man_path in manifest_paths:
        man = json.loads(man_path.read_text())
        if not any(is_member(f) for f in man["files"]):
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

    TWO PASSES over the tar's members, and neither holds the exposure twice:
    the first reads only each member's FITS HEADER — NAXIS2, the row count —
    and the second allocates the columns once at their exact final length and
    fills them slice by slice. Members are visited in sorted name order, so the
    row order is a function of the tar's contents alone.
    """
    with tarfile.open(tar_path) as tf:
        names = sorted(n for n in tf.getnames()
                       if Path(n).match(MEMBER_PATTERN))
        if not names:
            sys.exit(f"merge_star_cat: {tar_path} holds no {MEMBER_PATTERN}")

        # --- pass 1: row counts and dtypes, from headers alone --------------
        counts, dtypes, opt_dtypes, n_total = [], None, {}, 0
        for name in names:
            with fits.open(tf.extractfile(name), memmap=False,
                           ignore_missing_simple=True) as hdul:
                hdu = hdul[HDU]
                counts.append(hdu.header["NAXIS2"])
                # ColDefs.dtype describes the table without reading it. NOTE:
                # it is the RAW storage dtype and ignores TSCAL/TZERO, so a
                # scaled column would be allocated narrower than the values
                # .data returns. Latent, not live: no validation_psf column is
                # scaled. Read the dtype off .data if one ever is.
                cols = hdu.columns.dtype
                if dtypes is None:
                    dtypes = cols
                for col in OPTIONAL:
                    if col not in opt_dtypes and col in (cols.names or ()):
                        opt_dtypes[col] = cols[col]
            n_total += counts[-1]

        fields = [(c, dtypes[c]) for c in COLUMNS]
        # A column no file of this exposure carries still gets a column,
        # zero-filled, in the dtype the positional column X uses.
        fields += [(c, opt_dtypes.get(c, dtypes["X"])) for c in OPTIONAL]
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
                    # THIS file's schema, not the exposure's.
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
                         digest, "n_exposures")
    print(f"[merge_star_cat] {todo.describe()} -> {args.output} "
          f"({len(chosen)} exposure(s), {len(ALL_COLUMNS)} column(s), "
          f"campaign {args.campaign})")


if __name__ == "__main__":
    main()
