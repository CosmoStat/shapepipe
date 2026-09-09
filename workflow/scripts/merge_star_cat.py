#!/usr/bin/env python3
"""Concatenate the campaign's per-CCD PSF validation catalogues into ONE full_starcat.

Run as the shell of the campaign-level ``star_cat_merge`` rule, never by hand.

WHAT IT PRODUCES, AND FOR WHOM. ``<products_dir>/full_starcat-0000000.fits``:
every exposure's every CCD's ``validation_psf-<exp>-<ccd>.fits`` row, stacked,
with a ``CCD_NB`` column recording which CCD each row came from. It is the input
to the rho/tau statistics — sp_validation reads exactly this path
(``star_cat_path`` in its ``scripts/calibration/params.py``) and does no merging
of its own. Historically it was ``combine_runs.bash psf`` + a
``merge_starcat_runner`` pass; the workflow emitted neither, so the product set
was short one file. This script is that pass, driven by the DAG instead of by
bash.

IT DOES NOT REIMPLEMENT THE COLUMN LIST. The stacking, the column names and the
CCD_NB parse all live in ``MergeStarCatPSFEX``
(``shapepipe.modules.merge_starcat_package.merge_starcat``), which is what the
old runner called. This script only decides WHICH catalogues that class is
handed, and where the result lands. A column added to the module is a column
added here for free — which is the entire reason for the indirection.

IT READS THE TARS, IT DOES NOT UNPACK THEM. ``exp_persist`` packs each
exposure's keepers into one uncompressed tar on the persistent root
(``<products_dir>/exp/<shard>/<exp>/psf/<exp>.tar``) precisely because inodes,
not bytes, bind on /project. Unpacking ~20k tars × ~40 members to merge them
would materialise ~800k files on the filesystem that design exists to protect,
and then delete them. So members are read into memory
(``tarfile.extractfile(m).read()`` -> ``io.BytesIO``) one at a time and handed
to the merge class as ``[fileobj, member_name]`` pairs. The member NAME is what
the CCD_NB regex parses, which is why the pair carries it; the class takes the
name from the last element of the entry, so a plain ``[path]`` entry behaves
exactly as it always did.

WHICH EXPOSURES — AND WHY THE JOB DERIVES THE SET RATHER THAN BEING TOLD IT.
The set is the CAMPAIGN's: every exposure read by a tile that is both declared
in ``tile_list`` and present in the index, which is the Snakefile's TILES_READY
walked one edge further. This script rebuilds it from the same two files the
Snakefile started from (``--tile-list`` and ``--index-db``, both small, both on
the persistent root, both read through ``build_index.campaign_exposures`` so
there is one query and not two that can drift), and then takes the exposures
whose ``exp_persist`` manifest is on the persistent root.

It is derived rather than passed because at DR6 scale the set is ~20k paths, and
a shell command reaches ``execve`` as a SINGLE argv entry capped at 128 KiB by
``MAX_ARG_STRLEN``. Passing them would be a job that dies before it starts. So
the rule's ``input`` is the DAG EDGE — what must exist before this runs — and
the rule's ``params`` carries a FINGERPRINT of that same list, which is what
makes the merge rerun when the set changes.

THE TWO SETS ARE THE SAME SET, and that equality is the point of deriving it
this way rather than globbing the tree. The rule's input is ``star_cat_inputs()``
(Snakefile): for each exposure of TILES_READY whose PSF products are on the
persistent root, an edge — the ``exp_persist`` manifest for a live exposure, the
TAR for one whose scratch store reclamation already took (that function argues
the asymmetry, which is about not rebuilding a reclaimed chain from VOS).
Nothing at all for an exposure reclaimed before ``exp_persist`` existed, which
left neither and is unrecoverable short of that rebuild. What this script
selects is the same rule stated from the job's side: same tiles, same index,
manifest present — and by the time the job runs, every exposure with an edge has
one. A glob over ``<products_dir>/exp`` would NOT be the same set: it would
sweep in exposures of an earlier, larger tile list sharing the products root,
stacking rows the fingerprint never saw and no rerun trigger would notice.

THE MANIFEST, NOT THE TAR, IS WHAT IT READS FIRST: the manifest records what was
actually packed, pattern by pattern, member by member, with sizes. Selecting
members from it means this script never guesses at tar contents, and an exposure
whose keep list did not include the validation catalogues contributes nothing
visibly rather than silently.

BYTE-STABLE ON A NO-OP RERUN: written to a tmp path, compared, and moved only
if it differs (the pattern ``persist_exp.py`` and ``clean_exposure.py`` use).
An unconditional rewrite would move the output's mtime on every invocation.
Members are visited in sorted (exposure, member) order so the row order is a
function of the input set alone.

PSFEX ONLY, DELIBERATELY. ``PSF_MODEL`` is ``psfex`` in every campaign the
workflow has run; ``MergeStarCatMCCD`` and ``MergeStarCatSetools`` exist beside
it and take the same constructor, so the hook is the one-line class choice in
``merge_class()`` below — an implementation, not a design, away.
"""

import argparse
import filecmp
import io
import json
import logging
import shutil
import sys
import tarfile
import tempfile
from fnmatch import fnmatch
from pathlib import Path

from shapepipe.modules.merge_starcat_package import merge_starcat

# Same directory; the rule invokes this file by path, so it is sys.path[0].
import build_index

# The output name is not ours to choose: sp_validation hardcodes it
# (`star_cat_path = f"{data_dir}/full_starcat-0000000.fits"`), and
# MergeStarCatPSFEX writes exactly this basename into the output dir it is
# given. Kept here as the name this script promises to produce.
OUT_NAME = "full_starcat-0000000.fits"

# The keep-list pattern whose members this merge consumes. The rule refuses to
# exist unless `persist_exp:` contains a pattern matching this shape (the
# Snakefile does that check at parse time), so by the time we get here the
# members are expected to be present.
MEMBER_PATTERN = "validation_psf-*.fits"


def merge_class(psf_model: str):
    """The merge class for this PSF model — the one-line MCCD/setools hook."""
    try:
        return {"psfex": merge_starcat.MergeStarCatPSFEX,
                "mccd": merge_starcat.MergeStarCatMCCD,
                "setools": merge_starcat.MergeStarCatSetools}[psf_model]
    except KeyError:
        sys.exit(f"merge_star_cat: unknown psf_model {psf_model!r}")


def manifests(products_dir: Path, tile_list: Path, index_db: Path) -> list:
    """The campaign's exp_persist manifests that are on disk, in exposure order.

    Not a glob over the products root: see the module docstring on why the set
    is the campaign's and not the filesystem's.
    """
    out = []
    for exp in build_index.campaign_exposures(tile_list, index_db):
        path = (products_dir / "exp" / exp[:2] / exp / "manifests"
                / "exp_persist.json")
        if path.exists():
            out.append(path)
    return out


def entries(manifest_paths: list, pattern: str) -> tuple:
    """``[fileobj, member_name]`` for every matching member, and the tar count.

    One tar is opened at a time and its members are read into memory; the tars
    are never unpacked to disk (see the module docstring). The returned file
    objects are BytesIO, so nothing stays open on the filesystem — at ~50 KB per
    member and ~40 members per exposure this is ~2 MB per exposure held only for
    as long as the merge takes to consume it, but note that the merge class
    holds the whole stack in python lists regardless, which is the real memory
    term the rule's mem_mb is sized against.
    """
    out, n_tars, empty = [], 0, []
    for man_path in manifest_paths:
        man = json.loads(man_path.read_text())
        wanted = sorted(f["name"] for f in man["files"]
                        if fnmatch(f["name"], pattern))
        if not wanted:
            empty.append(man["unit"])
            continue
        tar_path = Path(man["tar"])
        if not tar_path.exists():
            sys.exit(f"merge_star_cat: {man_path} names a tar that is not "
                     f"there: {tar_path}")
        with tarfile.open(tar_path) as tf:
            for name in wanted:
                member = tf.extractfile(name)
                if member is None:
                    sys.exit(f"merge_star_cat: {tar_path} has no member "
                             f"{name}, which its manifest lists")
                out.append([io.BytesIO(member.read()), name])
        n_tars += 1
    return out, n_tars, empty


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--products-dir", required=True, type=Path,
                   help="the persistent root; exp_persist manifests and tars "
                        "are found beneath it")
    p.add_argument("--tile-list", required=True, type=Path,
                   help="the campaign's tile list (config tile_list)")
    p.add_argument("--index-db", required=True, type=Path,
                   help="the campaign's run index (config outputs.index_db)")
    p.add_argument("--output", required=True, type=Path,
                   help=f"the merged catalogue; its basename is {OUT_NAME}")
    p.add_argument("--psf-model", default="psfex")
    p.add_argument("--pattern", default=MEMBER_PATTERN,
                   help="tar-member glob to merge; default %(default)s")
    args = p.parse_args()

    if args.output.name != OUT_NAME:
        # The merge class writes OUT_NAME into a directory it is handed; a
        # differently-named declared output would silently never be produced.
        sys.exit(f"merge_star_cat: --output must be named {OUT_NAME} "
                 f"(got {args.output.name})")

    log = logging.getLogger("merge_star_cat")
    logging.basicConfig(format="[merge_star_cat] %(message)s",
                        level=logging.INFO, stream=sys.stdout)

    manifest_paths = manifests(args.products_dir, args.tile_list, args.index_db)
    file_list, n_tars, empty = entries(manifest_paths, args.pattern)
    if not file_list:
        # Not a no-op: an empty star catalogue would pass every downstream
        # existence check and produce meaningless rho statistics.
        sys.exit(f"merge_star_cat: no member matched {args.pattern!r} in any "
                 f"of {len(manifest_paths)} exp_persist manifest(s) for this "
                 f"campaign — is '{args.pattern}' in the persist_exp keep list?")
    if empty:
        log.info(f"{len(empty)} exposure(s) persisted no {args.pattern}: "
                 f"{', '.join(sorted(empty)[:5])}"
                 f"{' ...' if len(empty) > 5 else ''}")

    # tmp-then-cmp-then-mv. The merge class chooses its own basename inside the
    # directory it is given, so the tmp is a DIRECTORY, not a file path, and it
    # never outlives this process — an orphan on /project is an inode nothing
    # revisits.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(tempfile.mkdtemp(dir=args.output.parent,
                                    prefix=".star_cat_merge."))
    try:
        merge_class(args.psf_model)(file_list, str(tmp_dir), log).process()
        tmp = tmp_dir / OUT_NAME
        if not tmp.exists():
            sys.exit(f"merge_star_cat: the merge wrote no {OUT_NAME}")
        if args.output.exists() and filecmp.cmp(tmp, args.output, shallow=False):
            log.info(f"unchanged: {args.output}")
        else:
            tmp.replace(args.output)          # atomic: same filesystem
            log.info(f"{len(file_list)} catalogue(s) from {n_tars} exposure(s) "
                     f"-> {args.output}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
