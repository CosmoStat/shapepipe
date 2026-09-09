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

IT READS THE TARS, IT DOES NOT UNPACK THEM, AND IT STREAMS. ``exp_persist``
packs each exposure's keepers into one uncompressed tar on the persistent root
(``<products_dir>/exp/<shard>/<exp>/psf/<exp>.tar``) precisely because inodes,
not bytes, bind on /project. Unpacking ~20k tars × ~40 members to merge them
would materialise ~800k files on the filesystem that design exists to protect,
and then delete them. So members are read out of the tars in memory
(``tarfile.extractfile(m).read()`` -> ``io.BytesIO``) and handed to the merge
class as ``[fileobj, member_name]`` pairs — ONE AT A TIME, lazily, through
``TarMembers`` below, because materialising them all first is ~40 GB at DR6
scale. The member NAME is what the CCD_NB regex parses, which is why the pair
carries it; the class takes the name from the last element of the entry, so a
plain ``[path]`` entry behaves exactly as it always did.

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
import persist_exp

# The output name is not ours to choose: sp_validation hardcodes it
# (`star_cat_path = f"{data_dir}/full_starcat-0000000.fits"`), and
# MergeStarCatPSFEX writes exactly this basename into the output dir it is
# given. Kept here as the name this script promises to produce.
OUT_NAME = "full_starcat-0000000.fits"

# The members this merge consumes, named as the keep list names them and
# resolved through the same catalogue persist_exp packs by — so the glob has one
# definition and adding a product cannot leave the two disagreeing. The rule
# refuses to exist unless `persist_exp:` keeps something of this shape (the
# Snakefile checks at parse time), so the members are expected here.
MEMBER_PRODUCT = "psf_validation"
MEMBER_PATTERN = persist_exp.resolve(MEMBER_PRODUCT)


def merge_class(psf_model: str):
    """The merge class for this PSF model — the one-line MCCD/setools hook.

    Only psfex is exercised: it is what every campaign has run. MCCD reaches the
    tars unchanged (it takes its CCD numbers from the data, and it now reports
    by the entry's name like the others). SETOOLS would need one more thing —
    it passes ``input_file_list[0][0]`` to file_io as a template path, which a
    streamed entry is not — so wiring setools to this path is a change to that
    class, not a change here.
    """
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


def selection(manifest_paths: list, pattern: str) -> tuple:
    """``[(tar path, [member names])]`` for the merge, and the empty exposures.

    Reads the manifests only. Every tar is checked for existence HERE, so a
    products root missing a file fails before a single row is stacked rather
    than an hour in.
    """
    chosen, empty = [], []
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
        chosen.append((tar_path, wanted))
    return chosen, empty


class TarMembers:
    """The merge class's input list, materialised ONE TAR AT A TIME.

    ``MergeStarCatPSFEX`` wants something it can take the length of and iterate
    once, handing it ``[fileobj, name]`` entries; it never indexes and never
    rewinds. So it does not need a list, and a list is the one thing we cannot
    afford: reading every member up front is the whole campaign in memory at
    once — ~2 MB per exposure, so ~40 GB at DR6's ~20k exposures, against a
    rule asking for 16 GB. Read lazily, peak memory is ONE member's bytes plus
    the merge class's own accumulators, which are the real and unavoidable term.

    ``__len__`` comes from the manifests, so the class can log the count before
    a single tar is opened.

    IT IS ITERABLE MORE THAN ONCE, and must be: the merge makes two passes, one
    for row counts from the headers and one to fill. Each ``__iter__`` opens the
    archives afresh, so the second pass sees the same members in the same order.
    """

    def __init__(self, chosen):
        self._chosen = chosen

    def __len__(self):
        return sum(len(names) for _, names in self._chosen)

    def __iter__(self):
        for tar_path, names in self._chosen:
            with tarfile.open(tar_path) as tf:
                for name in names:
                    member = tf.extractfile(name)
                    if member is None:
                        sys.exit(f"merge_star_cat: {tar_path} has no member "
                                 f"{name}, which its manifest lists")
                    # The tar's own file object, not a BytesIO of the whole
                    # member: it is seekable (the archive is uncompressed by
                    # design) and astropy reads through it, so the merge's
                    # first pass costs a header rather than a member. The
                    # object is valid only until the next member is reached,
                    # which is exactly how the merge consumes it.
                    yield [member, name]


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
    chosen, empty = selection(manifest_paths, args.pattern)
    file_list = TarMembers(chosen)
    if not len(file_list):
        # Not a no-op: an empty star catalogue would pass every downstream
        # existence check and produce meaningless rho statistics.
        sys.exit(f"merge_star_cat: no member matched {args.pattern!r} in any "
                 f"of {len(manifest_paths)} exp_persist manifest(s) for this "
                 f"campaign — is '{MEMBER_PRODUCT}' in the persist_exp keep "
                 f"list?")
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
            log.info(f"{len(file_list)} catalogue(s) from {len(chosen)} "
                     f"exposure(s) -> {args.output}")
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
