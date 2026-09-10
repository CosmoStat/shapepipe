#!/usr/bin/env python3
"""Union the campaign's per-exposure defect fragments into ONE healsparse map.

Run as the shell of the campaign-level ``defect_map_merge`` rule, never by hand.

WHAT IT PRODUCES, AND FOR WHOM. ``<products_dir>/defect_map_<campaign>.hsp``:
a boolean ``HealSparseMap``, ``True`` where any exposure of the campaign flagged
the sky, at the mask ladder's own resolution (nside_sparse 131072 over
nside_coverage 128, ``bit_packed``). It is the pixel-domain half of the masking
the footprint could not otherwise see (CosmoStat/shapepipe#878;
``defect_map_exp.py`` argues why the map has to come from here at all).

IT IS A CAMPAIGN PRODUCT, NOT AN INPUT. Nothing in this workflow consumes it:
``config_tile_Mc.ini``'s ``MASK_EXT_PATHS`` names maps that exist before the run
starts, and this one exists only after it. Adding it to that ladder is a
deliberate, later config edit against a path a campaign has actually produced —
that file's header carries the recipe, and deliberately not the entry.

TRUE MEANS MASKED, matching every other map in the ladder. A position outside
the campaign's coverage reads ``False``, indistinguishable from clean, exactly
as the UNIONS bit maps behave outside theirs; the union's coverage is the
campaign's exposures and nothing else says where that is.

IT RECONCILES, AND THE ASYMMETRY IS THE WHOLE DESIGN. The output must be a
function of the input set — that is what makes the rule's fingerprint mean
anything — but a union is not invertible:

  * an exposure whose fragment is NEW is OR-ed into the map on the spot. This is
    the common case (a campaign grows by appending tiles) and it reads exactly
    the appended exposures;
  * an exposure that LEFT the campaign, or whose fragment CHANGED, forces a
    REBUILD from every fragment, because nothing can un-OR a pixel that two
    exposures both set. Rebuilding is honest about that rather than leaving a
    stale bit nothing would ever notice;
  * a plan with neither leaves the file UNTOUCHED — not rewritten identically,
    untouched, so its mtime cannot move. mtime is a rerun trigger.

Which exposures are already in the map is recorded in a SIDECAR beside it
(``defect_map_<campaign>.json``), each with its fragment's size and mtime — the
same "did this change since we read it" stamp ``merge_final_cat.py`` records on
each hdf5 dataset, and for the same reason. The map itself cannot carry that
record: a healsparse FITS header is no place for twenty thousand exposures. The
sidecar is therefore a DECLARED OUTPUT of the rule alongside the map; losing one
without the other would be a map nobody can reconcile, and snakemake removing
both on a failure is the correct recovery (the next run rebuilds).

MEMORY IS FLAT IN THE NUMBER OF EXPOSURES, which is the reason for the loop
below and not for a list comprehension over ``HealSparseMap.read``. Fragments
are never held together and are never OR-ed as maps: each is read, reduced to
its ``valid_pixels`` (~360k int64, ~3 MB on a real exposure), set into the
accumulator, and dropped. The job's footprint is therefore ONE accumulator plus
ONE fragment, whether the campaign is 127 exposures or 20k. The accumulator is
a function of the campaign's FOOTPRINT, not of its exposure count: a coverage
pixel costs (nside/nside_coverage)^2 / 8 bytes bit-packed — 128 KiB at the
ladder's resolution — so a DR6-scale footprint (~23k coverage pixels, measured
on the UNIONS ugriz maps) is ~3 GB resident and the rule is sized on exactly
that count.

WHICH EXPOSURES — AND WHY THE JOB DERIVES THE SET. The campaign's: every
exposure of every tile both declared in ``tile_list`` and present in the index,
through ``build_index.campaign_exposures`` so there is one definition and not
two that can drift. It is derived rather than passed because at DR6 scale the
set is ~20k paths and a shell command reaches ``execve`` as a single argv entry
capped at 128 KiB; the rule's ``input`` is the DAG edge and its ``params``
carries a fingerprint of the same ids.

A CAMPAIGN EXPOSURE WITH NO FRAGMENT IS SKIPPED, NOT AN ERROR, and that is the
one place this differs from ``merge_final_cat``. Fragments accumulate on the
persistent root across campaigns and survive reclamation, but an exposure
reclaimed by a workflow PREDATING this rule has none and never will without a
rebuild from VOS. Failing would make the map unbuildable for exactly the
campaigns that most want it; the count of skipped exposures is reported and
recorded in the sidecar instead.
"""

import argparse
import filecmp
import json
import sys
from pathlib import Path

import numpy as np

import healsparse as hsp

# Same directory; the rule invokes this file by path, so it is sys.path[0].
import build_index


def fragment_path(products_dir: Path, exp: str) -> Path:
    """Where ``exp_defect_map`` wrote this exposure's fragment."""
    return (products_dir / "exp" / exp[:2] / exp / "defect"
            / f"defect-{exp}.hsp")


def stamp(path: Path) -> list:
    """The fragment's identity, as recorded in the sidecar.

    Size and mtime, not a checksum: the fragment is ~2 MB, it is written
    byte-stably (so a no-op re-rasterization does not move its mtime), and the
    question is only "did this change since we read it".
    """
    st = path.stat()
    return [st.st_size, st.st_mtime_ns]


def fragments(products_dir: Path, tile_list: Path, index_db: Path) -> tuple:
    """``({exp: fragment path}, [exposures with no fragment])``, in exposure order."""
    have, missing = {}, []
    for exp in build_index.campaign_exposures(tile_list, index_db):
        path = fragment_path(products_dir, exp)
        if path.exists():
            have[exp] = path
        else:
            missing.append(exp)
    return have, missing


class Plan:
    """What reconciling this campaign into this map requires.

    ``append`` is the cheap path — OR these fragments into the map on disk.
    ``rebuild`` is the honest one: a union cannot drop a pixel, so a removal or
    a changed fragment means reading every fragment again.
    """

    def __init__(self, append, rebuild, reason):
        self.append, self.rebuild, self.reason = append, rebuild, reason

    def empty(self):
        return not (self.append or self.rebuild)

    def describe(self):
        if self.rebuild:
            return f"rebuilt from {len(self.rebuild)} fragment(s) ({self.reason})"
        return f"{len(self.append)} fragment(s) appended"


def read_sidecar(path: Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


def reconcile_plan(output: Path, sidecar: Path, have: dict) -> Plan:
    """Compare what is on disk with the campaign, WITHOUT writing anything.

    A missing map, or a sidecar that does not describe it, is a rebuild: the two
    are written together and either one alone is not evidence about the other.
    """
    record = read_sidecar(sidecar)
    known = record.get("exposures") or {}
    if not output.exists() or not known:
        return Plan([], sorted(have), "no map on disk")

    gone = sorted(set(known) - set(have))
    if gone:
        return Plan([], sorted(have),
                    f"{len(gone)} exposure(s) left the campaign")
    changed = sorted(exp for exp, path in have.items()
                     if exp in known and list(known[exp]) != stamp(path))
    if changed:
        return Plan([], sorted(have),
                    f"{len(changed)} fragment(s) changed on disk")
    return Plan(sorted(set(have) - set(known)), [], "")


def accumulate(target, paths, nside_coverage, nside) -> None:
    """OR each fragment into ``target``, ONE AT A TIME (see the docstring).

    A fragment at the wrong resolution is a hard error rather than a silent
    upgrade: the ladder's nside is a campaign-wide decision, and a fragment that
    disagrees with it was rasterized by a differently-configured run.
    """
    for exp, path in paths:
        frag = hsp.HealSparseMap.read(str(path))
        if (frag.nside_sparse, frag.nside_coverage) != (nside, nside_coverage):
            sys.exit(f"merge_defect_map: {path} is nside_sparse "
                     f"{frag.nside_sparse} / nside_coverage "
                     f"{frag.nside_coverage}, not {nside} / {nside_coverage}; "
                     f"re-rasterize {exp} before merging")
        pixels = frag.valid_pixels
        del frag
        if pixels.size:
            target[pixels] = True
        del pixels


def write_stable(tmp: Path, dest: Path) -> None:
    if dest.exists() and filecmp.cmp(tmp, dest, shallow=False):
        tmp.unlink()
    else:
        tmp.replace(dest)


def apply_plan(output: Path, sidecar: Path, plan: Plan, have: dict,
               missing: list, nside_coverage: int, nside: int) -> dict:
    """Carry the plan out on tmp copies, then move both files into place.

    Map and sidecar are moved together at the end, so a crash mid-merge leaves
    the previous PAIR intact rather than a map the record no longer describes.
    """
    if plan.rebuild:
        target = hsp.HealSparseMap.make_empty(
            nside_coverage, nside, np.bool_, bit_packed=True)
        todo = plan.rebuild
    else:
        target = hsp.HealSparseMap.read(str(output))
        todo = plan.append
    accumulate(target, [(exp, have[exp]) for exp in todo],
               nside_coverage, nside)

    record = {
        "campaign_exposures": len(have) + len(missing),
        "nside": nside,
        "nside_coverage": nside_coverage,
        "n_pixels": int(target.n_valid),
        "n_coverage_pixels": int(target.coverage_mask.sum()),
        # What the NEXT invocation reconciles against; sorted so the sidecar is
        # byte-stable for a given campaign state.
        "exposures": {exp: stamp(path) for exp, path in sorted(have.items())},
        # Recorded rather than merely printed: a map short of exposures should
        # say so on disk, not only in a job log nobody keeps.
        "exposures_without_fragment": missing,
    }

    map_tmp = output.with_name(output.name + ".tmp")
    side_tmp = sidecar.with_name(sidecar.name + ".tmp")
    try:
        target.write(str(map_tmp), clobber=True)
        side_tmp.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        write_stable(map_tmp, output)
        write_stable(side_tmp, sidecar)
    finally:
        map_tmp.unlink(missing_ok=True)
        side_tmp.unlink(missing_ok=True)
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--products-dir", required=True, type=Path,
                        help="the persistent root; fragments are found beneath it")
    parser.add_argument("--tile-list", required=True, type=Path,
                        help="the campaign's tile list (config tile_list)")
    parser.add_argument("--index-db", required=True, type=Path,
                        help="the campaign's run index (config outputs.index_db)")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sidecar", required=True, type=Path,
                        help="the reconciliation record written beside the map")
    parser.add_argument("--nside", type=int, default=131072)
    parser.add_argument("--nside-coverage", type=int, default=128)
    args = parser.parse_args()

    have, missing = fragments(args.products_dir, args.tile_list, args.index_db)
    if not have:
        # An empty map would satisfy every downstream existence check and mask
        # nothing anywhere.
        sys.exit(f"merge_defect_map: no campaign exposure in {args.tile_list} "
                 f"has a defect fragment under {args.products_dir}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plan = reconcile_plan(args.output, args.sidecar, have)
    if plan.empty():
        print(f"[merge_defect_map] unchanged: {args.output} "
              f"({len(have)} exposure(s))")
        return
    record = apply_plan(args.output, args.sidecar, plan, have, missing,
                        args.nside_coverage, args.nside)
    warn = (f"; {len(missing)} campaign exposure(s) have no fragment"
            if missing else "")
    print(f"[merge_defect_map] {plan.describe()} -> {args.output} "
          f"({record['n_pixels']} healpix pixel(s) over "
          f"{record['n_coverage_pixels']} coverage pixel(s)){warn}")


if __name__ == "__main__":
    main()
