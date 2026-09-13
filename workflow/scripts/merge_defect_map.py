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

WHAT IS AND IS NOT A FUNCTION OF THE INPUT SET, in the same terms
``hdf5_reconcile.py`` sets them. The map's CONTENT is: the same fragments give
the same valid pixels, the same counts and the same sidecar, whether they
arrived at once or one append at a time. Its BYTES are not — reaching a state by
append rather than by rebuild round-trips the map through healsparse's reader
and writer, which can lay the same pixels out in a different number of 2880-byte
FITS blocks (measured: 1,586,880 B rebuilt vs 1,589,760 B appended, identical
``valid_pixels``). That is the trade for not re-reading the campaign. It means
``write_stable`` below can move the map's mtime on a rebuild that changed
nothing — cheap while nothing consumes the map, and the thing to fix (by
rebuilding whenever the append path would rewrite anyway) if something ever
does. The no-op case is unaffected: it compares the PLAN, not the bytes, and
never opens the map at all.

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
# The two hdf5 merges' reconciler. This file cannot use its `plan`/`apply` — a
# union is not a group of independent datasets, so removing an exposure is a
# rebuild here and a `del` there — but "did this source change since we read
# it" is the SAME question, and answering it twice in two ways is how the two
# halves of a campaign's reconciliation drift apart. So the stamp is imported,
# not re-derived, and the vocabulary below (`stamp`, `Plan`, `empty`,
# `describe`, plan-then-apply, untouched-on-a-no-op) is deliberately theirs.
from hdf5_reconcile import stamp


def fragment_path(products_dir: Path, exp: str) -> Path:
    """Where ``exp_defect_map`` wrote this exposure's fragment."""
    return (products_dir / "exp" / exp[:2] / exp / "defect"
            / f"defect-{exp}.hsp")


def sidecar_stamp(path: Path) -> list:
    """``stamp`` as JSON round-trips it: a list, so a read record compares."""
    return list(stamp(path))


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
                     if exp in known and list(known[exp]) != sidecar_stamp(path))
    if changed:
        return Plan([], sorted(have),
                    f"{len(changed)} fragment(s) changed on disk")
    return Plan(sorted(set(have) - set(known)), [], "")


def union_coverage(paths, nside_coverage) -> np.ndarray:
    """The coverage pixels every fragment in ``paths`` touches, together.

    A PRE-PASS, so the accumulator is allocated ONCE. ``target[pixels] = True``
    into a coverage pixel the map has not seen yet makes healsparse GROW its
    sparse array, which copies it; at DR6 scale that array is gigabytes and a
    rebuild discovers coverage pixels all the way through the campaign, so the
    copies dominate everything the "reading fragments dominates" comment below
    models. Seeding the coverage up front turns O(fragments) reallocations of a
    growing array into one allocation of the final one.

    It costs a second read of each fragment's COVERAGE TABLE only —
    ``HealSparseCoverage.read`` never touches the sparse array — which is
    kilobytes against the megabytes the accumulation itself reads.

    MEASURED, on synthetic fragments at the campaign's own resolution (nside
    131072 / coverage 128), 400 fragments discovering 5131 coverage pixels — a
    656 MB accumulator: accumulation 8.4 s unseeded, 7.1 s seeded, with a 0.7 s
    coverage pre-pass. So the reallocations are ~16% of the accumulation here,
    not the dominant term a naive "copy the array once per new coverage pixel"
    reading predicts (healsparse grows the sparse array in blocks). The win
    grows with the accumulator; the pre-pass does not. Both paths produced
    identical maps.

    Only the rebuild path uses it: an append starts from the map on disk, whose
    coverage is already most of the footprint, and reads a handful of fragments.
    """
    mask = None
    for path in paths:
        cov = hsp.HealSparseCoverage.read(str(path))
        if cov.nside_coverage != nside_coverage:
            # accumulate() is the one place that reports a fragment built at the
            # wrong resolution, with the exposure id and what to do about it.
            # Here it is only a seed: give up on it and let that error stand.
            return None
        mask = (cov.coverage_mask.copy() if mask is None
                else mask | cov.coverage_mask)
    return None if mask is None else np.where(mask)[0]


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


def build_record(have: dict, missing: list, nside_coverage: int, nside: int,
                 n_pixels: int, n_coverage: int) -> dict:
    """The sidecar: what is in the map, and what the campaign wanted but lacked.

    Built apart from writing the map because it is not only the map's record.
    Two of its fields describe the CAMPAIGN — how many exposures it has and
    which of them have no fragment — and those can move while the map itself
    cannot: add tiles whose exposures were all reclaimed by a workflow
    predating this rule and there is nothing to append, nothing to rebuild, and
    a sidecar still reporting the previous campaign's counts. The docstring says
    a short map should say so ON DISK; that means the record has to be rewritten
    even when the map is untouched.
    """
    return {
        "campaign_exposures": len(have) + len(missing),
        "nside": nside,
        "nside_coverage": nside_coverage,
        "n_pixels": n_pixels,
        "n_coverage_pixels": n_coverage,
        # What the NEXT invocation reconciles against; sorted so the sidecar is
        # byte-stable for a given campaign state.
        "exposures": {exp: sidecar_stamp(path)
                      for exp, path in sorted(have.items())},
        # Recorded rather than merely printed: a map short of exposures should
        # say so on disk, not only in a job log nobody keeps.
        "exposures_without_fragment": sorted(missing),
    }


def write_sidecar(sidecar: Path, record: dict) -> None:
    tmp = sidecar.with_name(sidecar.name + ".tmp")
    try:
        tmp.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        write_stable(tmp, sidecar)
    finally:
        tmp.unlink(missing_ok=True)


def apply_plan(output: Path, sidecar: Path, plan: Plan, have: dict,
               missing: list, nside_coverage: int, nside: int) -> dict:
    """Carry the plan out on tmp copies, then move both files into place.

    Map and sidecar are moved together at the end, so a crash mid-merge leaves
    the previous PAIR intact rather than a map the record no longer describes.
    """
    if plan.rebuild:
        todo = plan.rebuild
        target = hsp.HealSparseMap.make_empty(
            nside_coverage, nside, np.bool_, bit_packed=True,
            cov_pixels=union_coverage(
                [have[exp] for exp in todo], nside_coverage))
    else:
        target = hsp.HealSparseMap.read(str(output))
        todo = plan.append
    accumulate(target, [(exp, have[exp]) for exp in todo],
               nside_coverage, nside)

    record = build_record(have, missing, nside_coverage, nside,
                          int(target.n_valid),
                          int(target.coverage_mask.sum()))

    map_tmp = output.with_name(output.name + ".tmp")
    try:
        target.write(str(map_tmp), clobber=True)
        write_stable(map_tmp, output)
    finally:
        map_tmp.unlink(missing_ok=True)
    write_sidecar(sidecar, record)
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
        # The MAP is untouched — that is what an empty plan means, and its mtime
        # must not move. The RECORD still can be stale: campaign_exposures and
        # exposures_without_fragment describe the campaign, not the map, so
        # tiles whose exposures all lack fragments change them without changing
        # a single bit of the union. Rewrite it alone when it differs;
        # write_stable drops the tmp when it does not.
        old_record = read_sidecar(args.sidecar)
        record = build_record(have, missing, args.nside_coverage, args.nside,
                              int(old_record.get("n_pixels", 0)),
                              int(old_record.get("n_coverage_pixels", 0)))
        stale = record != old_record
        if stale:
            write_sidecar(args.sidecar, record)
        print(f"[merge_defect_map] unchanged: {args.output} "
              f"({len(have)} exposure(s)"
              f"{'; sidecar refreshed' if stale else ''})")
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
