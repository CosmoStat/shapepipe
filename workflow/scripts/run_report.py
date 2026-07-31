#!/usr/bin/env python3
"""Standalone run report — NOT a DAG node.

A report rule that declares all tiles' outputs as inputs is a descendant of every
job, so any single hard failure (under --keep-going) poisons its cone and the
report never runs — the exact scenario it exists for (finding 8). So this is a
plain script: it disk-scans the tiles/ and exp/ trees against the count-floor
table (completeness.py) and the index, and emits run_report.json. It runs
automatically via the Snakefile's onsuccess/onerror hooks, and by hand any time
(``sp report``), mid-run included.

It distinguishes whole-exposure absence (all runners missing — a real gap or a
deletion/forest bug) from tolerated per-CCD attrition (counts between floor and
expect), so a future deletion bug cannot be silently absorbed (finding 5).

Attrition is counted at file granularity, not unit granularity: each stage's
``products`` block aggregates found-vs-expected file counts per runner across
all passing units, with a per-unit shortfall map where files are missing.
``warn`` runners (psfex_interp) are exempt from *failing* a unit, not from
being *reported* — the P0 validation found 46/760 psfex_interp CCDs missing
(6%) behind a unit-level "zero attrition" claim, which is the failure mode
this granularity exists to prevent.
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from completeness import COMPLETENESS, check_floor  # noqa: E402

# stage -> (level, run_sp_<prefix> dir name)
STAGE_DIR = {
    "tile_get_images":     ("tile", "run_sp_tile_Git"),
    "tile_uncompress":     ("tile", "run_sp_tile_Uz"),
    "tile_find_exposures": ("tile", "run_sp_tile_Fe"),
    "exp_get_images":      ("exp",  "run_sp_exp_Gie"),
    "exp_split":           ("exp",  "run_sp_exp_Sp"),
    "exp_mask":            ("exp",  "run_sp_exp_Ma"),
    "exp_psf":             ("exp",  "run_sp_exp_SxSePsfPi"),
    "tile_merge_headers":  ("tile", "run_sp_tile_Mh_exp"),
    "tile_mask":           ("tile", "run_sp_tile_Ma"),
    "tile_detect":         ("tile", "run_sp_tile_Sx"),
    "tile_vignets":        ("tile", "run_sp_tile_PiViVi"),
    "tile_merge_cats":     ("tile", "run_sp_Ms"),
    "tile_make_cat":       ("tile", "run_sp_Mc"),
}


def classify(stage, run_dir):
    """'absent' | 'under' | 'attrition' | 'complete' for one unit's stage."""
    if not run_dir.is_dir():
        return "absent", []
    ok, details = check_floor(stage, run_dir)
    if not ok:
        # all mandatory runners empty -> whole-unit absence; else below floor
        mand = [d for d in details if not d[4]]
        if mand and all(n == 0 for _, n, *_ in mand):
            return "absent", details
        return "under", details
    # warn runners are exempt from FAILING a unit (check_floor), not from
    # attrition reporting — psfex_interp is precisely the runner that attrits.
    if any(n < expect for _, n, _, expect, _ in details):
        return "attrition", details
    return "complete", details


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--index", required=True, type=Path)
    p.add_argument("--status", default="manual")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    tiles, exps = [], []
    if args.index.exists():
        con = sqlite3.connect(args.index)
        tiles = [r[0] for r in con.execute("SELECT tile_id FROM tiles")]
        exps = [r[0] for r in con.execute("SELECT exp_id FROM exposures")]
        con.close()

    missing_json = args.index.parent / "missing.json"
    missing = json.loads(missing_json.read_text()) if missing_json.exists() else []

    report = {"status": args.status, "n_tiles": len(tiles), "n_exposures": len(exps),
              "missing_tiles": missing, "stages": {}}
    for stage, (level, prefix) in STAGE_DIR.items():
        if stage not in COMPLETENESS:
            continue
        units = tiles if level == "tile" else exps
        sub = "tiles" if level == "tile" else "exp"
        tally = {"complete": 0, "attrition": 0, "under": [], "absent": []}
        # File-level aggregate per runner over PASSING units (complete/attrition);
        # under/absent units are enumerated by name above, not folded in here —
        # mixing them would double-count whole-unit absence as attrition.
        products = {r: {"found": 0, "expect": 0, "by_unit": {}}
                    for r in COMPLETENESS[stage]}
        for u in units:
            run_dir = args.run_dir / sub / u / "output" / prefix
            verdict, details = classify(stage, run_dir)
            if verdict in ("complete", "attrition"):
                tally[verdict] += 1
                for runner, n, _, expect, _ in details:
                    agg = products[runner]
                    agg["found"] += n
                    agg["expect"] += expect
                    if n < expect:
                        agg["by_unit"][u] = expect - n
            else:
                tally[verdict].append(u)
        for agg in products.values():
            if not agg["by_unit"]:
                del agg["by_unit"]
        tally["products"] = products
        report["stages"][stage] = tally

    finals = [t for t in tiles
              if (args.run_dir / "tiles" / t / f"final_cat-{t}.fits").exists()]
    report["final_cats"] = {"present": len(finals), "of": len(tiles),
                            "missing": [t for t in tiles if t not in finals]}

    out = args.out or (args.index.parent / "run_report.json")
    out.write_text(json.dumps(report, indent=2))
    print(f"[run_report] {report['final_cats']['present']}/{len(tiles)} final "
          f"cats -> {out}", file=sys.stderr)
    for stage, tally in report["stages"].items():
        for runner, agg in tally["products"].items():
            if agg["found"] < agg["expect"]:
                pct = 100 * (agg["expect"] - agg["found"]) / agg["expect"]
                print(f"[run_report] attrition {stage}/{runner}: "
                      f"{agg['found']}/{agg['expect']} files (-{pct:.1f}%) "
                      f"across {len(agg['by_unit'])} units", file=sys.stderr)


if __name__ == "__main__":
    main()
