#!/usr/bin/env python3
"""``sp report`` — the run's success/failure tables, read from the manifests.

This is a plain script, not a DAG node: depending on all tile outputs would
prevent reporting after a hard failure under ``--keep-going``. It can run
mid-run; the Snakefile's onsuccess/onerror hooks call it at invocation end.

It reads declared units and recorded verdicts:

  * the **index** (``run_index.sqlite``) — the units the run declared, and the
    tile->exposure edges that let an exposure failure be blamed on the tiles it
    blocks;
  * the **verdicts** written by ``completeness.py check`` — per-runner
    found/expect and log-scraped failure reasons — in the unit's
    ``logs/`` (every run) and ``manifests/`` (successes only; completeness.py
    defines the record placement).

Both directories are read and status comes from the file body. The manifest
glob also accepts ``<stage>.failed.json`` records. No record means "not run"
for reporting purposes; it is not evidence that no process executed.

Reclaimed units' records come from ``cleaned.json`` tombstones. See
``clean_exposure.py`` and ``clean_tile.py`` for reclamation and surviving
manifests, and ``absorb_tombstones`` for the generation boundary. Cleaned
exposures block no tile.

No disk scanning: counting products is the *check's* job, done once at the moment
the products were fresh. A unit with no manifest for a stage is "not run" — which
is a real and distinct answer from "ran and produced nothing".

Records are discovered by glob (``tiles/*/*/manifests/*.json`` and
``tiles/*/*/logs/*.json``), not by constructed path: units are found rather than
named, so a store holding a unit the index never heard of still reports. The
depth is fixed at two, matching ``tiles/<prefix>/<ID>/`` and avoiding a
recursive walk of unreclaimed ``output/`` trees. With an index, tallies cover
indexed units; without one, they cover units discovered from records.
"""

import argparse
import json
import os
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

# Stage order per level — the report's column order, and the definition of
# "expected" (a declared unit with no manifest for one of these is not run).
TILE_STAGES = ["tile_get_images", "tile_uncompress", "tile_find_exposures",
               "tile_merge_headers", "tile_detect", "tile_vignets",
               "tile_ngmix", "tile_merge_cats", "tile_make_cat"]
# The catalogue fetch exists only when the run joins the UNIONS catalogue
# (config.yaml's tile_detection); the callers pass the mode in the environment
# so a SExtractor run does not report the stage as not run.
if os.environ.get("SP_TILE_DETECTION") == "unions_catalogue":
    TILE_STAGES.insert(TILE_STAGES.index("tile_detect"), "tile_get_catalogue")
EXP_STAGES = ["exp_get_images", "exp_split", "exp_psf"]

# exp_persist is excluded: its manifest lives on products_dir, while this
# report reads only run_dir. Listing it would mark every exposure "not run".

# See clean_tile.py for ownership of these surviving manifests. Their
# presence alone is not evidence of a rebuilt chain (absorb_tombstones).
SURVIVING_TILE_STAGES = frozenset({"tile_vignets", "tile_find_exposures"})

STATUSES = ("complete", "warn", "failed", "not_run")


def _rank(m: dict) -> int:
    """Severity as a position in STATUSES; unknown statuses sort worst."""
    return STATUSES.index(m["status"]) if m.get("status") in STATUSES else len(STATUSES)


def keep_worst(records: dict, stage: str, m: dict) -> None:
    """Collapse a stage's records to the one with the worst status.

    @sc [label:operations] report-worst-stage-status
    Use the same severity reduction for on-disk and tombstone records so a
    warning or failure in any ngmix chunk survives reclamation. Chunks share
    ``stage: "tile_ngmix"`` even though their filenames differ.

    Product counts come from the retained record only, not the sum of chunks;
    ngmix attrition therefore describes one chunk, not a whole tile. Equal
    severity keeps the first record encountered.
    """
    prev = records.get(stage)
    if prev is None or _rank(m) > _rank(prev):
        records[stage] = m


def load_manifests(run_dir: Path, sub: str) -> dict:
    """``{unit: {stage: verdict}}`` for one store (``tiles`` or ``exp``).

    Reads both of a unit's record dirs: ``manifests/`` (the rules' declared
    outputs, success-only) and ``logs/`` (the rules' ``log:``, written every run
    and never deleted by snakemake, so this is where a failure survives).

    The unit key is the record dir's *parent directory name* — shard-depth
    agnostic, and the only form that joins to the index (the record's own
    ``unit`` field carries ``SP_UNIT_NUM``'s dashed form, ``210-282``, which is
    not the index's ``210.282``). The stage comes from the body, never the
    filename: ngmix chunks share a stage under per-chunk filenames, and a log
    names the same stage as the manifest beside it.

    Several files therefore map to one (unit, stage), and the worst status wins.
    That is what collapses the ngmix chunks to one entry, and it is why a
    successful stage's two byte-identical records cost nothing while a failure
    always speaks. A body with no ``stage`` field is skipped: it is not one of
    ours, which is what keeps a stray JSON deeper in the tree inert.
    """
    out: dict = defaultdict(dict)
    paths = sorted((run_dir / sub).glob("*/*/manifests/*.json")) \
        + sorted((run_dir / sub).glob("*/*/logs/*.json"))
    for path in paths:
        try:
            m = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[run_report] unreadable record {path}: {exc}", file=sys.stderr)
            continue
        if not isinstance(m, dict) or "stage" not in m:
            continue
        unit = path.parent.parent.name
        keep_worst(out[unit], m["stage"], m)
    return out


def absorb_tombstones(run_dir: Path, sub: str, manifests: dict,
                      survivors: frozenset = frozenset()) -> set:
    """Fill in reclaimed units from their ``cleaned.json``; return their ids.

    @sc [label:provenance] report-tombstone-generation-boundary
    Absorb a tombstone only if no on-disk record belongs to a stage outside
    ``survivors``. Such a record marks a rebuilt chain; mixing in tombstone
    stages could report a previous generation's success during a rerun or
    failure. Surviving manifests alone do not mark a rebuilt tile.

    This guard is all-or-nothing. A partially rebuilt tile can report prepare
    stages as "not run" even when they executed before reclamation. Per-stage
    recovery would need generation identifiers, which these records lack.
    """
    cleaned = set()
    for path in sorted((run_dir / sub).glob("*/*/cleaned.json")):
        unit = path.parent.name
        try:
            tomb = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[run_report] unreadable tombstone {path}: {exc}", file=sys.stderr)
            continue
        if any(s not in survivors for s in manifests.get(unit, {})):
            continue
        # Reduce chunk records by severity before merging with live records.
        absorbed: dict = {}
        for key, m in (tomb.get("manifests") or {}).items():
            if not isinstance(m, dict):
                continue
            keep_worst(absorbed, m.get("stage", key), m)
        for stage, m in absorbed.items():
            # Surviving on-disk records remain authoritative.
            manifests[unit].setdefault(stage, m)
        cleaned.add(unit)
    return cleaned


def shortfalls(m: dict) -> dict:
    """``{runner: (found, expect)}`` for every runner under expect."""
    return {r: (d["found"], d["expect"])
            for r, d in m.get("runners", {}).items() if d["found"] < d["expect"]}


def reasons(m: dict) -> list:
    """Flattened failure reasons, runner-tagged, for the report's why column."""
    out = []
    for f in m.get("failures", []):
        head = f"{f['runner']} {f['found']}/{f.get('expect', '?')}"
        out += [f"{head}: {r}" for r in f["reasons"]] or [head]
    return out


def tally_level(units, stages, manifests, cleaned=frozenset()) -> dict:
    """Per-stage counts + named unit lists, for one level.

    ``cleaned`` units are counted by the status their absorbed manifests carry
    and additionally listed under ``cleaned``, so a reclaimed campaign reads as
    reclaimed rather than as a campaign that never ran. The STATUS is preserved
    exactly, including a warn on any one of the eight ngmix chunks (keep_worst
    is what makes that true on the tombstone path as well as on disk).

    See ``keep_worst`` for the single-record product-count limitation.
    """
    per_stage = {}
    for stage in stages:
        # All five status keys hold unit-ID lists; printed counts use len().
        t = {"complete": [], "warn": [], "failed": [], "not_run": [], "cleaned": []}
        agg = defaultdict(lambda: {"found": 0, "expect": 0, "by_unit": {}})
        for u in units:
            m = manifests.get(u, {}).get(stage)
            if m is None:
                t["not_run"].append(u)
                continue
            if u in cleaned:
                t["cleaned"].append(u)
            status = m.get("status", "failed")
            status = status if status in ("complete", "warn") else "failed"
            t[status].append(u)
            if status == "failed":
                # Failed units are named above, never folded into the attrition
                # aggregate: a whole-unit failure is not per-CCD attrition, and
                # mixing them hides real deletion bugs behind a big denominator.
                continue
            for runner, d in m.get("runners", {}).items():
                a = agg[runner]
                a["found"] += d["found"]
                a["expect"] += d["expect"]
                if d["found"] < d["expect"]:
                    a["by_unit"][u] = d["expect"] - d["found"]
        for a in agg.values():
            if not a["by_unit"]:
                del a["by_unit"]
        t["products"] = dict(agg)
        per_stage[stage] = t
    return per_stage


def unit_rows(units, stages, manifests) -> list:
    """One row per non-clean unit: its first bad stage, shortfalls, why."""
    rows = []
    for u in units:
        got = manifests.get(u, {})
        bad = [s for s in stages
               if got.get(s) is None or got[s].get("status") != "complete"]
        if not bad:
            continue
        stage = bad[0]
        m = got.get(stage)
        rows.append({
            "unit": u,
            "stage": stage,
            "status": "not_run" if m is None else m.get("status", "failed"),
            "shortfalls": shortfalls(m) if m else {},
            "reasons": reasons(m) if m else [],
            "n_bad_stages": len(bad),
        })
    return rows


def print_table(title, rows, limit=25):
    print(f"\n{title}  ({len(rows)} affected)")
    if not rows:
        print("  — none")
        return
    print(f"  {'unit':<14} {'stage':<20} {'status':<8} why")
    for r in rows[:limit]:
        short = ", ".join(f"{k} {v[0]}/{v[1]}" for k, v in r["shortfalls"].items())
        why = (r["reasons"][0] if r["reasons"] else short) or "-"
        print(f"  {r['unit']:<14} {r['stage']:<20} {r['status']:<8} {why[:90]}")
    if len(rows) > limit:
        print(f"  … and {len(rows) - limit} more (see the JSON report)")


def print_stage_table(title, per_stage, n_units):
    print(f"\n{title}  ({n_units} units declared)")
    print(f"  {'stage':<20} {'ok':>6} {'warn':>6} {'fail':>6} {'not run':>8} "
          f"{'cleaned':>8}  attrition")
    for stage, t in per_stage.items():
        att = [f"{r} {a['found']}/{a['expect']}"
               for r, a in t["products"].items() if a["found"] < a["expect"]]
        print(f"  {stage:<20} {len(t['complete']):>6} {len(t['warn']):>6} "
              f"{len(t['failed']):>6} {len(t['not_run']):>8} "
              f"{len(t.get('cleaned', [])):>8}  {', '.join(att)[:60]}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", required=True, type=Path)
    p.add_argument("--index", required=True, type=Path)
    p.add_argument("--status", default="manual")
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--limit", type=int, default=25,
                   help="rows per stdout table (the JSON report is complete)")
    args = p.parse_args()

    tiles, exps, tile_exp = [], [], defaultdict(list)
    if args.index.exists():
        con = sqlite3.connect(args.index, timeout=60)
        tiles = [r[0] for r in con.execute("SELECT tile_id FROM tiles ORDER BY 1")]
        exps = [r[0] for r in con.execute("SELECT exp_id FROM exposures ORDER BY 1")]
        for tile_id, exp_id in con.execute("SELECT tile_id, exp_id FROM tile_exposures"):
            tile_exp[tile_id].append(exp_id)
        con.close()
    else:
        print(f"[run_report] no index at {args.index} — reporting manifests only",
              file=sys.stderr)

    tile_m = load_manifests(args.run_dir, "tiles")
    exp_m = load_manifests(args.run_dir, "exp")
    # See absorb_tombstones for the survivor-aware generation boundary.
    cleaned_exp = absorb_tombstones(args.run_dir, "exp", exp_m)
    cleaned_tiles = absorb_tombstones(args.run_dir, "tiles", tile_m,
                                      SURVIVING_TILE_STAGES)
    tiles = tiles or sorted(tile_m)
    exps = exps or sorted(exp_m)

    missing_json = args.index.parent / "missing.json"
    missing = json.loads(missing_json.read_text()) if missing_json.exists() else []

    report = {
        "status": args.status,
        "n_tiles": len(tiles), "n_exposures": len(exps),
        "missing_tiles": missing,
        "tile_stages": tally_level(tiles, TILE_STAGES, tile_m, cleaned_tiles),
        "exp_stages": tally_level(exps, EXP_STAGES, exp_m, cleaned_exp),
        "cleaned_exposures": sorted(cleaned_exp),
        "cleaned_tiles": sorted(cleaned_tiles),
        "tiles": unit_rows(tiles, TILE_STAGES, tile_m),
        "exposures": unit_rows(exps, EXP_STAGES, exp_m),
    }

    # An exposure blocks its consuming tiles if any stage is missing or failed,
    # not merely warned: per-CCD attrition is expected at production scale.
    # Check every stage so an early warning cannot hide a later failure.
    # Cleaned exposures never block; see clean_exposure.py for deletion gates.
    def _blocks(unit):
        if unit in cleaned_exp:
            return False
        for stage in EXP_STAGES:
            m = exp_m.get(unit, {}).get(stage)
            if m is None or m.get("status", "failed") == "failed":
                return True
        return False

    bad_exp = {e for e in exps if _blocks(e)}
    blocked = {t: sorted(set(tile_exp.get(t, [])) & bad_exp) for t in tiles}
    report["tiles_blocked_by_exposures"] = {t: e for t, e in blocked.items() if e}

    done = len(report["tile_stages"]["tile_make_cat"]["complete"])
    report["final_cats"] = {"present": done, "of": len(tiles)}

    out = args.out or (args.index.parent / "run_report.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"[run_report] status={args.status}  {done}/{len(tiles)} final cats"
          + (f"  ({len(missing)} tiles missing exposure lists)" if missing else "")
          + (f"  ({len(cleaned_exp)} exposures reclaimed)" if cleaned_exp else "")
          + (f"  ({len(cleaned_tiles)} tiles reclaimed)" if cleaned_tiles else ""))
    print_stage_table("EXPOSURES", report["exp_stages"], len(exps))
    print_stage_table("TILES", report["tile_stages"], len(tiles))
    print_table("exposures not complete", report["exposures"], args.limit)
    print_table("tiles not complete", report["tiles"], args.limit)
    nb = report["tiles_blocked_by_exposures"]
    if nb:
        print(f"\ntiles waiting on incomplete exposures  ({len(nb)})")
        for t, e in list(nb.items())[:args.limit]:
            print(f"  {t:<14} {', '.join(e[:6])}"
                  + (f" (+{len(e) - 6})" if len(e) > 6 else ""))
    print(f"\n[run_report] -> {out}")


if __name__ == "__main__":
    main()
