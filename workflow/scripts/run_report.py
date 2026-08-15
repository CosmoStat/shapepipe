#!/usr/bin/env python3
"""``sp report`` — the run's success/failure tables, read from the manifests.

NOT a DAG node. A report rule that declared all tiles' outputs as inputs would be
a descendant of every job, so one hard failure under --keep-going would poison
its cone and the report would never run — the exact scenario it exists for. So it
is a plain script, runnable at any time including mid-run; the Snakefile's
onsuccess/onerror hooks call it so every invocation ends with one.

It reads two things and nothing else (PRD D3):

  * the **index** (``run_index.sqlite``) — the units the run declared, and the
    tile->exposure edges that let an exposure failure be blamed on the tiles it
    blocks;
  * the **manifests** — one per unit per stage, written by
    ``completeness.py check``, carrying per-runner found/expect/floor and
    log-scraped failure reasons.

A stage that failed leaves ``<stage>.failed.json`` instead of ``<stage>.json``
(``completeness.py``): the success name is reserved for success, so the DAG
cannot build on a failure. Both names are read here, and the status comes from
the manifest BODY, so a stage is reported as failed either way — including
manifests written before that split, where the failure content sat under the
success name.

A reclaimed exposure has no manifests: ``clean_exposure`` deleted them after
copying them into ``<exp dir>/cleaned.json``. That tombstone is read as the
unit's record and the unit is reported as **cleaned** — not "not run", and it
blocks no tile.

No disk scanning: counting products is the *check's* job, done once at the moment
the products were fresh. A unit with no manifest for a stage is "not run" — which
is a real and distinct answer from "ran and produced nothing".

Manifests are discovered by glob (``tiles/**/manifests/*.json``), not by
constructed path: the unit stores are sharded (``tiles/<prefix>/<ID>/``) and the
sharding depth is not this script's business.
"""

import argparse
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

# Stage order per level — the report's column order, and the definition of
# "expected" (a declared unit with no manifest for one of these is not run).
TILE_STAGES = ["tile_get_images", "tile_uncompress", "tile_find_exposures",
               "tile_merge_headers", "tile_detect", "tile_vignets",
               "tile_ngmix", "tile_merge_cats", "tile_make_cat"]
EXP_STAGES = ["exp_get_images", "exp_star_cat", "exp_split", "exp_mask", "exp_psf"]

STATUSES = ("complete", "warn", "failed", "not_run")


def load_manifests(run_dir: Path, sub: str) -> dict:
    """``{unit: {stage: manifest}}`` for one store (``tiles`` or ``exp``).

    The unit key is the manifest dir's *parent directory name* — shard-depth
    agnostic, and the only form that joins to the index (the manifest's own
    ``unit`` field carries ``SP_UNIT_NUM``'s dashed form, ``210-282``, which is
    not the index's ``210.282``). The stage comes from the manifest body, never
    the filename: ngmix chunks share a stage under per-chunk filenames, and
    ``<stage>.failed.json`` names the same stage as ``<stage>.json``.

    Several files can therefore map to one (unit, stage), and the WORST status
    wins. That is what collapses the ngmix chunks to one entry, and what makes a
    failed manifest speak even in the transient window where its success-named
    counterpart has not yet been removed.
    """
    out: dict = defaultdict(dict)
    for path in sorted((run_dir / sub).glob("**/manifests/*.json")):
        try:
            m = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[run_report] unreadable manifest {path}: {exc}", file=sys.stderr)
            continue
        unit = path.parent.parent.name
        stage = m.get("stage", path.stem)
        prev = out[unit].get(stage)
        # Worst status wins when several manifests share a stage (ngmix chunks).
        rank = lambda d: STATUSES.index(d.get("status")) if d.get("status") in STATUSES else len(STATUSES)  # noqa: E731
        if prev is None or rank(m) > rank(prev):
            out[unit][stage] = m
    return out


def absorb_tombstones(run_dir: Path, sub: str, manifests: dict) -> set:
    """Fill in reclaimed units from their ``cleaned.json``; return their ids.

    A cleaned exposure has NO ``manifests/`` — ``clean_exposure`` deleted it,
    after copying every manifest verbatim into the tombstone. Read them back, or
    the report inverts the truth exactly when reclamation works: the exposure
    shows as "not run" and blocks the very tiles whose completion authorised the
    deletion.

    Manifests on disk win if both exist — that is a re-built chain, and the
    tombstone is then a stale record of the previous generation.
    """
    cleaned = set()
    for path in sorted((run_dir / sub).glob("**/cleaned.json")):
        unit = path.parent.name
        try:
            tomb = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            print(f"[run_report] unreadable tombstone {path}: {exc}", file=sys.stderr)
            continue
        if manifests.get(unit):
            continue
        for key, m in (tomb.get("manifests") or {}).items():
            if not isinstance(m, dict):
                continue
            manifests[unit][m.get("stage", key)] = m
        cleaned.add(unit)
    return cleaned


def shortfalls(m: dict) -> dict:
    """``{runner: (found, expect, floor)}`` for every runner under expect."""
    return {r: (d["found"], d["expect"], d["floor"])
            for r, d in m.get("runners", {}).items() if d["found"] < d["expect"]}


def reasons(m: dict) -> list:
    """Flattened failure reasons, runner-tagged, for the report's why column."""
    out = []
    for f in m.get("failures", []):
        head = f"{f['runner']} {f['found']}/{f.get('expect', '?')} (floor {f['floor']})"
        out += [f"{head}: {r}" for r in f["reasons"]] or [head]
    return out


def tally_level(units, stages, manifests, cleaned=frozenset()) -> dict:
    """Per-stage counts + named unit lists, for one level.

    ``cleaned`` units are counted by the status their absorbed manifests carry
    (complete or warn — attrition is preserved) and additionally listed under
    ``cleaned``, so a reclaimed campaign reads as reclaimed rather than as a
    campaign that never ran.
    """
    per_stage = {}
    for stage in stages:
        t = {"complete": 0, "warn": [], "failed": [], "not_run": [], "cleaned": []}
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
            if status == "complete":
                t["complete"] += 1
            else:
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
        print(f"  {stage:<20} {t['complete']:>6} {len(t['warn']):>6} "
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
    # Reclaimed exposures speak through their tombstones (D5, S5).
    cleaned_exp = absorb_tombstones(args.run_dir, "exp", exp_m)
    tiles = tiles or sorted(tile_m)
    exps = exps or sorted(exp_m)

    missing_json = args.index.parent / "missing.json"
    missing = json.loads(missing_json.read_text()) if missing_json.exists() else []

    report = {
        "status": args.status,
        "n_tiles": len(tiles), "n_exposures": len(exps),
        "missing_tiles": missing,
        "tile_stages": tally_level(tiles, TILE_STAGES, tile_m),
        "exp_stages": tally_level(exps, EXP_STAGES, exp_m, cleaned_exp),
        "cleaned_exposures": sorted(cleaned_exp),
        "tiles": unit_rows(tiles, TILE_STAGES, tile_m),
        "exposures": unit_rows(exps, EXP_STAGES, exp_m),
    }

    # Blame propagation: a BLOCKING exposure blocks every tile that reads it.
    # Without this, a tile stalled at tile_vignets looks like its own failure.
    #
    # Blocking means "failed" or "never ran" — NOT "warn". Warn is the expected
    # per-CCD attrition (setools rejecting a sparse CCD, psfex_interp short an
    # epoch); it is present in essentially every exposure at production scale, so
    # counting it here made every exposure block every tile and the table said
    # nothing.
    # Judged over ALL the exposure's stages, not just the first bad one, so an
    # exposure that warns early and fails late still blocks.
    # A CLEANED exposure never blocks: its store is gone precisely because every
    # consuming tile already had its vignets. Its absorbed manifests are read
    # above, so a cleaned exposure that genuinely failed still shows in the
    # tables — it just does not get to hold complete tiles hostage.
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

    done = report["tile_stages"]["tile_make_cat"]["complete"]
    report["final_cats"] = {"present": done, "of": len(tiles)}

    out = args.out or (args.index.parent / "run_report.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print(f"[run_report] status={args.status}  {done}/{len(tiles)} final cats"
          + (f"  ({len(missing)} tiles missing exposure lists)" if missing else "")
          + (f"  ({len(cleaned_exp)} exposures reclaimed)" if cleaned_exp else ""))
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
