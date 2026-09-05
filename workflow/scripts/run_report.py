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
  * the **verdicts** written by ``completeness.py check`` — per-runner
    found/expect/floor and log-scraped failure reasons — in the unit's
    ``logs/`` (every run) and ``manifests/`` (successes only; completeness.py
    argues the split).

Both dirs are read here and the status comes from the file BODY, so a failed
stage speaks through its log while a successful one is corroborated by two
identical files. A unit with neither ran nothing. Also read, for continuity with
stores written before this convention: ``manifests/<stage>.failed.json``, which
is where failure evidence used to go.

A reclaimed exposure has neither dir: ``clean_exposure`` deleted both, after
copying the manifests into ``<exp dir>/cleaned.json``. That tombstone is read as the
unit's record and the unit is reported as **cleaned** — not "not run", and it
blocks no tile.

A reclaimed TILE is the same story with one asymmetry: ``clean_tile`` deletes
``logs/`` but cannot empty ``manifests/``, because two of those manifests are
currency other mechanisms read (SURVIVING_TILE_STAGES). So the "records on disk
mean a rebuilt chain" test takes that survivor set explicitly — see
``absorb_tombstones``.

No disk scanning: counting products is the *check's* job, done once at the moment
the products were fresh. A unit with no manifest for a stage is "not run" — which
is a real and distinct answer from "ran and produced nothing".

Records are discovered by glob (``tiles/*/*/manifests/*.json`` and
``tiles/*/*/logs/*.json``), not by constructed path: units are found rather than
named, so a store holding a unit the index never heard of still reports. The
depth is FIXED at two, because the sharded layout is exactly two levels
(``tiles/<prefix>/<ID>/``) — a ``**`` walked the whole tree, including every
``output/`` a run has not reclaimed yet, to find records that can only ever be
at one depth.
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

# exp_persist and exp_footprint are DELIBERATELY NOT in that list. This report
# disk-scans the scratch run_dir, and those two are the exposure manifests that
# live on products_dir instead — that placement is what makes them survive
# clean_exposure. Listed here they would read as "not run" for every exposure
# in the campaign. Reporting on the persisted products means scanning the
# second root, which is a report this one does not yet do.

# The manifests clean_tile leaves on disk (workflow/scripts/clean_tile.py names
# the mechanism that owns each). Their presence is therefore NOT evidence that a
# tile's chain was rebuilt, which is the one thing absorb_tombstones has to know
# to read a reclaimed tile's record out of its tombstone.
SURVIVING_TILE_STAGES = frozenset({"tile_vignets", "tile_find_exposures"})

STATUSES = ("complete", "warn", "failed", "not_run")


def _rank(m: dict) -> int:
    """How BAD a record is, as a position in STATUSES; unknown sorts worst."""
    return STATUSES.index(m["status"]) if m.get("status") in STATUSES else len(STATUSES)


def keep_worst(records: dict, stage: str, m: dict) -> None:
    """Collapse the several records of one stage into the WORST of them.

    ONE rule, used by both readers in this module, and that is the point.
    Several files map to one (unit, stage) either way: on disk, a stage's
    manifest and its byte-identical log, plus the eight ngmix chunks, which all
    carry ``stage: "tile_ngmix"`` under per-chunk filenames; inside a tombstone,
    those same eight chunks again, keyed by file stem.

    absorb_tombstones used to resolve that collision by first-key-wins
    (``setdefault`` over sorted stems), so only ``tile_ngmix_1`` survived and a
    ``warn`` on any other chunk disappeared the moment the tile was reclaimed —
    ``tile_ngmix ok 0 warn 1`` before the clean, ``ok 1 warn 0`` after, and the
    tile silently left the "tiles not complete" table. The data was in the
    tombstone the whole time; only the reader dropped it (found in review,
    reproduced on 186.307 with chunk 5 flipped to warn).

    What this does NOT fix, because it is not a reclamation bug: a stage's
    per-runner ``products`` aggregate is taken from the single surviving record,
    so tile_ngmix attrition is counted over one chunk of eight. That is true
    before and after a clean, identically — it is the price of collapsing the
    chunks to one stage row, and it is the same on both paths by construction
    now.
    """
    prev = records.get(stage)
    if prev is None or _rank(m) > _rank(prev):
        records[stage] = m


def load_manifests(run_dir: Path, sub: str) -> dict:
    """``{unit: {stage: verdict}}`` for one store (``tiles`` or ``exp``).

    Reads BOTH of a unit's record dirs: ``manifests/`` (the rules' declared
    outputs, success-only) and ``logs/`` (the rules' ``log:``, written every run
    and never deleted by snakemake, so this is where a failure survives).

    The unit key is the record dir's *parent directory name* — shard-depth
    agnostic, and the only form that joins to the index (the record's own
    ``unit`` field carries ``SP_UNIT_NUM``'s dashed form, ``210-282``, which is
    not the index's ``210.282``). The stage comes from the body, never the
    filename: ngmix chunks share a stage under per-chunk filenames, and a log
    names the same stage as the manifest beside it.

    Several files therefore map to one (unit, stage), and the WORST status wins.
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

    A cleaned exposure has neither ``manifests/`` nor ``logs/``; every manifest
    was copied verbatim into the tombstone first. Read them back, or the report
    inverts the truth exactly when reclamation works: the exposure shows as "not
    run" and blocks the very tiles whose completion authorised the deletion.

    Manifests on disk win if both exist — that is a re-built chain, and the
    tombstone is then a stale record of the previous generation.

    ``survivors`` is what makes that test work for TILES. ``clean_tile`` cannot
    empty ``manifests/`` the way ``clean_exposure`` does: two of the manifests
    there are currency other mechanisms read (SURVIVING_TILE_STAGES below), so a
    reclaimed tile ALWAYS has records on disk. Without this argument the
    "manifests win" guard fires on every cleaned tile, the tombstone is ignored,
    and a fully reclaimed campaign reports as one that ran two stages and
    stopped. A unit counts as REBUILT — and keeps the guard — iff it has a record
    for some stage that is not a survivor; the default empty set reproduces the
    exposure test exactly.

    KNOWN AND DELIBERATE: the guard is all-or-nothing, so a PARTIALLY rebuilt
    cleaned tile loses its whole tombstone record. Force a rerun that restores
    tile_detect..tile_make_cat but not the prepare stages (which only the
    prepare invocation produces) and tile_get_images / tile_uncompress read
    "not run" though they did run and nothing invalidated them.

    Not fixed, because the obvious fix is wrong rather than long. Filling
    per-stage from the tombstone would, on a unit whose chain is rebuilding,
    report the PREVIOUS generation's "complete" for a stage whose manifest is
    absent precisely because it is mid-rerun or failed — a stale complete is
    invisibly wrong where "not run" is visibly incomplete, and the same hazard
    reaches exposures, where the mixing would be silent and campaign-wide. The
    all-or-nothing guard is a generation boundary and is worth more than the
    cosmetics. A correct fix needs a per-stage notion of which generation a
    record belongs to, which nothing here records today.
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
        # Collapse the tombstone's per-stem records to one per stage FIRST,
        # by the same rule the on-disk path uses (keep_worst), and only then
        # merge. Reversing those two steps is what lost a warning chunk.
        absorbed: dict = {}
        for key, m in (tomb.get("manifests") or {}).items():
            if not isinstance(m, dict):
                continue
            keep_worst(absorbed, m.get("stage", key), m)
        for stage, m in absorbed.items():
            # setdefault, not assignment: the surviving on-disk records are the
            # live ones and stay authoritative, even though the tombstone's
            # copies of them are byte-identical today.
            manifests[unit].setdefault(stage, m)
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
    and additionally listed under ``cleaned``, so a reclaimed campaign reads as
    reclaimed rather than as a campaign that never ran. The STATUS is preserved
    exactly, including a warn on any one of the eight ngmix chunks (keep_worst
    is what makes that true on the tombstone path as well as on disk).

    The per-runner ``products`` aggregate is NOT a per-chunk total: reading it
    as a whole-tile figure over-states completeness by 8x (see ``keep_worst``).
    """
    per_stage = {}
    for stage in stages:
        # All five status keys are unit-id LISTS, "complete" included: it used
        # to be a bare int, which made it the one key a caller had to special-
        # case. The emitted JSON gains the complete-unit list; the printed
        # counts are len() of it.
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
    # Reclaimed units speak through their tombstones (D5, S5). Tiles pass the
    # survivor set: clean_tile leaves two manifests on disk on purpose, and
    # without that argument they would read as a rebuilt chain.
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
