#!/usr/bin/env python
"""Extract the status-board metric set and UPSERT it into a longitudinal ``history.json``.

The status page (``build_status.py``) renders one snapshot per run; this script
is its complement across *time*. It reads the same run's result JSONs, distils
them to the small, well-defined **metric set v1** -- one scalar per metric, each
carrying a value, an optional 1-sigma error, and an optional tolerance -- and
folds those rows into a growing ``history.json`` keyed by ``(commit, metric)``.

The UPSERT is the whole contract: rebuilding the *same* commit updates its rows
in place; a *new* commit appends; nothing is ever truncated. So the file
accumulates a trend as commits land, and re-running a build is idempotent.

Three row sources feed the same UPSERT:

* **current build** -- the metric-set-v1 rows harvested from ``--results`` at
  ``--commit`` / ``--date`` (the live path);
* **``--backfill-status-json``** -- one historical ``status.json`` blob whose
  scalars predate any error bars, folded in as ``err=null`` rows (stored, but
  the plotter must not render them: no number without an uncertainty);
* **``--imsims-anchors``** -- the SKiLLS image-sim reference points, folded in as
  ``preview:true`` rows (rendered, but chip-labelled "SKiLLS . preliminary").

Rows with ``err==null`` are stored for provenance but the plotter skips them.
``preview`` rows are rendered and visually flagged. Design contract, mirrored
from ``build_status.py``: every number is read straight from a results JSON;
required keys are indexed (``d["key"]``), not ``.get(...,default)``, so a
missing field fails loud at the moment it goes missing -- *except* where a
source legitimately lacks a field (a degenerate grid cell whose ``m`` is
``null``), which is filtered out before the read.

Runs inside the ShapePipe sif (numpy + stdlib, no seaborn) via ``srun``::

    srun --jobid=<J> -n 1 -c 4 apptainer exec --bind ... <sif> \\
        python scripts/python/build_history.py \\
            --results <run>/results/board \\
            --out <proto>/history.json \\
            --commit v1.4.0-639-ga59bed9a --date 2026-07-10 \\
            --backfill-status-json <proto>/hist_0706_status.json \\
            --imsims-anchors
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# Resolved-cell gate for the breakdown grid: a cell counts toward a band's mean
# only if its intrinsic resolution clears this ratio *and* its degeneracy guard
# did not fire (m not null). Mirrors the status-page "resolved non-degenerate"
# set. res_ratio 0.3657 sits just under the gate, so the resolved set starts at
# the 0.492 grid point.
RES_RATIO_MIN = 0.37

# The four star-response arms, mapped from their result subdir to the short arm
# label the metric name carries.
STAR_ARMS = {
    "star_response_grid_wcs": "wcs",
    "star_response_grid_hsm": "hsm",
    "star_response_grid_gauss_wcs": "gauss_wcs",
    "star_response_grid_gauss_hsm": "gauss_hsm",
}

# The imsims / SKiLLS preview anchors: 3 build labels x (m1, m2), rendered but
# chip-labelled preliminary. commit == the build-label string; date per row.
IMSIMS_ANCHORS = [
    ("2026-07-10", "im_sims pre-correction", +0.27436, 0.028, +0.23122, 0.022),
    ("2026-07-10", "amplitude-corrected", +0.01949, 0.02237, -0.01502, 0.01747),
    ("2026-07-11", "gate766", +0.01966, 0.02120, -0.01582, 0.01746),
]


def load(results: Path, arm: str, fname: str) -> dict | None:
    """Load ``<results>/<arm>/<fname>``; ``None`` (not an error) if absent."""
    p = results / arm / fname
    return json.loads(p.read_text()) if p.is_file() else None


def row(date: str, commit: str, metric: str, value, err, tol, preview: bool = False) -> dict:
    """One history row. ``value``/``err``/``tol`` coerced to plain floats (or
    ``None``); ``preview`` only set when true, to keep the schema minimal."""
    r = {
        "date": date,
        "commit": commit,
        "metric": metric,
        "value": None if value is None else float(value),
        "err": None if err is None else float(err),
        "tol": None if tol is None else float(tol),
    }
    if preview:
        r["preview"] = True
    return r


def ivw(values: np.ndarray, errs: np.ndarray) -> tuple[float, float]:
    """Inverse-variance-weighted mean and its combined 1-sigma error.

    ``w = 1/err^2``; ``mean = sum(w x)/sum(w)``; ``err = sqrt(1/sum(w))``."""
    w = 1.0 / errs**2
    return float(np.sum(w * values) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


# --------------------------------------------------------------------------- #
# Current-build extraction: results dir -> metric-set-v1 rows
# --------------------------------------------------------------------------- #
def rows_mbias(d: dict, date: str, commit: str) -> list[dict]:
    """Tier-1 m-bias arm: the multiplicative bias and the metacal response. Both
    err=null (single-point in-memory sim, no bootstrap on the arm scalar)."""
    return [
        row(date, commit, "mbias_m", d["multiplicative_bias_m"], None, 0.005),
        row(date, commit, "mbias_R11", d["metacal_response_R11"], None, None),
    ]


def rows_breakdown(d: dict, date: str, commit: str) -> list[dict]:
    """Breakdown grid (v2 cells): per-band IVW m1/m2 over the resolved
    non-degenerate set, the max-|c1| order statistic, and alpha at the fiducial
    cell."""
    cells = d["cells"]
    m_tol = d["meta"]["m_tol_resolved"]
    out: list[dict] = []

    for side in ("m1", "m2"):
        for band in ("high", "moderate", "low"):
            # Resolved, non-degenerate, and carrying a real value+err. The
            # degenerate guard (m null) is a legitimate missing field -- filtered
            # here, so the value/err reads below are always safe.
            sel = [c for c in cells
                   if c["s2n_band"] == band
                   and c["res_ratio"] >= RES_RATIO_MIN
                   and not c[f"{side}_degenerate"]
                   and c[side] is not None and c[f"{side}_err"] is not None]
            if not sel:
                continue
            vals = np.array([c[side] for c in sel], dtype=float)
            errs = np.array([c[f"{side}_err"] for c in sel], dtype=float)
            mean, err = ivw(vals, errs)
            out.append(row(date, commit, f"breakdown_{side}_{band}", mean, err, m_tol))

    # max|c1| over ALL cells (order statistic -> err null, tol null).
    c1s = [c["c1"] for c in cells if c["c1"] is not None]
    out.append(row(date, commit, "breakdown_max_abs_c1", max(np.abs(c1s)), None, None))

    # Fiducial alpha: band=moderate, nearest res_ratio to 1.0 among cells that
    # are non-degenerate on m1 and carry alpha+alpha_err. Skip (do not invent) if
    # none qualify.
    fid_pool = [c for c in cells
                if c["s2n_band"] == "moderate"
                and not c["m1_degenerate"]
                and c["alpha"] is not None and c["alpha_err"] is not None]
    if fid_pool:
        fid = min(fid_pool, key=lambda c: abs(c["res_ratio"] - 1.0))
        out.append(row(date, commit, "breakdown_alpha_fid", fid["alpha"], fid["alpha_err"], None))
        print(f"  [breakdown_alpha_fid] fiducial cell = index {fid['cell_index']} "
              f"(band=moderate, res_ratio={fid['res_ratio']:.4f}, nearest 1.0)")
    else:
        print("  [breakdown_alpha_fid] no qualifying fiducial cell -- metric skipped")

    return out


def rows_star(results: Path, date: str, commit: str) -> list[dict]:
    """Star-response arms: mean R1/R2 with jackknife errors, per arm present."""
    out: list[dict] = []
    for subdir, arm in STAR_ARMS.items():
        d = load(results, subdir, "star_response.json")
        if d is None:
            continue
        tol = d["tolerance_R"]
        for comp in ("R1", "R2"):
            out.append(row(date, commit, f"star_{comp}_mean_{arm}",
                           d[comp]["mean"], d[comp]["jk_err"], tol))
    return out


def build_current(results: Path, date: str, commit: str) -> list[dict]:
    """All metric-set-v1 rows for the current build under ``results``."""
    out: list[dict] = []

    mbias = load(results, "mbias_sim", "mbias.json")
    if mbias is not None:
        out += rows_mbias(mbias, date, commit)

    breakdown = load(results, "breakdown_grid", "breakdown_grid.json")
    if breakdown is not None:
        out += rows_breakdown(breakdown, date, commit)

    out += rows_star(results, date, commit)
    return out


# --------------------------------------------------------------------------- #
# Backfill sources
# --------------------------------------------------------------------------- #
def build_backfill_status(blob: dict) -> list[dict]:
    """Historical 07-06 status blob -> err=null rows (stored, not rendered).

    Synthetic commit label (no git_describe existed at that snapshot)."""
    date, commit = "2026-07-06", "status-2026-07-06"
    lvl1 = blob["level1a_mbias"]
    grid = blob["level2_grid"]
    out = [
        row(date, commit, "mbias_m", lvl1["m"], None, 0.005),
        row(date, commit, "mbias_R11", lvl1["R11"], None, None),
    ]
    for arm in ("wcs", "hsm", "gauss_wcs", "gauss_hsm"):
        a = grid[arm]
        tol = a["tol"]
        out.append(row(date, commit, f"star_R1_mean_{arm}", a["R1_mean"], None, tol))
        out.append(row(date, commit, f"star_R2_mean_{arm}", a["R2_mean"], None, tol))
    return out


def build_imsims() -> list[dict]:
    """The 3x2 imsims preview anchors: preview=true, rendered, chip-labelled."""
    out: list[dict] = []
    for date, label, m1, m1_err, m2, m2_err in IMSIMS_ANCHORS:
        out.append(row(date, label, "imsims_ref_m1", m1, m1_err, None, preview=True))
        out.append(row(date, label, "imsims_ref_m2", m2, m2_err, None, preview=True))
    return out


# --------------------------------------------------------------------------- #
# UPSERT + IO
# --------------------------------------------------------------------------- #
def upsert(existing: list[dict], new: list[dict]) -> list[dict]:
    """Fold ``new`` rows into ``existing``, keyed by ``(commit, metric)``: a
    matching key updates in place, a fresh key appends. Never truncates."""
    index = {(r["commit"], r["metric"]): i for i, r in enumerate(existing)}
    merged = list(existing)
    for r in new:
        key = (r["commit"], r["metric"])
        if key in index:
            merged[index[key]] = r
        else:
            index[key] = len(merged)
            merged.append(r)
    return merged


def load_history(path: Path) -> list[dict]:
    """Existing rows from ``history.json`` (schema 1), or [] if absent."""
    if not path.is_file():
        return []
    doc = json.loads(path.read_text())
    assert doc["schema"] == 1, f"unexpected history schema {doc['schema']}"
    return doc["rows"]


def print_summary(rows: list[dict]) -> None:
    """Print the metric table for the current invocation's rows."""
    print(f"\n  {'metric':<28} {'value':>12} {'err':>12} {'tol':>8}  preview")
    print("  " + "-" * 76)
    for r in sorted(rows, key=lambda r: (r["date"], r["metric"])):
        v = "None" if r["value"] is None else f"{r['value']:+.5f}"
        e = "null" if r["err"] is None else f"{r['err']:.5f}"
        t = "null" if r["tol"] is None else f"{r['tol']:.4f}"
        p = "yes" if r.get("preview") else ""
        print(f"  {r['metric']:<28} {v:>12} {e:>12} {t:>8}  {p}")
    n_err = sum(r["err"] is not None for r in rows)
    print(f"\n  {n_err}/{len(rows)} rows carry an error; {len(rows) - n_err} have err=null.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", type=Path, required=True,
                    help="board results dir (per-arm subdirs of result JSONs)")
    ap.add_argument("--out", type=Path, required=True, help="history.json path")
    ap.add_argument("--commit", default=None, help="git_describe label for the current build")
    ap.add_argument("--date", default=None, help="YYYY-MM-DD of the current build")
    ap.add_argument("--backfill-status-json", type=Path, default=None,
                    help="historical status.json blob to fold in as err=null rows")
    ap.add_argument("--imsims-anchors", action="store_true",
                    help="fold in the 3x2 imsims preview anchors")
    a = ap.parse_args(argv)

    new: list[dict] = []

    if a.commit or a.date:
        assert a.commit and a.date, "--commit and --date must be given together"
        assert a.results.is_dir(), f"--results dir not found: {a.results}"
        current = build_current(a.results, a.date, a.commit)
        print(f"\ncurrent build [{a.date} . {a.commit}]: {len(current)} rows")
        new += current

    if a.backfill_status_json is not None:
        blob = json.loads(a.backfill_status_json.read_text())
        backfill = build_backfill_status(blob)
        print(f"backfill status blob [2026-07-06 . status-2026-07-06]: {len(backfill)} rows (err=null)")
        new += backfill

    if a.imsims_anchors:
        imsims = build_imsims()
        print(f"imsims preview anchors: {len(imsims)} rows (preview=true)")
        new += imsims

    assert new, "nothing to do: give --commit/--date, --backfill-status-json, or --imsims-anchors"

    existing = load_history(a.out)
    merged = upsert(existing, new)
    merged.sort(key=lambda r: (r["date"], r["metric"], r["commit"]))

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({"schema": 1, "rows": merged}, indent=1))

    print_summary(new)
    print(f"\nwrote {a.out}: {len(existing)} existing + {len(new)} upserted -> {len(merged)} total rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
