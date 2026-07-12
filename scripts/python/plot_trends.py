#!/usr/bin/env python
"""Small-multiple trends figure from the board ``history.json``.

The companion to ``build_history.py``: that script accumulates the status-board
metric set across commits; this one draws it as a grid of small multiples, one
panel per metric *family*, x = build (date + short commit), y = value with its
1-sigma error bar. Where a family carries a tolerance the panel shades the
``+/-tol`` band as a sliver at 0.

**Page rule, enforced here: no number without an uncertainty.** Rows whose
``err`` is null are stored in the history for provenance but never plotted -- so
the m-bias family (both single-point, err null) renders empty and says so, and
the historical 07-06 star scalars drop out while their 07-10 jackknifed
counterparts remain. ``preview`` rows (the imsims / SKiLLS anchors) *are* drawn,
but their panel is chip-labelled ``SKiLLS . preliminary`` so a reader never
mistakes a preliminary reference point for a vetted board metric.

The caption states what is shown; it offers no interpretation.

Matplotlib only (Agg backend, board palette by hand) so it runs inside the
ShapePipe sif -- seaborn is not installed there::

    srun --jobid=<J> -n 1 -c 4 apptainer exec --bind ... <sif> \\
        python scripts/python/plot_trends.py \\
            --history <proto>/history.json \\
            --output <proto>/trends.png
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Board palette (mirrored from the status-page CSS). One hue per series within a
# panel; families reuse the palette across panels.
CINNABAR = "#BC4538"
TEAL = "#3F8278"
COBALT = "#3D5BA0"
OCHRE = "#C49333"
INK = "#1B1611"
MUTED = "#7A6F5C"
RULE = "#C5B89E"

# A metric family = one panel. Each maps a title + y-label to the ordered series
# it draws: (series-label, metric-name, colour). Series that carry no err-bearing
# row simply do not appear (page rule), which is how the m-bias panel renders
# empty and the star panels collapse to their jackknifed build.
FAMILIES = [
    ("m-bias (in-memory tripwire)", "value",
     [("m", "mbias_m", CINNABAR)]),
    ("Breakdown $m_1$ by S/N band", "IVW mean $m_1$",
     [("high", "breakdown_m1_high", TEAL),
      ("moderate", "breakdown_m1_moderate", COBALT),
      ("low", "breakdown_m1_low", CINNABAR)]),
    ("Breakdown $m_2$ by S/N band", "IVW mean $m_2$",
     [("high", "breakdown_m2_high", TEAL),
      ("moderate", "breakdown_m2_moderate", COBALT),
      ("low", "breakdown_m2_low", CINNABAR)]),
    (r"Star $\langle R_1\rangle$ by arm", r"$\langle R_1\rangle$",
     [("wcs", "star_R1_mean_wcs", CINNABAR),
      ("hsm", "star_R1_mean_hsm", OCHRE),
      ("gauss_wcs", "star_R1_mean_gauss_wcs", TEAL),
      ("gauss_hsm", "star_R1_mean_gauss_hsm", COBALT)]),
    (r"Star $\langle R_2\rangle$ by arm", r"$\langle R_2\rangle$",
     [("wcs", "star_R2_mean_wcs", CINNABAR),
      ("hsm", "star_R2_mean_hsm", OCHRE),
      ("gauss_wcs", "star_R2_mean_gauss_wcs", TEAL),
      ("gauss_hsm", "star_R2_mean_gauss_hsm", COBALT)]),
    ("imsims reference $m_1$, $m_2$", "reference $m$",
     [("$m_1$", "imsims_ref_m1", CINNABAR),
      ("$m_2$", "imsims_ref_m2", COBALT)]),
]


# Build-narrative order within a shared date, where alphabetical commit-sort
# would otherwise scramble the physical sequence. The imsims anchors are the one
# multi-build same-date case: pre-correction is the *starting* bias, then the
# amplitude correction, so it must sit left of amplitude-corrected. Any commit
# not listed keeps its natural (alphabetical / git-describe) order after these.
BUILD_ORDER = {
    "im_sims pre-correction": 0,
    "amplitude-corrected": 1,
    "gate766": 2,
}


def build_x_axis(rows: list[dict]) -> tuple[dict, list[str]]:
    """Shared x-axis for the whole figure: every distinct ``(date, commit)`` in
    the history, ordered by date then by build-narrative order (``BUILD_ORDER``,
    falling back to commit name), mapped to an integer position and an
    ``MM-DD\\n<short-commit>`` tick label."""
    seen: list[tuple[str, str]] = []
    for r in rows:
        key = (r["date"], r["commit"])
        if key not in seen:
            seen.append(key)
    order = sorted(seen, key=lambda k: (k[0], BUILD_ORDER.get(k[1], len(BUILD_ORDER)), k[1]))
    pos = {k: i for i, k in enumerate(order)}
    labels = [f"{d[5:]}\n{c[:12]}" for d, c in order]   # MM-DD + short commit
    return pos, labels


def render_only(rows: list[dict]) -> list[dict]:
    """The page rule: keep only rows that carry an uncertainty."""
    return [r for r in rows if r["err"] is not None]


def draw_panel(ax, title: str, ylabel: str, series, by_metric, xpos,
               n_x: int, is_preview: bool) -> int:
    """Draw one metric family. Returns the number of points plotted so the
    caller can flag an empty (err-null-only) family."""
    n_pts = 0
    tol = None
    for label, metric, colour in series:
        rows = sorted(by_metric.get(metric, []), key=lambda r: xpos[(r["date"], r["commit"])])
        if not rows:
            continue
        x = [xpos[(r["date"], r["commit"])] for r in rows]
        y = [r["value"] for r in rows]
        e = [r["err"] for r in rows]
        ax.errorbar(x, y, yerr=e, color=colour, marker="o", ms=6, lw=1.6,
                    capsize=3, elinewidth=1.3, label=label, zorder=3)
        n_pts += len(rows)
        tol = tol or next((r["tol"] for r in rows if r["tol"] is not None), None)

    ax.set_title(title, fontsize=12, color=INK, loc="left")
    ax.set_ylabel(ylabel, fontsize=10, color=INK)
    ax.axhline(0, color=MUTED, lw=0.8, zorder=0)
    if tol is not None:
        ax.axhspan(-tol, tol, color=RULE, alpha=0.35, zorder=0,
                   label=fr"$\pm${tol:g} tol.")
    ax.set_xlim(-0.5, n_x - 0.5)
    ax.margins(y=0.25)
    ax.tick_params(labelsize=8, colors=INK)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    if n_pts == 0:
        ax.text(0.5, 0.5, "err = null\n(not rendered)", transform=ax.transAxes,
                ha="center", va="center", color=MUTED, fontsize=10, style="italic")
    else:
        ax.legend(frameon=False, fontsize=7.5, loc="best", ncol=2)

    if is_preview and n_pts:
        # Chip: preview panels carry a preliminary flag so a reader never reads
        # them as a vetted board metric.
        ax.text(0.015, 0.965, "SKiLLS . preliminary", transform=ax.transAxes,
                ha="left", va="top", fontsize=8, color="#8a6a1f", family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", fc=(0.77, 0.58, 0.20, 0.16),
                          ec=OCHRE, lw=0.8))
    return n_pts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--history", type=Path, required=True, help="history.json path")
    ap.add_argument("--output", type=Path, required=True, help="output PNG path")
    a = ap.parse_args(argv)

    doc = json.loads(a.history.read_text())
    assert doc["schema"] == 1, f"unexpected history schema {doc['schema']}"
    all_rows = doc["rows"]

    # X axis spans every build in the file (so panels share a timeline); the
    # series themselves render only err-bearing rows.
    xpos, xlabels = build_x_axis(all_rows)
    rows = render_only(all_rows)
    by_metric: dict[str, list[dict]] = {}
    for r in rows:
        by_metric.setdefault(r["metric"], []).append(r)

    n = len(FAMILIES)
    ncol = 2
    nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(12.5, 3.1 * nrow),
                             facecolor="white")
    axes = axes.ravel()

    empty: list[str] = []
    for ax, (title, ylabel, series) in zip(axes, FAMILIES):
        is_preview = series[0][1].startswith("imsims")
        n_pts = draw_panel(ax, title, ylabel, series, by_metric, xpos,
                           len(xlabels), is_preview)
        if n_pts == 0:
            empty.append(title)
        ax.set_xticks(range(len(xlabels)))
        ax.set_xticklabels(xlabels, fontsize=7, family="monospace")
    for ax in axes[n:]:
        ax.set_visible(False)

    suptitle = ("ShapePipe status board -- metric trends across builds")
    fig.suptitle(suptitle, fontsize=15, color=INK, x=0.5, y=0.995, ha="center")

    n_render = len(rows)
    n_preview = sum(bool(r.get("preview")) for r in rows)
    caption = (
        f"x = build (date . short commit); y = metric value with 1-sigma error bar; "
        f"shaded band = tolerance where defined. {n_render} of {len(all_rows)} history rows "
        f"rendered -- rows with err = null are stored but not shown (page rule: no number "
        f"without an uncertainty), so the m-bias panel is empty and the 07-06 star scalars "
        f"drop out. The imsims panel ({n_preview} preview points) is flagged "
        f"SKiLLS . preliminary."
    )
    fig.text(0.5, 0.005, caption, ha="center", va="bottom", fontsize=8.5,
             color=MUTED, wrap=True)

    fig.tight_layout(rect=(0, 0.045, 1, 0.975))
    fig.savefig(a.output, dpi=140, facecolor="white", bbox_inches="tight")
    print(f"wrote {a.output}")
    if empty:
        print(f"  empty panels (err=null only, per page rule): {', '.join(empty)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
