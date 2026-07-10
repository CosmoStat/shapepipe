# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib", "seaborn"]
# ///
"""S4 comparison figure: metacal breakdown-grid ablations vs baseline.

Reads the v2 breakdown_grid JSON (see ``run_breakdown_grid.py`` / ``plot_breakdown_grid.py``)
for a baseline run and one or more ablation runs, and draws a 3 (S/N band) x 3 (m1, m2, c1)
grid of panels vs resolution ratio, baseline in dark gray and ablations in husl hues, so a
reader can see at a glance whether any ablation moves the bias trends outside the tolerance
bands. Cells with a null estimator (the degeneracy guard fired) are skipped -- there is no
value to plot. Grids may be partial (a single cell, as in the S4 smoke tests); the script
draws whatever is present. Seaborn is not in the ShapePipe sif, so this runs on the host:

    uv run scripts/python/plot_s4_ablations.py \\
        --baseline <run>/results/baseline/breakdown_grid \\
        --ablations <run>/results/wcscen/breakdown_grid <run>/results/jac/breakdown_grid \\
        --labels "WCS centroid" "Jacobian" \\
        --output <figures>

Add ``--paired <ablation-dir>`` when an ablation's images are bit-identical to baseline
(e.g. the WCS-centroid arm) to also recompute a paired Delta-m1 per cell straight from the
raw metacal records in ``<dir>/cells/*.npz``, which cancels shape noise and gives an error
roughly two orders of magnitude tighter than the independent-run comparison.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter
import seaborn as sns

BANDS = [("high", r"S/N $\approx$ 1000"), ("moderate", r"S/N $\approx$ 50"),
         ("low", r"S/N $\approx$ 15")]
COLUMNS = [("m1", r"$m_1$", 5e-3), ("m2", r"$m_2$", 5e-3), ("c1", r"$c_1$", 1e-3)]
GAMMA = 0.02


def load_grid(grid_dir):
    """{band: {res_ratio: cell}} for a breakdown_grid_v2.json, skipping cells with any
    null estimator -- a null means the degeneracy guard fired and there is nothing to plot."""
    d = json.load(open(Path(grid_dir) / "breakdown_grid_v2.json"))
    by_band = {}
    for c in d["cells"]:
        if c["m1"] is None or c["m2"] is None or c["c1"] is None:
            continue
        by_band.setdefault(c["s2n_band"], {})[c["res_ratio"]] = c
    return by_band, d["meta"]


def series(by_band, band, key):
    """(res_ratio, value, err) for one band and one estimator key, sorted by res_ratio."""
    cells = by_band.get(band, {})
    rows = sorted(cells.values(), key=lambda c: c["res_ratio"])
    r = np.array([c["res_ratio"] for c in rows])
    val = np.array([c[key] for c in rows], dtype=float)
    err = np.array([c[f"{key}_err"] for c in rows], dtype=float)
    return r, val, err


def draw_panel(ax, runs, band, key, tol, pal, dodge_frac=0.015):
    """One (band, key) panel: baseline dark gray, ablations husl-hued, slight x-dodge so
    error bars at the same res_ratio don't overlap."""
    n = len(runs)
    for i, (label, by_band, color) in enumerate(runs):
        r, val, err = series(by_band, band, key)
        if r.size == 0:
            continue
        dodge = (i - (n - 1) / 2) * dodge_frac
        ax.errorbar(r * (1 + dodge), val, yerr=err, color=color, marker="o", ms=6,
                    lw=1.6, capsize=3, elinewidth=1.2, label=label, zorder=3 if i == 0 else 2)
    ax.axhspan(-tol, tol, color="0.75", alpha=0.35, zorder=0)
    ax.axhline(0, color="0.5", lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter())  # dense log-minor labels collide otherwise
    ax.xaxis.set_major_formatter(ScalarFormatter())


def make_figure(baseline, ablations, labels, out_png):
    base_by_band, _ = load_grid(baseline)
    abl_grids = [load_grid(d)[0] for d in ablations]
    pal = sns.color_palette("husl", len(ablations))
    runs = [("baseline", base_by_band, "0.25")] + list(zip(labels, abl_grids, pal))

    sns.set_theme(style="ticks", context="talk")
    fig, axes = plt.subplots(len(BANDS), len(COLUMNS), figsize=(4.4 * len(COLUMNS), 3.6 * len(BANDS)),
                              sharex="col", squeeze=False)

    for row, (band, band_lab) in enumerate(BANDS):
        for col, (key, key_lab, tol) in enumerate(COLUMNS):
            ax = axes[row][col]
            draw_panel(ax, runs, band, key, tol, pal)
            if row == 0:
                ax.set_title(key_lab)
            if col == 0:
                ax.set_ylabel(f"{band_lab}\nvalue")
            if row == len(BANDS) - 1:
                ax.set_xlabel(r"resolution ratio  $r$")

    handles = [Line2D([0], [0], color=color, marker="o", lw=1.6, label=label)
               for label, _, color in runs]
    fig.legend(handles=handles, frameon=False, loc="upper center",
               ncol=min(len(runs), 5), bbox_to_anchor=(0.5, 1.02))
    sns.despine(fig)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    print("wrote", out_png)


def print_ablation_table(baseline, ablations, labels):
    """Per (ablation, band): max |m1|/|m2| over resolved cells (res_ratio >= 0.35), its
    error, and the max |Delta| vs baseline with combined error and significance."""
    base_by_band, _ = load_grid(baseline)
    for label, abl_dir in zip(labels, ablations):
        abl_by_band, _ = load_grid(abl_dir)
        print(f"\n=== {label} ===")
        header = f"{'band':<10}{'stat':<6}{'max|val|':>10}{'err':>10}{'max|Δ|':>10}{'σ_Δ':>10}{'|Δ|/σ':>8}"
        print(header)
        for band, _ in BANDS:
            for key in ("m1", "m2"):
                base_r, base_val, base_err = series(base_by_band, band, key)
                abl_r, abl_val, abl_err = series(abl_by_band, band, key)
                common = sorted(set(base_r[base_r >= 0.35]) & set(abl_r[abl_r >= 0.35]))
                if not common:
                    continue
                base_idx = {r: i for i, r in enumerate(base_r)}
                abl_idx = {r: i for i, r in enumerate(abl_r)}
                i_max = max(range(len(abl_r)), key=lambda i: abs(abl_val[i]) if abl_r[i] >= 0.35 else -1)
                max_val, max_err = abs(abl_val[i_max]), abl_err[i_max]

                deltas = np.array([abl_val[abl_idx[r]] - base_val[base_idx[r]] for r in common])
                sigmas = np.array([np.sqrt(abl_err[abl_idx[r]]**2 + base_err[base_idx[r]]**2)
                                   for r in common])
                j = int(np.argmax(np.abs(deltas)))
                max_delta, sigma_delta = abs(deltas[j]), sigmas[j]
                sig = max_delta / sigma_delta if sigma_delta > 0 else np.inf
                flag = " ***" if sig > 3 else ""
                print(f"{band:<10}{key:<6}{max_val:>10.4f}{max_err:>10.4f}"
                      f"{max_delta:>10.4f}{sigma_delta:>10.4f}{sig:>8.2f}{flag}")


def paired_delta_m1(baseline_dir, ablation_dir, B=1000, seed=0):
    """Paired Delta-m1 per cell: recompute m1 directly from raw metacal records so shape
    noise -- common to both runs, since the images are bit-identical -- cancels. Bootstraps
    over seeds with one shared index resample per replicate across all arrays."""
    base_cells = sorted(Path(baseline_dir, "cells").glob("*.npz"))
    abl_cells = sorted(Path(ablation_dir, "cells").glob("*.npz"))
    base_by_name = {p.name: p for p in base_cells}
    abl_by_name = {p.name: p for p in abl_cells}
    common_names = sorted(set(base_by_name) & set(abl_by_name))
    if not common_names:
        raise ValueError(f"no matching cell files between {baseline_dir} and {ablation_dir}")

    rng = np.random.default_rng(seed)
    results = []
    for name in common_names:
        base = np.load(base_by_name[name])
        abl = np.load(abl_by_name[name])
        ok = (base["g1p__flags"] == 0) & (base["g1p__finite"]) \
            & (base["g1m__flags"] == 0) & (base["g1m__finite"]) \
            & (abl["g1p__flags"] == 0) & (abl["g1p__finite"]) \
            & (abl["g1m__flags"] == 0) & (abl["g1m__finite"])
        if not ok.any():
            raise ValueError(f"no jointly-valid seeds in {name}")

        base_g1p, base_g1m = base["g1p__g1"][ok], base["g1m__g1"][ok]
        base_R11p, base_R11m = base["g1p__R11"][ok], base["g1m__R11"][ok]
        abl_g1p, abl_g1m = abl["g1p__g1"][ok], abl["g1m__g1"][ok]

        Rbar_base = 0.5 * (base_R11p.mean() + base_R11m.mean())
        delta_m1 = ((abl_g1p - base_g1p).mean() - (abl_g1m - base_g1m).mean()) / (2 * GAMMA * Rbar_base)

        n = ok.sum()
        boot = np.empty(B)
        for b in range(B):
            idx = rng.integers(0, n, size=n)
            Rbar_b = 0.5 * (base_R11p[idx].mean() + base_R11m[idx].mean())
            boot[b] = ((abl_g1p[idx] - base_g1p[idx]).mean()
                       - (abl_g1m[idx] - base_g1m[idx]).mean()) / (2 * GAMMA * Rbar_b)
        results.append((name, float(delta_m1), float(boot.std())))
    return results


def print_paired_table(baseline_dir, ablation_dir):
    print(f"\n=== paired Δm1: {ablation_dir} vs {baseline_dir} ===")
    for name, delta_m1, err in paired_delta_m1(baseline_dir, ablation_dir):
        print(f"{name:<24}Δm1 = {delta_m1:+.6f} ± {err:.6f}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--baseline", required=True, help="baseline breakdown_grid dir.")
    p.add_argument("--ablations", required=True, nargs="+", help="ablation breakdown_grid dirs.")
    p.add_argument("--labels", required=True, nargs="+", help="one label per ablation dir.")
    p.add_argument("--output", required=True, help="output figure directory.")
    p.add_argument("--paired", default=None,
                    help="ablation dir (must also be in --ablations) to additionally run "
                         "the paired same-noise Δm1 estimator against --baseline.")
    a = p.parse_args()

    if len(a.ablations) != len(a.labels):
        raise ValueError(f"--ablations has {len(a.ablations)} entries but --labels has "
                          f"{len(a.labels)}; need one label per ablation dir.")

    out_dir = Path(a.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    make_figure(a.baseline, a.ablations, a.labels, out_dir / "s4_ablations.png")
    print_ablation_table(a.baseline, a.ablations, a.labels)

    if a.paired is not None:
        print_paired_table(a.baseline, a.paired)


if __name__ == "__main__":
    main()
