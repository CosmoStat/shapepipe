"""Shear bias from one detector defect against its distance from the object,
under the production defect treatment and under noise fill.

Production interpolates the defects ``interpolable_defects`` selects
(columns, 3-px bleeds, finite bleeds, single pixels) and noise-fills the
rest (4-column clusters, edge bands). The counterfactual noise-fills every defect, by making
``interpolated_defects`` select none. The recovery is the full-matrix
metacal measurement of ``tests/helpers/defect_response``; the figure plots
the worse axis of |m| and |c| with the veto radii and the bounds
|m| < 1%, |c| < 5e-4.

Run from the shapepipe checkout root inside the shapepipe container:
  PYTHONPATH=src:. python scripts/validation/masked_pixels/bias_vs_distance.py \
      scan OUT.json NPROC [KIND,KIND...]
  PYTHONPATH=src:scripts/validation/masked_pixels \
    python scripts/validation/masked_pixels/bias_vs_distance.py \
      plot OUT.png IN.json [IN.json ...]
"""
import json
import sys

import numpy as np

N, CENTRE = 51, 25
DISTANCES = range(3, 17)
SEEDS = range(6)
GALAXIES = [(0.3, 0.7), (0.5, 0.7), (0.7, 0.9)]  # (hlr, PSF FWHM) arcsec
INTERPOLATED = ("column", "bleed", "finite_bleed", "pixel")
KINDS = INTERPOLATED + ("cluster4", "edge")
BOUND_M, BOUND_C = 0.01, 5e-4


def geometry(kind, distance):
    """A defect whose nearest pixel is ``distance`` px from the stamp
    centre (the geometries of tests/science/test_defect_recovery.py)."""
    bad = np.zeros((N, N), dtype=bool)
    near = CENTRE + distance
    if kind == "pixel":
        bad[CENTRE, near] = True
    elif kind == "column":
        bad[:, near] = True
    elif kind == "bleed":
        bad[:, near:near + 3] = True
    elif kind == "finite_bleed":
        bad[CENTRE - 5:CENTRE + 6, near:near + 3] = True
    elif kind == "cluster4":
        bad[:, near:near + 4] = True
    elif kind == "edge":
        bad[:, near:] = True
    return bad


def one(job):
    import time
    kind, fill, distance, hlr, psf, *extra = job
    psf_shear = extra[0] if extra else (0.0, 0.0)
    from tests.helpers.defect_response import defect_response

    row = dict(kind=kind, fill=fill, d=distance, hlr=hlr, psf=psf,
               psf_shear=psf_shear)
    start = time.monotonic()
    try:
        # Use the public option: the historical three-argument monkeypatch
        # no longer matches interpolated_defects' four-argument signature.
        # This helper fits directly, bypassing BOTH epoch cuts deliberately.
        r = defect_response(geometry(kind, distance), hlr=hlr, psf=psf,
                            seeds=SEEDS, psf_shear=psf_shear,
                            options={"defect_fill": "noise" if fill == "noise"
                                     else "interpolate"})
        row.update(r)
    except Exception as e:  # noqa: BLE001
        row["error"] = repr(e)
    row["seconds"] = time.monotonic() - start
    return row


def napoli_jobs():
    """Same galaxy grid, round AND elliptical PSFs, inside the old radii.

    Production narrow defects are interpolated (0..7 inclusive); the
    counterfactual and production wide holes are noise-filled (0..13).
    These are isolated high-S/N fixtures, not a survey calibration.
    """
    return [(k, f, d, h, p, e) for h, p in GALAXIES
            for e in ((0.0, 0.0), (0.05, 0.02)) for k in KINDS
            for f in (("production", "noise") if k in INTERPOLATED
                      else ("production",))
            for d in range(8 if f == "production" and k in INTERPOLATED else 14)]


def napoli_scan(out_path, nproc):
    """One pool/allocation; journal each completed case, including failures."""
    from multiprocessing import get_context
    jobs = napoli_jobs()
    rows = []
    with open(out_path + ".jsonl", "w", buffering=1) as journal:
        with get_context("spawn").Pool(nproc, maxtasksperchild=1) as pool:
            for row in pool.imap_unordered(one, jobs, chunksize=1):
                rows.append(row)
                journal.write(json.dumps(row) + "\n")
                print(f'{len(rows)}/{len(jobs)} {row["kind"]} {row["fill"]} '
                      f'd={row["d"]} {row.get("error", "ok")}', flush=True)
    rows.sort(key=lambda r: (r["hlr"], r["psf_shear"], r["kind"], r["fill"], r["d"]))
    with open(out_path, "w") as stream:
        json.dump(rows, stream, indent=1)


def napoli_plot(out_path, in_path):
    """Signed diagonal m with bootstrap errors and historical veto radii."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = json.load(open(in_path))
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True, sharey=True)
    colors = dict(zip(KINDS, ("#0072B2", "#D55E00", "#CC79A7", "#009E73",
                              "#E69F00", "#333333")))
    for j, (hlr, psf) in enumerate(GALAXIES):
        for i, e in enumerate(((0., 0.), (0.05, 0.02))):
            ax = axs[i, j]
            ax.axhspan(-1, 1, color="#eeeeee")
            ax.axvline(7, color="black", linestyle=":", label="old interpolation radius")
            ax.axvline(13, color="black", linestyle="--", label="old noise radius")
            for kind in KINDS:
                for fill in ("production", "noise"):
                    sel = sorted([r for r in rows if r["hlr"] == hlr and r["kind"] == kind
                                  and r["fill"] == fill and tuple(r["psf_shear"]) == e
                                  and "error" not in r], key=lambda r: r["d"])
                    for axis, marker in ((0, "o"), (1, "s")):
                        ax.errorbar([r["d"] for r in sel],
                                    [100 * r["m"][axis] for r in sel],
                                    yerr=[100 * r["m_err"][axis] for r in sel],
                                    color=colors[kind], linestyle="-" if fill == "production" else "--",
                                    marker=marker, ms=2, lw=.8, alpha=.8,
                                    label=f'{kind} {fill} m{axis+1}{axis+1}' if i == j == 0 else None)
            ax.set_yscale("symlog", linthresh=1)
            ax.set_title(f'HLR {hlr}″, PSF {psf}″, e={e}')
            ax.grid(alpha=.2)
            ax.set_xticks(range(0, 14, 2))
            if j == 0:
                ax.set_ylabel("signed m [%] (symmetric log scale)")
            if i == 1:
                ax.set_xlabel("nearest defect pixel [px]")
    fig.suptitle("No-veto fixtures: 6 paired seeds; bars = bootstrap 1σ; grey = |m| < 1%\n"
                 "Solid: production (wide holes still noise-filled); dashed: noise counterfactual; o: m11, s: m22")
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=7)
    fig.tight_layout(rect=(0, .17, 1, .93))
    fig.savefig(out_path, dpi=160)


def scan(out_path, nproc, kinds=KINDS):
    from multiprocessing import get_context

    jobs = [(k, f, d, h, p) for h, p in GALAXIES for d in DISTANCES
            for k in kinds
            for f in (("production", "noise") if k in INTERPOLATED
                      else ("production",))]
    # spawn: each worker imports ngmix afresh, so the noise-fill patch
    # cannot leak into a production job.
    with get_context("spawn").Pool(nproc, maxtasksperchild=1) as pool:
        rows = pool.map(one, jobs, chunksize=1)
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=1)


def worst(row, q):
    return max(abs(x) for x in row[q])


def smallest_passing(rows):
    """Smallest distance from which every larger distance passes both
    bounds (None if the largest fails)."""
    rows = sorted(rows, key=lambda r: r["d"])
    passing = None
    for r in reversed(rows):
        if "error" in r or worst(r, "m") >= BOUND_M or worst(r, "c") >= BOUND_C:
            break
        passing = r["d"]
    return passing


def plot(out_path, *in_paths):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from figstyle import AQUA, BLUE, INK, INK2, ORANGE, VIOLET
    from shapepipe.modules.ngmix_package.ngmix import (
        EPOCH_CENTRAL_DEFECT_RADIUS as RN,
        EPOCH_INTERPOLATED_DEFECT_RADIUS as RI,
        EPOCH_MASKED_FRACTION_CUT as FRAC_CUT,
    )

    rows = [r for path in in_paths for r in json.load(open(path))]
    grid, fail, bad_ink = "#e4e3df", "#f6e3dc", "#a8462a"
    GOLD = "#c49a1c"
    series = {
        "column": ("bad column", BLUE),
        "bleed": ("3-px bleed", ORANGE),
        "finite_bleed": ("finite 3-px bleed (11 rows)", VIOLET),
        "pixel": ("single pixel", AQUA),
        "edge": ("edge band (noise-filled)", INK),
        "cluster4": ("4-column cluster (noise-filled)", GOLD),
    }
    plt.rcParams.update({"font.size": 10.5, "axes.edgecolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2,
                         "axes.labelcolor": INK})
    fig, axs = plt.subplots(
        2, len(GALAXIES), figsize=(13, 7.4), dpi=150, sharex=True,
        sharey="row", gridspec_kw=dict(hspace=0.1, wspace=0.05))
    for j, (hlr, psf) in enumerate(GALAXIES):
        for i, (q, bound, scale, ylabel) in enumerate([
                ("m", BOUND_M, 100, "worse-axis |m|  [%]"),
                ("c", BOUND_C, 1, "worse-axis |c|")]):
            ax = axs[i, j]
            b = bound * scale
            ax.axhspan(b, 1e4, color=fail, zorder=0, lw=0)
            ax.axhline(b, color=bad_ink, lw=0.9, zorder=1)
            for r, txt, ha, dx in (
                    (RI, f"interpolated\nveto {RI:g} px", "right", -0.2),
                    (RN, f"noise-fill\nveto {RN:g} px", "left", 0.2)):
                ax.axvline(r, color=INK2, ls=":", lw=1.3, zorder=1)
                if i == 0:
                    ax.text(r + dx, 8e2, txt, color=INK2, fontsize=8.5,
                            va="top", ha=ha, linespacing=1.1)
            for kind, (label, col) in series.items():
                for fill, ls, lw, alpha in (("production", "-", 2.0, 1.0),
                                            ("noise", (0, (3, 2)), 1.4, 0.8)):
                    sel = sorted(
                        (r for r in rows if r["kind"] == kind
                         and r["fill"] == fill and r["hlr"] == hlr
                         and r["psf"] == psf and "error" not in r),
                        key=lambda r: r["d"])
                    if not sel:
                        continue
                    d = np.array([r["d"] for r in sel])
                    v = np.array([worst(r, q) * scale for r in sel])
                    z = 3 if fill == "production" else 2
                    ax.plot(d, v, color=col, ls=ls, lw=lw, alpha=alpha,
                            zorder=z)
                    if fill != "production":
                        continue
                    cut = np.array([geometry(kind, x).mean() > FRAC_CUT
                                    for x in d])
                    ax.plot(d[~cut], v[~cut], "o", color=col, ms=3.5,
                            zorder=z)
                    ax.plot(d[cut], v[cut], "o", mfc="white", mec=col,
                            mew=1.1, ms=4, zorder=z)
            ax.set_yscale("log")
            ax.set_ylim((1e-3, 1e3) if q == "m" else (1e-6, 0.5))
            ax.grid(axis="y", color=grid, lw=0.6, which="major")
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.set_xticks(range(4, 17, 2))
            if j == 0:
                ax.set_ylabel(ylabel)
                ax.text(16.5, b * 1.3, "|m| > 1%" if q == "m" else
                        "|c| > 5e-4", color=bad_ink, fontsize=9,
                        ha="right", va="bottom")
            if i == 0:
                ax.set_title(f'{hlr}″ galaxy, {psf}″ PSF', color=INK,
                             fontsize=11, loc="left")
    fig.supxlabel("distance of the defect's nearest pixel from the object "
                  "centre  [px]", fontsize=10.5, color=INK, y=0.045)
    fig.suptitle("Shear bias from one defect: production treatment (solid) "
                 "vs noise fill (dashed); 6 seeds, worse of axes 1, 2",
                 color=INK2, fontsize=10, y=0.965)
    handles = [Line2D([], [], color=c, lw=2, marker="o", ms=3.5, label=l)
               for l, c in series.values()]
    handles += [
        Line2D([], [], color=INK2, lw=2, marker="o", ms=3.5,
               label="production (narrow defects interpolated)"),
        Line2D([], [], color=INK2, lw=1.4, ls=(0, (3, 2)),
               label="noise fill (not used for narrow defects)"),
        Line2D([], [], color=INK2, lw=0, marker="o", mfc="white", mec=INK2,
               ms=4, label=f"epoch dropped by the {FRAC_CUT:.0%} masked-fraction cut"),
    ]
    fig.legend(handles=handles, frameon=False, fontsize=9.5, ncol=3,
               loc="lower center", bbox_to_anchor=(0.5, -0.075))
    fig.savefig(out_path, bbox_inches="tight")


if __name__ == "__main__":
    if sys.argv[1] == "napoli":
        napoli_scan(sys.argv[2], int(sys.argv[3]))
    elif sys.argv[1] == "napoli_plot":
        napoli_plot(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "scan":
        scan(sys.argv[2], int(sys.argv[3]),
             *([tuple(sys.argv[4].split(","))] if len(sys.argv) > 4 else []))
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], *sys.argv[3:])
    elif sys.argv[1] == "summary":
        rows = [r for path in sys.argv[2:] for r in json.load(open(path))]
        for hlr, psf in GALAXIES:
            for kind in KINDS:
                for fill in ("production", "noise"):
                    sel = [r for r in rows if r["kind"] == kind
                           and r["fill"] == fill and r["hlr"] == hlr
                           and r["psf"] == psf]
                    if sel:
                        print(f"{hlr}/{psf} {kind:13s} {fill:10s} "
                              f"passes from {smallest_passing(sel)} px")
