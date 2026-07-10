#!/usr/bin/env python
"""Statistically rigorous resolution x noise breakdown grid (v2).

The Tier-1 companion of ``run_mbias.py``. Where that recipe reports the
multiplicative bias at a single well-resolved, high-S/N operating point, this
one sweeps two axes to *map* where the metacal shear calibration stays within
tolerance and where it breaks down: galaxy size relative to the PSF (the
resolution ratio ``r = gal_hlr / psf_fwhm``) and image noise (three S/N bands).

**What changed from v1.** v1 ran a single-signed arm per axis and reported bare
``m = <g>/<R>/gamma - 1`` with no uncertainty; at moderate/low S/N the per-seed
shear noise gives sigma_m ~ 0.5, so v1's scary S/N trends were *unmeasured*. v2
makes every number carry an honest bootstrap error and cancels measurement noise
by pairing signed arms that share an identical noise realisation:

* **Paired shear response.** Six arms per seed -- (g1+, g1-), (g2+, g2-),
  (psf+, psf-) -- driven so the noise draw and sub-pixel shift are BIT-IDENTICAL
  across arms of a seed (same ``RandomState(seed+1000)`` into ``make_data``;
  galsim ``.shear()`` adds no rng draw). The +/- difference then cancels the
  shared shape noise that dominated v1. ``m1 = (<g1[g1p]> - <g1[g1m]>)/(2 gamma
  Rbar11) - 1`` with ``Rbar11`` the mean R11 over both g1 arms.
* **Additive bias from a PSF-shear pair.** ``c1 = (<g1[psfp]> - <g1[psfm]>)/(2
  Rbar11_psf)``, leakage ``alpha = c1 / 0.05``. c2 (expected 0) is a null check.
* **Bootstrap-over-seeds errors** (B=1000) on every cell estimator.
* **Degeneracy guard.** When ``|Rbar| < 5 se(Rbar)`` the m is set null with
  ``degenerate: true`` rather than dividing by a response consistent with zero.
* **Legacy estimator** (v1's ratio-of-averages) carried alongside for continuity.

Estimator math is exact; the finite-difference response step is hardcoded to
``do_ngmix_metacal``'s mainline value (Spec-Code Invariant: if that line moves,
``STEP`` moves with it).

Runs ONLY inside the ShapePipe apptainer container with the live source tree
bound over the baked copy -- same invocation as ``run_mbias.py`` (two binds: the
live shapepipe src AND the live cs_util, since mainline ngmix.py imports
``cs_util.size``). Cells are embarrassingly parallel (``--jobs``); the full
24-cell production grid is driven by the orchestrator on SLURM, not here.

    apptainer exec \\
      --bind /home,/scratch,/automnt,/n17data,/n23data1,/n09data \\
      --bind /n17data/cdaley/unions/code/shapepipe/src/shapepipe:/app/src/shapepipe \\
      --bind /n17data/cdaley/unions/code/cs_util/cs_util:/app/.venv/lib/python3.12/site-packages/cs_util \\
      /n17data/cdaley/containers/shapepipe_ngmix_v2.0-dev.sif \\
      python scripts/python/run_breakdown_grid.py \\
      --output <run>/results/baseline/breakdown_grid

Emits ``<output>/breakdown_grid_v2.json`` (meta + provenance + per-cell summaries
with all errors) plus a mirror ``breakdown_grid.json`` for the status-page reader,
and per-cell ``<output>/cells/<band>_r<ratio>.npz`` holding the raw per-(seed,arm)
arrays so a cell can be re-analysed without refitting. The column figure is drawn
from the JSON by the host-side ``plot_breakdown_grid.py`` (seaborn, not in the sif).
"""
import argparse
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor

import galsim
import ngmix
import numpy as np
from modopt.math.stats import sigma_mad

from shapepipe.modules.ngmix_package.ngmix import (
    Postage_stamp,
    do_ngmix_metacal,
    get_prior,
)
from shapepipe.testing.simulate import make_data

# Metacal finite-difference step: R = (g1p - g1m) / (2*step). Hardcoded to
# do_ngmix_metacal's mainline value (ngmix.py); the Spec-Code Invariant says if
# that line moves, this must move with it.
STEP = 0.01
# PSF-shear amplitude for the additive-bias (leakage) arms; c1 = dgn / (2*this).
PSF_SHEAR = 0.05
# The three noise levels map to nominal S/N ~ 1000 / 50 / 15, calibrated at the
# grid midpoint (r = 0.5). S/N drifts with galaxy size, so mean_s2n is recorded
# per cell rather than trusted from the label.
NOISE = {"high": 0.08, "moderate": 1.6, "low": 5.0}
# Default seeds per band: moderate/low need many more to beat their shear noise
# down even after pairing. CLI-overridable.
N_SEEDS_DEFAULT = {"high": 100, "moderate": 400, "low": 1600}
# The six signed arms: (shear, psf_shear). g1/g2 probe multiplicative response;
# psf probes additive leakage. Arms of one seed share a noise realisation.
ARMS = {
    "g1p": ((0.02, 0.0), (0.0, 0.0)), "g1m": ((-0.02, 0.0), (0.0, 0.0)),
    "g2p": ((0.0, 0.02), (0.0, 0.0)), "g2m": ((0.0, -0.02), (0.0, 0.0)),
    "psfp": ((0.0, 0.0), (PSF_SHEAR, 0.0)), "psfm": ((0.0, 0.0), (-PSF_SHEAR, 0.0)),
}
# The per-(seed,arm) scalars persisted to the .npz for refit-free re-analysis.
RECORD_KEYS = ["g1", "g2", "R11", "R12", "R21", "R22", "s2n", "T", "T_psf_orig",
               "flags", "finite", "sigma_mad_est", "noise_true"]


def one_stamp(noise, gal_hlr, psf_fwhm, img_size, n_epochs, shear, psf_shear,
              seed, true_noise):
    """One injected-truth realisation -> per-arm record dict.

    Uses ``RandomState(seed)`` for the fit rng and ``RandomState(seed+1000)`` for
    make_data, so every arm of a given seed sees the IDENTICAL noise draw and
    sub-pixel shift (galsim ``.shear()`` consumes no rng) -- the property the
    paired estimator rests on. ``true_noise`` routes the per-pixel true-inverse-
    variance path (see ``--true-noise``)."""
    rng = np.random.RandomState(seed)
    prior = get_prior(0.1857, rng)
    gals, psfs, _, weights, flags, jacobs = make_data(
        rng=np.random.RandomState(seed + 1000), shear=shear, noise=noise,
        n_epochs=n_epochs, gal_hlr=gal_hlr, psf_fwhm=psf_fwhm, img_size=img_size,
        psf_shear=psf_shear,
    )
    stamp = Postage_stamp(bkg_sub=False, megacam_flip=False)
    stamp.gals, stamp.psfs, stamp.weights, stamp.flags, stamp.jacobs = (
        gals, psfs, weights, flags, jacobs,
    )
    if true_noise:
        # True per-pixel sigma into bkg_rms -> prepare_ngmix_weights takes the
        # inverse-variance path and metacal fixnoise gets the true sigma.
        stamp.bkg_rms = [np.full((img_size, img_size), noise) for _ in gals]
    res = do_ngmix_metacal(stamp, prior, 1.0, rng)
    d, orig = res.resdict, res.orig
    ns = d["noshear"]
    finite = bool(np.all(np.isfinite(ns["g"])))
    return dict(
        g1=ns["g"][0] if finite else np.nan, g2=ns["g"][1] if finite else np.nan,
        R11=response(d, 1, 1), R12=response(d, 1, 2),
        R21=response(d, 2, 1), R22=response(d, 2, 2),
        s2n=ns.get("s2n", np.nan), T=ns.get("T", np.nan),
        T_psf_orig=orig["T_psf"], flags=int(ns["flags"]), finite=finite,
        # sigma_mad on the epoch-0 galaxy image is the noise the estimator would
        # infer from the data; its ratio to noise_true is the contamination
        # diagnostic (>1 when object flux inflates the MAD, the moderate/low path).
        sigma_mad_est=float(sigma_mad(gals[0])), noise_true=float(noise),
    )


def response(d, i, j):
    """R_ij = (g_i[<axis j>p] - g_i[<axis j>m]) / (2*step)."""
    p, m = {1: "1p", 2: "2p"}[j], {1: "1m", 2: "2m"}[j]
    gp, gm = d[p]["g"], d[m]["g"]
    if not (np.all(np.isfinite(gp)) and np.all(np.isfinite(gm))):
        return np.nan
    return (gp[i - 1] - gm[i - 1]) / (2 * STEP)


BOOTSTRAP_B = 1000


def _boot_idx(n, B, rng):
    """One seed-index array per bootstrap replicate.

    Every quantity inside a replicate must be computed from the SAME resampled
    seed set: resampling each arm independently silently destroys the ±γ
    pairing and inflates the reported error to the unpaired level (the bug this
    replaces read paired σ_m = 0.23 where the per-seed pairs give 0.009)."""
    return rng.integers(0, n, size=(B, n))


def harvest_cell(args):
    """Six-arm paired estimators for one (resolution, noise) cell.

    ``args`` is a flat tuple so the cell is a single ProcessPoolExecutor task.
    Returns (summary dict, raw per-(seed,arm) arrays for the .npz)."""
    (cell_index, res_ratio, band, noise, gal_hlr, psf_fwhm, img_size, n_epochs,
     gamma, n_seeds, true_noise) = args
    # Seeds independent per cell, never recycled across cells.
    seeds = [10_000_000 + cell_index * 100_000 + s for s in range(n_seeds)]

    # Per-arm columns of per-seed records (aligned by seed index).
    cols = {arm: {k: [] for k in RECORD_KEYS} for arm in ARMS}
    for seed in seeds:
        for arm, (shear, psf_shear) in ARMS.items():
            rec = one_stamp(noise, gal_hlr, psf_fwhm, img_size, n_epochs,
                            shear, psf_shear, seed, true_noise)
            for k in RECORD_KEYS:
                cols[arm][k].append(rec[k])
    arr = {arm: {k: np.asarray(v) for k, v in cols[arm].items()} for arm in ARMS}

    def ok(arm):
        """Boolean mask: this arm's fit succeeded (flags==0 and finite g)."""
        return (arr[arm]["flags"] == 0) & arr[arm]["finite"].astype(bool)

    def pair_mask(a, b):
        """Both arms of the pair usable for this seed (pair-drop policy)."""
        return ok(a) & ok(b)

    fail_frac = {arm: float(np.mean(~ok(arm))) for arm in ARMS}
    rng = np.random.default_rng(cell_index)  # deterministic bootstrap per cell

    summary = dict(
        cell_index=cell_index, res_ratio=float(res_ratio), s2n_band=band,
        gal_hlr=float(gal_hlr), psf_fwhm=float(psf_fwhm), noise=float(noise),
        n_seeds=int(n_seeds), fail_frac=fail_frac,
    )

    # --- multiplicative bias, paired, per axis ---
    for axis, (ap, am, ri, comp) in {
        1: ("g1p", "g1m", "R11", "g1"), 2: ("g2p", "g2m", "R22", "g2"),
    }.items():
        mask = pair_mask(ap, am)
        npair = int(mask.sum())
        gp, gm = arr[ap][comp][mask], arr[am][comp][mask]
        Rp, Rm = arr[ap][ri][mask], arr[am][ri][mask]
        Rbar = float(0.5 * (Rp.mean() + Rm.mean()))
        # Joint bootstrap: one seed resample per replicate, all arms together —
        # the pairing (identical noise across ±γ) must survive the resampling.
        idx = _boot_idx(npair, BOOTSTRAP_B, rng)
        bR = 0.5 * (Rp[idx].mean(axis=1) + Rm[idx].mean(axis=1))
        seR = float(np.std(bR))
        degenerate = abs(Rbar) < 5 * seR
        bm_boot = (gp[idx].mean(axis=1) - gm[idx].mean(axis=1)) / (2 * gamma * bR) - 1.0
        m = None if degenerate else float((np.mean(gp) - np.mean(gm)) / (2 * gamma * Rbar) - 1.0)
        m_err = None if degenerate else float(np.std(bm_boot))
        summary[f"m{axis}"] = m
        summary[f"m{axis}_err"] = m_err
        summary[f"m{axis}_degenerate"] = bool(degenerate)
        summary[f"m{axis}_n_pairs"] = npair
        summary[f"Rbar{axis}{axis}"] = Rbar
        summary[f"Rbar{axis}{axis}_err"] = seR

    # --- additive bias c1/c2 from the PSF-shear pair, leakage alpha ---
    mask = pair_mask("psfp", "psfm")
    npair = int(mask.sum())
    Rpsf_p, Rpsf_m = arr["psfp"]["R11"][mask], arr["psfm"]["R11"][mask]
    Rbar_psf = float(0.5 * (Rpsf_p.mean() + Rpsf_m.mean()))
    idx_c = _boot_idx(npair, BOOTSTRAP_B, rng)
    bRpsf = 0.5 * (Rpsf_p[idx_c].mean(axis=1) + Rpsf_m[idx_c].mean(axis=1))
    seR_psf = float(np.std(bRpsf))
    degenerate_c = abs(Rbar_psf) < 5 * seR_psf
    for axis, comp in ((1, "g1"), (2, "g2")):
        gp, gm = arr["psfp"][comp][mask], arr["psfm"][comp][mask]
        bc = (gp[idx_c].mean(axis=1) - gm[idx_c].mean(axis=1)) / (2 * bRpsf)
        c = None if degenerate_c else float((np.mean(gp) - np.mean(gm)) / (2 * Rbar_psf))
        c_err = None if degenerate_c else float(np.std(bc))
        summary[f"c{axis}"] = c
        summary[f"c{axis}_err"] = c_err
    summary["c_degenerate"] = bool(degenerate_c)
    summary["c_n_pairs"] = npair
    summary["Rbar_psf"] = Rbar_psf
    summary["Rbar_psf_err"] = seR_psf
    summary["alpha"] = None if summary["c1"] is None else float(summary["c1"] / PSF_SHEAR)
    summary["alpha_err"] = None if summary["c1_err"] is None else float(summary["c1_err"] / PSF_SHEAR)

    # --- legacy ratio-of-averages (v1 continuity): m = <g[+arm]>/<R[+arm]>/gamma - 1 ---
    for axis, (ap, ri, comp) in {1: ("g1p", "R11", "g1"), 2: ("g2p", "R22", "g2")}.items():
        m = ok(ap)
        g, R = arr[ap][comp][m], arr[ap][ri][m]
        idx_l = _boot_idx(len(g), BOOTSTRAP_B, rng)  # joint: g and R share the seed resample
        summary[f"m{axis}_legacy"] = float(np.mean(g) / np.mean(R) / gamma - 1.0)
        summary[f"m{axis}_legacy_err"] = float(
            np.std(g[idx_l].mean(axis=1) / R[idx_l].mean(axis=1) / gamma - 1.0))

    # --- diagnostics ---
    g1ok = ok("g1p")
    summary["mean_s2n"] = float(np.nanmean(arr["g1p"]["s2n"][g1ok]))
    summary["mean_T_ratio"] = float(np.nanmean(
        arr["g1p"]["T"][g1ok] / arr["g1p"]["T_psf_orig"][g1ok]))
    # contamination diagnostic: >1 means the data-inferred sigma exceeds the true
    # noise because object flux inflates the whole-stamp MAD.
    summary["mean_sigma_ratio"] = float(np.mean(
        arr["g1p"]["sigma_mad_est"] / arr["g1p"]["noise_true"]))

    return summary, {f"{arm}__{k}": arr[arm][k] for arm in ARMS for k in RECORD_KEYS}


def _git_describe():
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--always", "--dirty"],
            cwd=os.path.dirname(os.path.abspath(__file__)), text=True,
        ).strip()
    except Exception:
        return "unknown"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output", required=True, help="lc output directory.")
    p.add_argument("--gamma", type=float, default=0.02, help="injected shear.")
    p.add_argument("--psf-fwhm", type=float, default=0.55, help="Moffat PSF FWHM (arcsec).")
    p.add_argument("--img-size", type=int, default=51)
    p.add_argument("--n-epochs", type=int, default=1)
    p.add_argument("--n-seeds", type=json.loads, default=None,
                   help='JSON per-band seed counts, e.g. \'{"high":100,"moderate":400,"low":1600}\'.')
    p.add_argument("--n-res", type=int, default=8, help="resolution-ratio grid points.")
    p.add_argument("--res-lo", type=float, default=0.15)
    p.add_argument("--res-hi", type=float, default=1.2)
    p.add_argument("--bands", nargs="+", default=list(NOISE), choices=list(NOISE),
                   help="subset of S/N bands to run (default all three).")
    p.add_argument("--res-index", type=int, nargs="+", default=None,
                   help="subset of resolution-grid indices to run (default all).")
    p.add_argument("--true-noise", action="store_true",
                   help="route the per-pixel true-inverse-variance weight path.")
    p.add_argument("--jobs", type=int, default=4, help="parallel cells (ProcessPoolExecutor).")
    a = p.parse_args()

    n_seeds = a.n_seeds or N_SEEDS_DEFAULT
    res_grid = np.geomspace(a.res_lo, a.res_hi, a.n_res)
    res_idx = a.res_index if a.res_index is not None else list(range(a.n_res))
    # cell_index enumerates the FULL 8x3 grid in (res, band) order so seeds are
    # stable regardless of which subset a given run selects.
    all_cells = [(ri, bi) for ri in range(a.n_res) for bi, _ in enumerate(NOISE)]
    cell_of = {(ri, bi): k for k, (ri, bi) in enumerate(all_cells)}
    bands_ordered = list(NOISE)

    tasks = []
    for ri in res_idx:
        for band in a.bands:
            bi = bands_ordered.index(band)
            cell_index = cell_of[(ri, bi)]
            gal_hlr = float(res_grid[ri] * a.psf_fwhm)
            tasks.append((cell_index, res_grid[ri], band, NOISE[band], gal_hlr,
                          a.psf_fwhm, a.img_size, a.n_epochs, a.gamma,
                          n_seeds[band], a.true_noise))

    os.makedirs(os.path.join(a.output, "cells"), exist_ok=True)
    cells, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for summary, raw in ex.map(harvest_cell, tasks):
            cells.append(summary)
            np.savez_compressed(
                os.path.join(a.output, "cells",
                             f"{summary['s2n_band']}_r{summary['res_ratio']:.3f}.npz"),
                **raw)
            m1 = "  null" if summary["m1"] is None else f"{summary['m1']:+.4f}"
            m1e = "     " if summary["m1_err"] is None else f"{summary['m1_err']:.4f}"
            c1 = "  null" if summary["c1"] is None else f"{summary['c1']:+.5f}"
            print(f"r={summary['res_ratio']:.2f} {summary['s2n_band']:9s} "
                  f"s2n={summary['mean_s2n']:7.1f} m1={m1}+-{m1e} c1={c1} "
                  f"Rbar11={summary['Rbar11']:.3f} "
                  f"deg={summary['m1_degenerate']} fail={summary['fail_frac']['g1p']:.2f}",
                  flush=True)

    out = dict(
        meta=dict(
            version=2, step=STEP, gamma=a.gamma, psf_shear=PSF_SHEAR,
            psf_fwhm=a.psf_fwhm, img_size=a.img_size, n_epochs=a.n_epochs,
            n_seeds=n_seeds, res_grid=res_grid.tolist(), noise_map=NOISE,
            bands=a.bands, res_index=res_idx, true_noise=a.true_noise,
            bootstrap_B=1000, degeneracy_k=5, m_tol_resolved=5e-3,
            runtime_s=time.time() - t0, args=vars(a),
        ),
        provenance=dict(
            ngmix_version=ngmix.__version__, galsim_version=galsim.__version__,
            shapepipe_source=os.path.dirname(__import__("shapepipe").__file__),
            git_describe=_git_describe(),
        ),
        cells=cells,
    )
    # v2 JSON is primary; mirror to breakdown_grid.json so the status-page reader
    # (build_status.py loads that exact name) stays green with no code change.
    for fname in ("breakdown_grid_v2.json", "breakdown_grid.json"):
        with open(os.path.join(a.output, fname), "w") as fh:
            json.dump(out, fh, indent=2)
    print(f"DONE {len(cells)} cells in {time.time() - t0:.1f}s "
          f"-> {a.output}/breakdown_grid_v2.json (+ breakdown_grid.json mirror)")


if __name__ == "__main__":
    main()
