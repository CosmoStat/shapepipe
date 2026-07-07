#!/usr/bin/env python
"""Runnable resolution x noise breakdown grid for the ShapePipe digital twin.

The Tier-1 companion of ``run_mbias.py``. Where that recipe reports the
multiplicative bias at a single well-resolved, high-S/N operating point, this
one sweeps two axes to *map* where the metacal shear calibration stays within
tolerance and where it breaks down: galaxy size relative to the PSF (the
resolution ratio ``r = gal_hlr / psf_fwhm``) and image noise (three S/N bands).

Each grid cell is one self-contained in-memory GalSim simulation: an exponential
galaxy convolved with a Moffat (beta = 2.5) PSF, sampled at 0.1857" pixels, with
a known shear injected. The stamp is passed through the *live* pipeline shear
estimator (``do_ngmix_metacal``) -- the same code the real survey uses -- and the
recovered shear, the metacal response ``R``, and the additive bias ``c`` are read
back. Averaging over ``n_seeds`` per cell gives ``m``, ``R11`` and ``c1`` for that
cell. Estimator: *ratio of averages* -- ``m = <g_noshear> / <R> / gamma - 1``,
``c1 = <g_noshear> / <R>`` -- the same shape estimator the twin grounds against.

Runs ONLY inside the ShapePipe apptainer container with the live source tree
bound over the baked copy -- same invocation as ``run_mbias.py`` (two binds: the
live shapepipe src AND the live cs_util, since mainline ngmix.py imports
``cs_util.size``). ~2 min on a candide login node for the default 8x3 grid,
16 seeds/cell:

    apptainer exec \\
      --bind /home,/scratch,/automnt,/n17data,/n23data1,/n09data \\
      --bind /n17data/cdaley/unions/code/shapepipe/src/shapepipe:/app/src/shapepipe \\
      --bind /n17data/cdaley/unions/code/cs_util/cs_util:/app/.venv/lib/python3.12/site-packages/cs_util \\
      /n17data/cdaley/containers/shapepipe_ngmix_v2.0-dev.sif \\
      python scripts/python/run_breakdown_grid.py \\
      --output <run>/results/baseline/breakdown_grid

Emits ``<output>/breakdown_grid.json``; the column figure is drawn from that JSON
by the host-side ``plot_breakdown_grid.py`` (seaborn, not in the sif).
"""
import argparse
import json
import os
import time

import galsim
import ngmix
import numpy as np

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
# The three noise levels map to nominal S/N ~ 1000 / 50 / 15, calibrated at the
# grid midpoint (r = 0.5). S/N drifts with galaxy size, so mean_s2n is recorded
# per cell rather than trusted from the label.
NOISE = {"high": 0.08, "moderate": 1.6, "low": 5.0}


def one_stamp(noise, gal_hlr, psf_fwhm, img_size, n_epochs, shear, psf_shear, seed):
    """One injected-truth realisation -> (metacal resdict, original-fit dict)."""
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
    res = do_ngmix_metacal(stamp, prior, 1.0, rng)
    return res.resdict, res.orig


def response(d, i, j):
    """R_ij = (g_i[<axis j>p] - g_i[<axis j>m]) / (2*step)."""
    p, m = {1: "1p", 2: "2p"}[j], {1: "1m", 2: "2m"}[j]
    return (d[p]["g"][i - 1] - d[m]["g"][i - 1]) / (2 * STEP)


def harvest_cell(noise, gal_hlr, psf_fwhm, img_size, n_epochs, gamma, seeds):
    """Ratio-of-averages m, R11, c1 for one (resolution, noise) cell.

    Three arms per seed: a g1-shear arm (m1, R11), a g2-shear arm (m2, R22), and
    a PSF-shear-only arm (the additive-bias null c1). Failed fits (nonzero flags
    or non-finite shear) are dropped and counted in fail_frac."""
    g1ns, g2ns, R11, R22, c1ns, cR, s2n, trat = ([] for _ in range(8))
    nfail = ntot = 0
    arms = {
        "g1": ((gamma, 0.0), (0.0, 0.0)),
        "g2": ((0.0, gamma), (0.0, 0.0)),
        "psf": ((0.0, 0.0), (0.05, 0.0)),
    }
    for seed in seeds:
        for arm, (shear, psf_shear) in arms.items():
            ntot += 1
            d, orig = one_stamp(noise, gal_hlr, psf_fwhm, img_size, n_epochs,
                                shear, psf_shear, seed)
            ns = d["noshear"]
            if ns["flags"] != 0 or not np.all(np.isfinite(ns["g"])):
                nfail += 1
                continue
            if arm == "g1":
                g1ns.append(ns["g"][0]); R11.append(response(d, 1, 1))
                s2n.append(ns["s2n"]); trat.append(ns["T"] / orig["T_psf"])
            elif arm == "g2":
                g2ns.append(ns["g"][1]); R22.append(response(d, 2, 2))
            else:
                c1ns.append(ns["g"][0]); cR.append(response(d, 1, 1))
    mR11, mR22 = np.mean(R11), np.mean(R22)
    return dict(
        gal_hlr=gal_hlr, psf_fwhm=psf_fwhm, noise=noise,
        m1=float(np.mean(g1ns) / mR11 / gamma - 1.0),
        m2=float(np.mean(g2ns) / mR22 / gamma - 1.0),
        c1=float(np.mean(c1ns) / np.mean(cR)),
        R11=float(mR11), R22=float(mR22),
        mean_s2n=float(np.mean(s2n)), mean_T_ratio=float(np.mean(trat)),
        fail_frac=float(nfail / ntot), n_seeds=len(seeds),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output", required=True, help="lc output directory.")
    p.add_argument("--gamma", type=float, default=0.02, help="injected shear.")
    p.add_argument("--psf-fwhm", type=float, default=0.55, help="Moffat PSF FWHM (arcsec).")
    p.add_argument("--img-size", type=int, default=51)
    p.add_argument("--n-epochs", type=int, default=1)
    p.add_argument("--n-seeds", type=int, default=16, help="realisations averaged per cell.")
    p.add_argument("--n-res", type=int, default=8, help="resolution-ratio grid points.")
    p.add_argument("--res-lo", type=float, default=0.15)
    p.add_argument("--res-hi", type=float, default=1.2)
    a = p.parse_args()

    seeds = list(range(a.n_seeds))
    res_grid = np.geomspace(a.res_lo, a.res_hi, a.n_res)
    cells, t0 = [], time.time()
    for r_ratio in res_grid:
        gal_hlr = float(r_ratio * a.psf_fwhm)
        for band, noise in NOISE.items():
            cell = harvest_cell(noise, gal_hlr, a.psf_fwhm, a.img_size,
                                a.n_epochs, a.gamma, seeds)
            cell.update(res_ratio=float(r_ratio), s2n_band=band)
            cells.append(cell)
            print(f"r={r_ratio:.2f} {band:9s} s2n={cell['mean_s2n']:7.1f} "
                  f"m1={cell['m1']:+.4f} m2={cell['m2']:+.4f} "
                  f"c1={cell['c1']:+.5f} R11={cell['R11']:.3f} "
                  f"fail={cell['fail_frac']:.2f}", flush=True)

    out = dict(
        meta=dict(
            step=STEP, gamma=a.gamma, psf_fwhm=a.psf_fwhm, img_size=a.img_size,
            n_epochs=a.n_epochs, n_seeds=a.n_seeds, seeds=seeds,
            res_grid=res_grid.tolist(), noise_map=NOISE,
            m_tol_resolved=5e-3, runtime_s=time.time() - t0,
        ),
        provenance=dict(
            ngmix_version=ngmix.__version__, galsim_version=galsim.__version__,
            shapepipe_source=os.path.dirname(__import__("shapepipe").__file__),
        ),
        cells=cells,
    )
    os.makedirs(a.output, exist_ok=True)
    with open(os.path.join(a.output, "breakdown_grid.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"DONE {len(cells)} cells in {time.time() - t0:.1f}s "
          f"-> {a.output}/breakdown_grid.json")


if __name__ == "__main__":
    main()
