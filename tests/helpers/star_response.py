"""Star shear-response metric — Fabian's R-function, faithfully.

The deconvolution step in metacalibration should leave stars with *no* net
shear response: applied to point sources, the per-object metacal response
``R = dg/dgamma`` must average to zero. A non-zero ``<R>`` means PSF handling
is leaking shape into the response — the very failure being debugged (see the
Fabian handoff, 2026-06-12).

This module reproduces Fabian's exact computation as a library function so it
can be both unit-tested and run on demand:

* ngmix metacal catalogs at
  ``output_*/run_sp_tile_ngmix_Ng1u_*/ngmix_runner/output/ngmix-*.fits``
* HDU map (load-bearing): ``[2]=1p, [1]=1m, [4]=2p, [3]=2m``
* mask ``20 < mag < 26``
* ``R1 = (g1_1p - g1_1m) / 0.02``,  ``R2 = (g2_2p - g2_2m) / 0.02``
* 100-group jackknife for the error on ``<R>``

The defaults match the handoff so a bare call reproduces Fabian's number.
"""

import glob
import os

import numpy as np
from astropy.io import fits


DEFAULT_BASE_DIR = "/home/hervas/n25/SP_simu_fab/SP_1z2z_star_grid/outputs/"
METACAL_STEP = 0.02
MAG_LO, MAG_HI = 20.0, 26.0
N_JACKKNIFE = 100

# HDU index -> metacal type. Load-bearing: this is the on-disk ordering of the
# ngmix_runner output (verified against ngmix-*.fits: [1]=1M [2]=1P [3]=2M [4]=2P).
HDU = {"1p": 2, "1m": 1, "2p": 4, "2m": 3}


def find_ngmix_files(base_dir=DEFAULT_BASE_DIR):
    """Latest ngmix FITS per ``output_*`` tile dir, in Fabian's order.

    Returns ``[(tile_id, path), ...]``. Tiles with no run dir / no FITS are
    skipped silently — matching the handoff's ``continue`` behaviour.
    """
    found = []
    for output_dir in sorted(glob.glob(os.path.join(base_dir, "output_*"))):
        run_dirs = sorted(glob.glob(
            os.path.join(output_dir, "run_sp_tile_ngmix_Ng1u_*")
        ))
        if not run_dirs:
            continue
        ngmix_files = sorted(glob.glob(
            os.path.join(run_dirs[-1], "ngmix_runner/output/ngmix-*.fits")
        ))
        if not ngmix_files:
            continue
        ngmix_file = ngmix_files[-1]
        tile_id = os.path.basename(ngmix_file).replace("ngmix-", "").replace(".fits", "")
        found.append((tile_id, ngmix_file))
    return found


def tile_responses(ngmix_file):
    """Per-object ``(R1, R2)`` for one tile, masked to ``20 < mag < 26``."""
    with fits.open(ngmix_file) as cat:
        c1p, c1m = cat[HDU["1p"]].data, cat[HDU["1m"]].data
        c2p, c2m = cat[HDU["2p"]].data, cat[HDU["2m"]].data
        mask = (c1m["mag"] < MAG_HI) & (c1m["mag"] > MAG_LO)
        R1 = (c1p["g1"][mask] - c1m["g1"][mask]) / METACAL_STEP
        R2 = (c2p["g2"][mask] - c2m["g2"][mask]) / METACAL_STEP
    return R1, R2


def jackknife_error(values, n_groups=N_JACKKNIFE):
    """Delete-d jackknife std of the mean over ``n_groups`` contiguous groups."""
    n = len(values)
    groups = np.array_split(np.arange(n), n_groups)
    jk = np.empty(n_groups)
    for i, idx in enumerate(groups):
        keep = np.ones(n, dtype=bool)
        keep[idx] = False
        jk[i] = np.average(values[keep])
    return np.sqrt((n_groups - 1) / n_groups * np.sum((jk - jk.mean()) ** 2))


def compute_star_response(base_dir=DEFAULT_BASE_DIR, n_jackknife=N_JACKKNIFE):
    """Full star shear-response metric across all tiles under ``base_dir``.

    Returns
    -------
    dict
        ``R1`` / ``R2`` (per-object arrays), ``R1_mean`` / ``R2_mean``,
        ``R1_jk_err`` / ``R2_jk_err``, ``n_obj``, ``n_tiles``, ``base_dir``,
        ``tolerance``. Raises ``FileNotFoundError`` if no FITS are found and
        ``ValueError`` if fewer objects than jackknife groups.
    """
    files = find_ngmix_files(base_dir)
    if not files:
        raise FileNotFoundError(
            f"no ngmix-*.fits under {base_dir}/output_*/run_sp_tile_ngmix_Ng1u_*"
        )

    R1_all = np.array([])
    R2_all = np.array([])
    for _tile_id, ngmix_file in files:
        R1, R2 = tile_responses(ngmix_file)
        R1_all = np.append(R1_all, R1)
        R2_all = np.append(R2_all, R2)

    n_obj = len(R1_all)
    if n_obj < n_jackknife:
        raise ValueError(
            f"only {n_obj} objects; need >= {n_jackknife} for jackknife"
        )

    return {
        "R1": R1_all,
        "R2": R2_all,
        "R1_mean": float(np.average(R1_all)),
        "R2_mean": float(np.average(R2_all)),
        "R1_jk_err": float(jackknife_error(R1_all, n_jackknife)),
        "R2_jk_err": float(jackknife_error(R2_all, n_jackknife)),
        "n_obj": int(n_obj),
        "n_tiles": len(files),
        "base_dir": str(base_dir),
        "tolerance": 0.03,
    }
