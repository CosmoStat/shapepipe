#!/usr/bin/env python
"""Monitor m-bias convergence as pipeline_all tiles complete.

For tiles that have finished all ShapePipe jobs (bit 2048) in every sim,
runs merge/extract/calibrate in a monitoring subdirectory, then computes
m1, m2, c1, c2 and reports them alongside the tile count.

Run repeatedly while pipeline_all is in progress: merge is incremental
(adds new tiles, skips tiles already in the HDF5), so each run picks up
newly finished tiles automatically.

Usage (from run directory):
    python ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/monitor_mbias.py [-c config.yaml] [-v]
"""

import sys
import os
import subprocess
import shutil
import yaml
import argparse


CONFIGFILE = "config_image_sims.yaml"
JOB_LAST   = 2048

HOME     = os.path.expanduser("~")
PATH_SH  = f"{HOME}/astro/repositories/github/shapepipe/scripts/sh"
PATH_PY  = f"{HOME}/astro/repositories/github/shapepipe/scripts/python"
PATH_NB  = f"{HOME}/astro/repositories/github/sp_validation/notebooks"
PATH_SPV = f"{HOME}/astro/repositories/github/sp_validation"
SPV_SRC  = f"{PATH_SPV}/src"
SPV_SIF  = "/n17data/mkilbing/sp_validation_im_sims.sif"
APPTAINER = f"{PATH_SH}/apptainer_noslurm.sh"

SPV_PREFIX = [
    APPTAINER,
    "--bind", "/n17data,/n09data,/home",
    "--env", f"PYTHONPATH={SPV_SRC}",
    SPV_SIF,
    "python",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "-c", "--config", default=CONFIGFILE,
        help=f"Snakemake config YAML (default: {CONFIGFILE})"
    )
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args()


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def sims_from_config(cfg):
    num = cfg["num"]
    str_type = f"_{cfg['type']}" if cfg["type"] == "grid" else ""
    return [
        f"1{d1}2{d2}{str_type}_{num}"
        for d1, d2 in zip(["m", "p", "z", "z", "z"],
                           ["z", "z", "m", "p", "z"])
    ]


def done_tiles(grids, sims, tile_ids):
    """Return tiles that have log_job_{tile}_2048.txt in every sim."""
    result = []
    for tile in tile_ids:
        if all(
            os.path.isfile(f"{grids}/{sim}/logs/log_job_{tile}_{JOB_LAST}.txt")
            for sim in sims
        ):
            result.append(tile)
    return result


def setup_sim_dir(mon_sim_dir, real_sim_dir, verbose):
    """Prepare monitoring sim dir with symlinks and copied params.py."""
    os.makedirs(mon_sim_dir, exist_ok=True)

    for name in ("cfis", "config_mask.yaml"):
        link = os.path.join(mon_sim_dir, name)
        if not os.path.exists(link):
            target = os.path.join(real_sim_dir, name)
            if os.path.exists(target):
                os.symlink(target, link)
            elif verbose:
                print(f"  Warning: {target} not found, skipping symlink")

    params_src = os.path.join(real_sim_dir, "params.py")
    params_dst = os.path.join(mon_sim_dir, "params.py")
    if os.path.exists(params_src):
        shutil.copy2(params_src, params_dst)
    elif verbose:
        print(f"  Warning: {params_src} not found")


def run_spv(script_args, cwd, label, verbose):
    """Run a script inside the SPV container; return True on success."""
    cmd = SPV_PREFIX + script_args
    if verbose:
        print(f"    cmd: {' '.join(cmd)}")
        print(f"    cwd: {cwd}")
    result = subprocess.run(
        cmd, cwd=cwd,
        capture_output=(not verbose),
        text=True,
    )
    if result.returncode != 0:
        print(f"  {label} FAILED (exit {result.returncode})")
        if not verbose and result.stderr:
            print(result.stderr[-2000:])
        return False
    return True


def run_merge(mon_sim_dir, grids, sim, verbose):
    return run_spv(
        [
            f"{PATH_PY}/create_final_cat.py",
            "-I",
            "-m", f"final_cat_{sim}.hdf5",
            "-i", grids,
            "-P", sim,
            "-p", f"{grids}/{sim}/cfis/final_cat.param",
            "-o", "n_tiles_monitor.txt",
        ] + (["-v"] if verbose else []),
        cwd=mon_sim_dir,
        label="merge",
        verbose=verbose,
    )


def run_extract(mon_sim_dir, verbose):
    return run_spv(
        [f"{PATH_NB}/extract_info.py"],
        cwd=mon_sim_dir,
        label="extract",
        verbose=verbose,
    )


def run_calibrate(mon_sim_dir, verbose):
    return run_spv(
        [f"{PATH_NB}/calibrate_comprehensive_cat.py", "-s", "calibrate"],
        cwd=mon_sim_dir,
        label="calibrate",
        verbose=verbose,
    )


def run_mbias(mon_dir, num, verbose):
    cfg_path = os.path.join(mon_dir, "image_sims_m_bias_monitor.yaml")
    out_path = os.path.join(mon_dir, "m_bias_monitor.yaml")
    mbias_cfg = {
        "base":             mon_dir,
        "num":              num,
        "catalog_name":     "shape_catalog_cut_ngmix.fits",
        "shear_amplitude":  0.02,
        "match_radius_deg": 0.0002,
        "w_col":            "w_des",
        "n_bootstrap":      500,
        "output_path":      out_path,
    }
    with open(cfg_path, "w") as f:
        yaml.dump(mbias_cfg, f, default_flow_style=False)

    ok = run_spv(
        [f"{PATH_SPV}/scripts/compute_m_bias_image_sims.py", "-c", cfg_path]
        + (["-v"] if verbose else []),
        cwd=mon_dir,
        label="mbias",
        verbose=verbose,
    )
    return ok, out_path


def n_tiles_from_file(path):
    """Read tile count written by create_final_cat.py -o."""
    try:
        with open(path) as f:
            return int(f.read().strip())
    except Exception:
        return None


def main():
    args = parse_args()
    cfg = load_config(args.config)

    base     = cfg["base"]
    num      = cfg["num"]
    tile_ids = cfg.get("tile_IDs", [cfg.get("tile_ID", "254.286")])
    grids    = f"{base}/{cfg['type']}s"
    sims     = sims_from_config(cfg)
    mon_dir  = os.path.join(grids, "monitoring")

    # -- Detect done tiles -----------------------------------------------
    done = done_tiles(grids, sims, tile_ids)
    n_done  = len(done)
    n_total = len(tile_ids)

    print(f"Tiles finished in all sims: {n_done} / {n_total}")
    if done:
        print(f"  {done}")

    if n_done == 0:
        print("No tiles complete yet. Nothing to compute.")
        return 0

    # -- Set up monitoring directory ------------------------------------
    print(f"\nMonitoring dir: {mon_dir}")
    os.makedirs(mon_dir, exist_ok=True)

    # -- Per-sim: merge → extract → calibrate --------------------------
    for sim in sims:
        real_sim_dir = os.path.join(grids, sim)
        mon_sim_dir  = os.path.join(mon_dir, sim)
        print(f"\n{sim}")

        setup_sim_dir(mon_sim_dir, real_sim_dir, args.verbose)

        print("  merge    ...", end="", flush=True)
        if not run_merge(mon_sim_dir, grids, sim, args.verbose):
            return 1
        n_found = n_tiles_from_file(os.path.join(mon_sim_dir, "n_tiles_monitor.txt"))
        print(f" ok  ({n_found} tiles in HDF5)")

        print("  extract  ...", end="", flush=True)
        if not run_extract(mon_sim_dir, args.verbose):
            return 1
        print(" ok")

        print("  calibrate...", end="", flush=True)
        if not run_calibrate(mon_sim_dir, args.verbose):
            return 1
        print(" ok")

    # -- m-bias --------------------------------------------------------
    print("\nComputing m-bias...")
    ok, results_path = run_mbias(mon_dir, num, args.verbose)
    if not ok:
        return 1

    if os.path.isfile(results_path):
        with open(results_path) as f:
            res = yaml.safe_load(f)
        print()
        print(f"  n_tiles = {n_done} / {n_total}")
        print(f"  m1 = {res['m1']:+.4f} ± {res['m1_err']:.4f}")
        print(f"  c1 = {res['c1']:+.4f} ± {res['c1_err']:.4f}")
        print(f"  m2 = {res['m2']:+.4f} ± {res['m2_err']:.4f}")
        print(f"  c2 = {res['c2']:+.4f} ± {res['c2_err']:.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
