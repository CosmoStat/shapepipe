"""Info script for image-simulation Snakemake pipeline."""

import sys
import os
import math
import yaml
import subprocess


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SNAKEFILE  = os.path.join(SCRIPT_DIR, "Snakefile")
CONFIGFILE = "config.yaml"
BASE_CMD   = f"snakemake -s {SNAKEFILE} --configfile {CONFIGFILE}"


def load_config():
    with open(CONFIGFILE) as f:
        return yaml.safe_load(f)


def load_tile_ids(cfg):
    """Load tile IDs from config, handling both lists and file paths."""
    num = cfg.get("num", "")
    tile_ids = cfg.get("tile_IDs", None)

    if tile_ids is None:
        # Auto-construct from num
        tile_ids = os.path.expanduser(
            f"~/shapepipe/auxdir/CFIS/im_sims_202606/tile_numbers_{num}.txt"
        )

    if isinstance(tile_ids, str):
        # Expand {num} placeholder
        tile_ids = tile_ids.replace("{num}", str(num))
        if os.path.isfile(tile_ids):
            with open(tile_ids) as f:
                tile_ids = [line.strip() for line in f if line.strip()]
        else:
            tile_ids = [tile_ids]

    return tile_ids if isinstance(tile_ids, list) else [tile_ids]


def get_hdf5_tile_count(hdf5_path):
    """Query tile count from HDF5 file using create_final_cat.py -l.

    Returns the count as a string, or "." if file doesn't exist or query fails.
    """
    if not os.path.isfile(hdf5_path):
        return "."

    try:
        HOME = os.path.expanduser("~")
        script = f"{HOME}/astro/repositories/github/shapepipe/scripts/python/create_final_cat.py"
        result = subprocess.run(
            [script, "-I", "-l", "-m", hdf5_path, "-v"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        if result.returncode != 0:
            return "?"

        for line in result.stdout.split('\n'):
            if line.startswith("Total:"):
                parts = line.split()
                if len(parts) >= 2:
                    return parts[1]

        return "?"
    except Exception:
        return "?"


def get_extract_tile_count(sim_path):
    """Read tile count from shape_catalog_comprehensive_ngmix.hdf5 attributes.

    Returns the count as a string, or "." if file doesn't exist or count not found.
    """
    try:
        import h5py
    except ImportError:
        return "?"

    hdf5_path = os.path.join(sim_path, "shape_catalog_comprehensive_ngmix.hdf5")
    if not os.path.isfile(hdf5_path):
        return "."

    try:
        with h5py.File(hdf5_path, 'r') as hf:
            if 'n_tiles' in hf.attrs:
                return str(hf.attrs['n_tiles'])
        return "."
    except Exception:
        return "?"


def monitor(cfg, verbose=0):
    base      = cfg["base"]
    num       = cfg["num"]
    sim_type  = cfg["sims_type"]
    tile_ids  = load_tile_ids(cfg)
    n_tiles   = len(tile_ids)
    job_seq   = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
    job_hex   = [format(int(math.log2(j)), 'X') for j in job_seq]

    str_type = f"_{sim_type}" if sim_type == "grid" else ""
    sims = [
        f"1{d1}2{d2}{str_type}_{num}"
        for d1, d2 in zip(["m","p","z","z"], ["z","z","m","p"])
    ]
    grids = base

    col_w = max(len(s) for s in sims) + 2
    downstream = [
        ("merge",  "final_cat_{sim}.hdf5"),
        ("extr",   "shape_catalog_comprehensive_ngmix.hdf5"),
        ("calib",  "shape_catalog_cut_ngmix.fits"),
        ("diag",   "results/diagnostics.log"),
        ("bias",   "results/m_bias_results.yaml"),
    ]
    ds_w = max(len(s) for s, _ in downstream)

    tile_label = f"({n_tiles} tile{'s' if n_tiles > 1 else ''})"
    if verbose >= 1:
        mid_head = " ".join(f"{h:>2}" for h in job_hex)
    else:
        mid_head = "".join(job_hex)
    bits_w    = len(mid_head)
    header1   = f"  {'sim':<{col_w}}  {'#jobs':>12}  {'2^n':<{bits_w}}  " + \
                "  ".join(f"{s:<{ds_w}}" for s, _ in downstream)
    header2   = f"  {'':^{col_w}}  {'':>12}  {mid_head}"
    print(f"  {tile_label}")
    print(header1)
    print(header2)
    print("-" * len(header1))

    for sim in sims:
        logs_dir = os.path.join(grids, sim, "logs")

        # For each job bit, count how many tiles have the log file
        tile_counts = {}
        for j in job_seq:
            tile_counts[j] = sum(
                1 for t in tile_ids
                if os.path.isfile(os.path.join(logs_dir, f"log_job_{t}_{j}.txt"))
            )

        jobs_sum  = sum(1 for j in job_seq if tile_counts[j] == n_tiles)
        bits_str  = "".join("o" if tile_counts[j] == n_tiles else "." for j in job_seq)

        ds_status = []
        for idx, (name, pattern) in enumerate(downstream):
            path = os.path.join(grids, sim, pattern.format(sim=sim))
            sim_path = os.path.join(grids, sim)
            if name == "merge":
                ds_status.append(get_hdf5_tile_count(path))
            elif name == "extr":
                ds_status.append(get_extract_tile_count(sim_path))
            else:
                ds_status.append("o" if os.path.isfile(path) else "")

        if verbose >= 1:
            counts_str = " ".join(f"{tile_counts[j]:>2}" for j in job_seq)
            row = f"  {sim:<{col_w}}  {jobs_sum:>12}  {counts_str}  " + \
                  "  ".join(f"{s:>{ds_w}}" for s in ds_status)
        else:
            row = f"  {sim:<{col_w}}  {jobs_sum:>12}  {bits_str}  " + \
                  "  ".join(f"{s:>{ds_w}}" for s in ds_status)
        print(row)


def main():
    monitor_mode = "-m" in sys.argv or "--monitor" in sys.argv
    verbose = sys.argv.count("-v") + sys.argv.count("--verbose") + \
              (1 if "-vv" in sys.argv else 0) * 2

    cfg = load_config()

    if monitor_mode:
        monitor(cfg, verbose=verbose)
        return

    base      = cfg["base"]
    tile_ids  = load_tile_ids(cfg)
    sample    = cfg["sample"]
    num       = cfg["num"]
    sim_type  = cfg["sims_type"]
    job_seq   = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

    str_type = f"_{sim_type}" if sim_type == "grid" else ""
    sims = [
        f"1{d1}2{d2}{str_type}_{num}"
        for d1, d2 in zip(["m","p","z","z"], ["z","z","m","p"])
    ]
    grids = base

    print("=" * 70)
    print("  Image simulations pipeline")
    print("=" * 70)
    print()
    print(f"  Config   : {CONFIGFILE}")
    print(f"  Base     : {base}")
    print(f"  Tile IDs : {tile_ids}")
    print(f"  Sample   : {sample}  (mask config mask_v1.X.{sample}_im_sim.yaml)")
    print(f"  Run num  : {num}")
    print(f"  Sims     : {', '.join(sims)}")
    print(f"  Job seq  : {job_seq}")
    print()

    print("# Basic command")
    print()
    print(f"  {BASE_CMD} -j N_CPU <target>")
    print()

    print("# Target order")
    print()
    stages = [
        ("init_all", "initialise run directories"),
        ("pipeline_all", f"run full job sequence {job_seq} in order"),
        ("merge_all", "merge ShapePipe output into final_cat_<sim>.hdf5"),
        ("extract_all", "extract comprehensive shape catalogue (HDF5)"),
        ("calibrate_all", "apply masks and calibrate (= all)"),
        ("diagnostics", "compute and plot diagnostics"),
        ("mbias", "compute multiplicative biases")
    ]
    for target, desc in stages:
        print(f"  {target:<16}  {desc}")
    print()

    print("# Useful options")
    print()
    print(f"  -p                            print shell commands")
    print(f"  --config 'force=1'            pass --force to run_job_sp_canfar_v2.0.bash")
    print(f"  --forcerun pipeline_all       force Snakemake to rerun all pipeline jobs")
    print(f"  --config 'job=J'              target a specific job bit J")
    print(f"  --rerun-incomplete            if errors occur due to previous interrupted run")
    print(f"  pipeline --config 'job=J' --rerun-incomplete check and fix runs")
    print()

    print("# Monitor progress")
    print()
    print(f"  python info_image_sims.py -m")
    print()

    print("# Run a single shapepipe job for one sim")
    print()
    print(f"  {BASE_CMD} -j 1 \\")
    print(f"    {grids}/<sim>/logs/log_job_J.txt")
    print()



if __name__ == "__main__":
    main()
