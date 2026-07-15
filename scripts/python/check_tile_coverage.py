#!/usr/bin/env python3
"""Check input-tile coverage for image simulations.

Scans the weight maps of the input simulation tiles and computes the
fraction of nonzero pixels (subsampled). Tiles below the coverage
threshold are written to an exclude list (YAML), which the orchestrating
workflow reads at DAG-build time to drop them from the run. Orchestration
now lives in sp_validation; this script is the standalone coverage check
it calls.

The weight-map location is read from the ShapePipe get-images config
{base}/{sim}/cfis/config_tile_Git_symlink.ini (created by init), for the
first sheared variant (1p2z); coverage is identical across variants.

(Near-)empty tiles occur at the edge of the simulated footprint; they
crash the pipeline downstream (vignetmaker on 0 objects, ngmix on empty
PSF vignettes) and their edge objects bias the m-bias estimate.

Usage:
    check_tile_coverage.py -c config.yaml [-o output.yaml] [-t THRESH]
                           [-s SUBSAMPLE] [-v]
"""

import argparse
import configparser
import os
import re
import sys

import numpy as np
import yaml
from astropy.io import fits


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-c", "--config", required=True, help="run config YAML file")
    p.add_argument(
        "-o", "--output",
        help="output YAML path (default {base}/tile_coverage_grid_{num}.yaml)",
    )
    p.add_argument(
        "-t", "--threshold", type=float,
        help="minimum coverage fraction (default: min_tile_coverage from"
             " config, or 0.5)",
    )
    p.add_argument(
        "-s", "--subsample", type=int, default=20,
        help="read every Nth pixel per axis (default: %(default)s)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="verbose output")
    return p.parse_args()


def get_tile_ids(config):
    """Read the tile ID list the same way the orchestrating workflow does."""
    num = config["num"]
    path = config["tile_IDs"].replace("{num}", str(num))
    with open(path) as f:
        tile_ids = [line.strip() for line in f if line.strip()]

    exclusive = config.get("tile_IDs_exclusive")
    if exclusive and isinstance(exclusive, list):
        tile_ids = exclusive

    return tile_ids


def get_weight_template(config):
    """Weight-map path template from the get-images ShapePipe config.

    Reads {base}/{sim}/cfis/config_tile_Git_symlink.ini of the first
    sheared variant (1p2z) and returns (directory, pattern, extension,
    numbering regex) of the weight entry.
    """
    num = config["num"]
    sims_type = config.get("sims_type", "grid")
    str_type = f"_{sims_type}" if sims_type == "grid" else ""
    sim_dir = os.path.join(config["base"], f"1p2z{str_type}_{num}")
    ini_path = os.path.join(sim_dir, "cfis", "config_tile_Git_symlink.ini")

    cp = configparser.ConfigParser(interpolation=None)
    if not cp.read(ini_path):
        raise FileNotFoundError(
            f"{ini_path} not found; run 'snakemake init_all' first"
        )
    sec = cp["GET_IMAGES_RUNNER"]
    paths = [p.strip() for p in sec["INPUT_PATH"].split(",")]
    patterns = [p.strip() for p in sec["INPUT_FILE_PATTERN"].split(",")]
    exts = [e.strip() for e in sec["INPUT_FILE_EXT"].split(",")]
    numbering = sec["INPUT_NUMBERING"].split(",")[0].strip()

    idx = next(
        (i for i, p in enumerate(patterns) if "weight" in p.lower()), None
    )
    if idx is None:
        raise ValueError(f"no weight entry in INPUT_FILE_PATTERN of {ini_path}")

    weight_dir = paths[idx].replace("$SP_DIR", sim_dir)
    return weight_dir, patterns[idx], exts[idx], numbering


def coverage_fraction(path, subsample):
    """Fraction of nonzero pixels in a weight map, subsampled."""
    with fits.open(path) as hdul:
        data = hdul[0].data[::subsample, ::subsample]
    return float(np.count_nonzero(data) / data.size)


def main():
    args = parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    threshold = args.threshold
    if threshold is None:
        threshold = config.get("min_tile_coverage", 0.5)

    num = config["num"]
    tile_ids = get_tile_ids(config)
    weight_dir, pattern, ext, numbering = get_weight_template(config)
    out_path = args.output or os.path.join(
        config["base"], f"tile_coverage_grid_{num}.yaml"
    )

    print(f"Scanning {len(tile_ids)} tiles in {weight_dir}")
    print(f"Coverage threshold: {threshold}")

    coverage = {}
    exclude = []
    missing = []
    for tile in tile_ids:
        fname = re.sub(numbering, tile.replace(".", "-"), pattern) + ext
        path = os.path.join(weight_dir, fname)
        if not os.path.isfile(path):
            missing.append(tile)
            exclude.append(tile)
            print(f"  {tile}: WEIGHT FILE MISSING -> exclude")
            continue
        frac = coverage_fraction(path, args.subsample)
        coverage[tile] = round(frac, 4)
        if frac < threshold:
            exclude.append(tile)
        if args.verbose or frac < threshold:
            flag = "  <-- exclude" if frac < threshold else ""
            print(f"  {tile}: coverage {frac:.3f}{flag}")

    result = {
        "num": num,
        "min_coverage": threshold,
        "coverage": coverage,
        "missing": missing,
        "exclude": sorted(exclude),
    }
    with open(out_path, "w") as f:
        yaml.dump(result, f, default_flow_style=False, sort_keys=False)

    n_keep = len(tile_ids) - len(exclude)
    print(f"\n{len(exclude)} of {len(tile_ids)} tiles below threshold"
          f" ({n_keep} remain)")
    print(f"Written to {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
