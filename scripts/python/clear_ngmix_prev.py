#!/usr/bin/env python3
"""
Clean up previous ngmix run directories for completed tiles.
Removes *_prev directories when the current run has finished.
"""

import argparse
import os
import re
import shutil
import glob
from pathlib import Path
from astropy.io import fits


def get_fits_info(fits_file):
    """Get number of objects (NAXIS2 in HDU 1) and number of HDUs."""
    try:
        with fits.open(fits_file) as hdul:
            n_hdus = len(hdul)
            # Get NAXIS2 from HDU 1 (second HDU, index 1)
            if len(hdul) > 1:
                n_objects = hdul[1].header.get("NAXIS2", "N/A")
            else:
                n_objects = "N/A (no HDU 1)"
            return n_objects, n_hdus
    except Exception as e:
        return f"Error: {e}", "N/A"


def check_finished(log_file):
    """Check if 'finished' appears in the log file."""
    try:
        with open(log_file, "r") as f:
            content = f.read()
            return "finished" in content
    except Exception as e:
        print(f"Warning: Could not read {log_file}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Clean up previous ngmix run directories for completed tiles."
    )
    parser.add_argument(
        "-n", "--dry-run",
        action="store_true",
        help="Show what would be removed without actually removing files"
    )
    args = parser.parse_args()

    if args.dry_run:
        print("dry run mode: No files will be removed\n")

    base_dir = Path("tile_runs")

    # Pattern to match tile IDs like 123.456
    tile_id_pattern = re.compile(r"^\d+\.\d+$")

    # Find all tile ID directories
    if not base_dir.exists():
        print(f"Error: {base_dir} does not exist")
        return

    tile_dirs = [
        d
        for d in base_dir.iterdir()
        if d.is_dir() and tile_id_pattern.match(d.name)
    ]

    print(f"Found {len(tile_dirs)} tile directories")

    removed_count = 0
    finished_count = 0

    for tile_dir in sorted(tile_dirs):
        tile_id = tile_dir.name

        # Construct path to ngmix_runner directory
        ngmix_dir = (
            tile_dir / "output" / "run_sp_tile_ngmix_Ng1u" / "ngmix_runner"
        )

        if not ngmix_dir.exists():
            continue

        # Check for process-*.log files
        log_pattern = str(ngmix_dir / "logs" / "process-*.log")
        log_files = glob.glob(log_pattern)

        if not log_files:
            continue

        # Check if any log file contains "finished"
        is_finished = any(check_finished(log_file) for log_file in log_files)

        if is_finished:
            finished_count += 1

            # Path to the _prev directory to remove
            prev_dir = tile_dir / "output" / "run_sp_tile_ngmix_Ng1u_prev"

            if prev_dir.exists():

                # File info of new
                ngmix_out_dir = ngmix_dir / "output"
                fits_files = list(ngmix_out_dir.glob("ngmix-*.fits"))
                if fits_files:
                    for fits_file in fits_files:
                        n_objects, n_hdus = get_fits_info(fits_file)
                else:
                    n_objects, n_hdus = (-1, -1)

                # File info of prev
                prev_out_dir = prev_dir / "output"
                fits_files = list(prev_out_dir.glob("ngmix-*.fits"))
                if fits_files:
                    for fits_file in fits_files:
                        n_objects_prev, n_hdus_prev = get_fits_info(fits_file)
                else:
                    n_objects_prev, n_hdus_prev = (-1, -1)

                action = "Would remove" if args.dry_run else "Removing"
                print(
                    f"Tile {tile_id}: {action} {prev_dir} (new: {n_objects}"
                    + f" {n_hdus}) (prev: {n_objects_prev} {n_hdus_prev})"
                )

                if n_objects < n_objects_prev:
                    print(f"Warning: {tile_id}")

                if not args.dry_run:
                    try:
                        shutil.rmtree(prev_dir)
                        removed_count += 1
                    except Exception as e:
                        print(f"Error removing {prev_dir}: {e}")
                else:
                    removed_count += 1
            else:
                pass
        else:
            print(f"Tile {tile_id}: Not finished")

    print(f"\nSummary:")
    print(f"  Tiles with finished runs: {finished_count}")
    if args.dry_run:
        print(f"  Previous directories that would be removed: {removed_count}")
    else:
        print(f"  Previous directories removed: {removed_count}")


if __name__ == "__main__":
    main()
