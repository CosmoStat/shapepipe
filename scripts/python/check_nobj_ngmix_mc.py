#!/usr/bin/env python

import os
import glob
from astropy.io import fits

# Base directory for the tile runs
tile_runs_dir = "tile_runs"


def get_naxis2_value(file_path):
    """Get the NAXIS2 value from HDU #1 of a FITS file."""
    try:
        hdu_list = fits.open(file_path)
    except:
        raise OSError(f"problem with opening {file_path}")

    try:
        return hdu_list[1].header.get("NAXIS2", None)
    except:
        raise IndexError(f"problem with {file_path} header")

    hdu_list.close()


# Iterate through directories
for ID in os.listdir(tile_runs_dir):
    id_path = os.path.join(tile_runs_dir, ID, "output")

    # Check if ID has the correct format
    if not os.path.isdir(id_path): # or not ID.replace(".", "", 1).isdigit():
        print(f"{ID} -- path {id_path} not found")
        continue

    IDx = ID.replace(".", "-")

    # Define file paths
    ngmix_file_pattern = os.path.join(
        id_path,
        f"run_sp_Ms_202*/merge_sep_cats_runner/output/ngmix-{IDx}.fits",
    )
    mc_file_pattern = os.path.join(
        id_path,
        f"run_sp_Mc_202*/make_cat_runner/output/final_cat-{IDx}.fits",
    )

    # Resolve file names using glob to handle wildcard
    ngmix_files = [f for f in glob.glob(ngmix_file_pattern)]
    mc_files = [f for f in glob.glob(mc_file_pattern)]

    if not ngmix_files:
        print(f"{ID} -- ngmix cats {ngmix_file_pattern} not found")
        continue
    if not mc_files:
        print(f"{ID} -- mc cats not found {mc_file_pattern}")
        continue

    # Get newest files
    ngmix_file = max(ngmix_files, key=lambda f: os.path.getmtime(f))
    mc_file = max(mc_files, key=lambda f: os.path.getmtime(f))

    # Get NAXIS2 values
    ngmix_naxis2 = get_naxis2_value(ngmix_file)
    mc_naxis2 = get_naxis2_value(mc_file)

    # Look at separate ngmix cats
    sep_ngmix_file_pattern = os.path.join(
        id_path,
        f"run_sp_tile_ngmix_Ng*/ngmix_runner/output/ngmix-{IDx}.fits",
    )
    sep_ngmix_files = [f for f in glob.glob(sep_ngmix_file_pattern)]
    sep_ngmix_file = max(sep_ngmix_files, key=lambda f: os.path.getmtime(f))

    if not sep_ngmix_files:
        raise ValueError(f"No separate ngmix files found for ID {ID}")
    sep_naxis2 = 0
    for sep_ngmix_file in sep_ngmix_files:
        sep_naxis2 += get_naxis2_value(sep_ngmix_file)

    if ngmix_naxis2 is None:
        print(f"{ID} -- no NAXIS2 keyword found in ngmix cat {ngmix_file}")
        continue
    if mc_naxis2 is None:
        print(f"{ID} -- no NAXIS2 keyword found in mc cat")
        continue

    # Check run time of 128, 256, and 512 jobs
    if (
        os.path.getmtime(ngmix_file) > os.path.getmtime(mc_file)
        or os.path.getmtime(mc_file) > os.path.getmtime(sep_ngmix_file)
    ):
        print("{ID} -- mc, sep, and/or merged ngmix file out of order")

    # Compare values and output ID if they differ
    elif ngmix_naxis2 != sep_naxis2:
        print(
            f"{ID} -- ngmix differ between sep {sep_naxis2} and merged"
            + f" {ngmix_naxis2}"
        )

    elif ngmix_naxis2 != mc_naxis2:
        ratio = float(ngmix_naxis2) / float(mc_naxis2)
        if ratio < 0.75:

            print(
                f"{ID} -- NAXIS2 differ {ngmix_naxis2} / {mc_naxis2} ="
                + f" {ratio:.3f} (sum {len(sep_ngmix_files)} sep_ng ="
                + f" {sep_naxis2})"
            )
        else:
            print(f"{ID} -- similar NAXIS2, {ngmix_naxis2} / {mc_naxis2} = {ratio:.3f}")
    else:
        print(f"{ID} -- identcal NAXIS2 {ngmix_naxis2} {mc_naxis2}")
