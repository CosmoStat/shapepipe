#!/usr/bin/env python

import sys
import os
import re
import numpy as np

from tqdm import tqdm


def get_input_IDs_patch(tiles, patch):

    # Read in files from summary output with missing CCDs.
    # Create and update CCD name set.
    fname = f"{patch}/tile_numbers.txt"
    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            ID = line.rstrip("\n")
            tiles.append(
                {
                    "ID": ID,
                    "patch": patch,
                }
            )

    return tiles


def check_zero_weight(tiles, patches):

    for patch in tqdm(patches, desc="check zero weight"):

        fname = f"{patch}/summary/special_job_16_sextractor_runner_5.txt"

        if os.path.exists(fname):

            with open(fname) as f:
                lines = f.readlines()
                for line in lines:
                    m = re.search("process-(\d{3})-(\d{3})\.log", line)
                    if m:
                        ID = f"{m[1]}.{m[2]}"
                        my_tile = next((tile for tile in tiles if tile.get("ID") == ID), None)
                        my_tile["zero_weight"] = 1

    return tiles


def check_exceptions(tiles, patches):

    for patch in tqdm(patches, desc="check processing exceptions"):

        fname = f"{patch}/exceptions.txt"

        if os.path.exists(fname):

            with open(fname) as f:
                lines = f.readlines()
                for line in lines:
                    m = re.search("(\d{3})\.(\d{3})\s+\t+(.*)", line)
                    if m:
                        ID = f"{m[1]}.{m[2]}"
                        my_tile = next((tile for tile in tiles if tile.get("ID") == ID), None)
                        my_tile[m[3]] = 1

    return tiles


def get_final(tiles, patches):

    n_final = {}

    for patch in tqdm(patches, desc="check processing exceptions"):

        n_final[patch] = 0

        fname = f"{patch}/n_tiles_final.txt"
        if os.path.exists(fname):
            with open(fname) as f:
                lines = f.readlines()
            for line in lines:
                num = line.rstrip("\n")
            n_final[patch] = int(num)

    n_final["all"] = sum(n_final.values())

    return n_final


def fill_True(tiles, keys):


    for tile in tiles:
        for key in keys:
            tile[key] = 0

    return tiles


def summary(patches, keys):

    tiles = []

    n_final = get_final(tiles, patches)

    # Loop over patches
    for patch in tqdm(patches, desc="read input tile lists"):
        tiles = get_input_IDs_patch(tiles, patch)

    tiles = fill_True(tiles, keys)

    tiles = check_zero_weight(tiles, patches)

    tiles = check_exceptions(tiles, patches)

    return tiles, n_final

def print_tile(tile, f_out, keys):

    for key in keys:
        print(f"{tile[key]}", file=f_out, end=" ")
    print(file=f_out)


def output_tiles(tiles, path, keys):

    with open(path, "w") as f_out:

        for tile in tiles:
            print_tile(tile, f_out, keys)


def print_summary(tiles, patches, n_final, keys, keys_sum, nform):

    for key in keys_sum:
        print(key.rjust(nform, " "), end="")
    print("tot_miss".rjust(nform, " "))

    for patch in patches:
        tiles_patch = [entry for entry in tiles if entry.get("patch") == patch]
        print_summary_patch(tiles_patch, patch, n_final[patch], keys, nform)
    print_summary_patch(tiles, "all", n_final["all"], keys, nform)

def print_summary_patch(tiles, patch, n_final, keys, nform):

    n_tiles = len(tiles)
    print(patch.rjust(nform, " "), end="")
    print(str(n_tiles).rjust(nform, " "), end="")
    print(str(n_final).rjust(nform, " "), end="")
    n_all = 0
    for key in keys:
        n = sum(v[key] for v in tiles)
        print(str(n).rjust(nform, " "), end="")
        n_all += n
    print(str(n_all).rjust(nform, " "), f"{n_all / n_tiles:10.2%}", end="")
    print()


def main(argv):

    n_patch = 8
    patches = [f'P{x}' for x in np.arange(n_patch) + 1]

    keys = ["zero_weight", "aborted", "large_num_det"]

    nform = 15

    tiles, n_final = summary(patches, keys)

    output_path = "tile_summary.txt"

    keys_ext = ["ID", "patch"] + keys
    output_tiles(tiles, output_path, keys_ext)

    keys_sum = ["patch", "n_tiles", "n_final"] + keys
    print_summary(tiles, patches, n_final, keys, keys_sum, nform)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
