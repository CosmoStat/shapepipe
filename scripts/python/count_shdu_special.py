#!/usr/bin/env python

import os
import re
import sys

ID = sys.argv[1]
print(f"Checking tile ID {ID}...")

# Base directory for the tile runs
tile_runs_dir = "tile_runs"

sum_dir = "summary"
specials = ["special_job_32_setools_runner_0", "special_job_32_setools_runner_1"]


shdu_IDs = []
print("Going through links...")
#3for ID in os.listdir(tile_runs_dir):
if True:
    id_path = os.path.join(tile_runs_dir, ID, "output")

    for filename in os.listdir(id_path):
        if (
            filename.startswith("run_sp_exp_SxSePsf")
            and os.path.islink(os.path.join(id_path, filename))
        ):
            symlink_path = os.path.join(id_path, filename)
            target_path = os.readlink(symlink_path)

            pattern = r"/exp_runs/([^/]+)/output/"

            match = re.search(pattern, target_path)
            if match:
                shdu_IDs.append(match.group(1))
            else:
                print(f"{filename} -> No match in target path: {target_path}")

if len(shdu_IDs) == 0:
    print("No shdus found, exiting")
    sys.exit(1)
else:
    print(f"Found {len(shdu_IDs)} links in dir {id_path}")


print("Going through special files...")
count = 0
for special in specials:
    print(f"special {special}...")
    path = f"{sum_dir}/{special}.txt"
    with open(path, "r") as f:
        content = f.read()
        for shdu_ID in shdu_IDs:
            if re.search(shdu_ID, content):
                count += 1

print(f"Found {count} / {len(shdu_IDs)} = {count / len(shdu_IDs):.2%}")
