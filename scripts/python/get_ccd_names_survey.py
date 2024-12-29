import sys
import numpy as np

from tqdm import tqdm


def main(argv):

    job = 32

    n_CCD = 40

    n_patch = 7
    patches = [f'P{x}' for x in np.arange(n_patch) + 1]

    missing_CCDs = {}
    for patch in patches:
        missing_CCDs[patch] = set()

    # Loop over patches
    for patch in tqdm(patches, desc="read summaries"):

        # Read in files from summary output with missing CCDs.
        # Create and update CCD name set.
        fname = f"{patch}/summary/missing_job_{job}_all.txt"
        with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                missing_CCDs[patch].update([line.rstrip("\n")])

    # Create CCD name lists
    all_CCDs = {}
    used_CCDs = {}
    for patch in patches:
        all_CCDs[patch] = []
        used_CCDs[patch] = set()

    # Read in exposure name files.
    for patch in tqdm(patches, desc="read exposure names"):
        fname = f"{patch}/exp_numbers.txt"
        with open(fname) as f:
            lines = f.readlines()
            for line in lines:
                exp = line.rstrip("\n")
                for idx in range(n_CCD):
                    all_CCDs[patch].append(f"{exp}-{idx}")

        used_CCDs[patch] = set(all_CCDs[patch]) - missing_CCDs[patch]

    # Combine all patches
    used_all = set()
    missing_all = set()
    all_all = []
    n_missing = 0
    n_used = 0
    for patch in patches:
        used_all.update(used_CCDs[patch])
        missing_all.update(missing_CCDs[patch])
        all_all.extend(all_CCDs[patch])
        n_missing += len(missing_CCDs[patch])
        n_used += len(used_CCDs[patch])
    n_all = len(set(all_all))

    print("patch    missing      total       used")
    for patch in patches:
        print(f"{patch:5s} {len(missing_CCDs[patch]):10d} {len(all_CCDs[patch]):10d} {len(used_CCDs[patch]):10d}")

    print(f"{'sum':5s} {n_missing:10d} {len(all_all):10d} {n_used:10d}")
    print(f"{'all':5s} {len(missing_all):10d} {n_all:10d} {len(used_all):10d}")

    print(f"Multiple CCDs across patches = {(len(all_all) - n_all) / len(all_all):.1%}")

    for patch in tqdm(patches, desc="write output files"):
        out_path = f"ccds_with_psf_{patch}.txt"
        with open(out_path, "w") as f_out:
            f_out.writelines(f"{ccd}\n" for ccd in used_CCDs[patch])
    out_path = f"ccds_with_psf_all.txt"
    with open(out_path, "w") as f_out:
        f_out.writelines(f"{ccd}\n" for ccd in used_all)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
