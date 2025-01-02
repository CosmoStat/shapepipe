### sqlite_get_ccd_names.py

# Extract names of single-exposure single-HDU (CCD) names
# from sqlite PSF interpolation information file.



import os
import sys
import glob
from tqdm import tqdm

from sqlitedict import SqliteDict


def get_ccds_from_file(fname):
    """Return set of CCD names from sqlite file.
    """

    # Open sqlite file
    psf_vign_cat = SqliteDict(fname)
    
    # Initialise set of CCD names
    ccds = set()

    # Lazy iteration (do not store list in variable).
    # Keys in postage stamp file are object (galaxy) IDs for which
    # the PSF postage stamp is used elsehere,
    #for obj_ID in tqdm(psf_vign_cat.keys(), desc=f"Processing {fname}"):

    keys = psf_vign_cat.keys()
    for obj_ID in tqdm(keys, desc=f"Processing {fname}"):

        # If vignet is not empty, add all used CCDs (=keys) to name set
        if psf_vign_cat[obj_ID] != "empty":
            ccds.update(psf_vign_cat[obj_ID].keys())

    print(f"{len(ccds)} non-empty CCDs found")

    return ccds


def get_all_ccds(filenames):
    """Return set of CCD names from list of sqlite files.
    """

    all_ccds = set()
    for fname in tqdm(filenames, desc="Overall progress"):
        all_ccds.update(get_ccds_from_file(fname))
    return all_ccds


def main(argv):

    filenames = []
    out_path = "ccds_from_interp_psf.txt"

    basedir = "tile_runs"
    iterable = os.scandir(basedir)
    iterable = tqdm(iterable, desc="tile IDs")
    for ID in iterable:
        subdir = f"{ID.path}/output"
        if os.path.exists(subdir):
            #iterable_sub = tqdm(iterable_sub, desc="module_outdirs", leave=False)
            #with iterable_sub as entries:
            with os.scandir(subdir) as entries:
                for entry in entries:
                    if entry.name.startswith("run_sp_tile_PsViSmVi_20"):
                        modoutdir = f"{entry.path}/psfex_interp_runner/output"
                        with os.scandir(modoutdir) as entries_out:
                            for entry_out in entries_out:
                                if entry_out.name.startswith("galaxy_psf"):
                                    filenames.append(entry_out.path)

    print(f"{len(filenames)} galaxy psf sqlite files found")

    all_ccds = get_all_ccds(filenames)
    print(all_ccds)

    with open(out_path, "w") as f_out:
        f_out.writelines(f"{ccd}\n" for ccd in all_ccds)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
    

