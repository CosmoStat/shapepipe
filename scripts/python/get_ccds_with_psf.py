#!/usr/bin/env python3

"""GET_CCDS_WITH_PSF

Obtain list of CCDs (single-exposure single-HDU files) for which valid PSF information
is available. This can serve to create a footprint coverage mask.

Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""


import sys
import numpy as np

from shapepipe.utilities import summary


def get_lines(fname):
    """Get Lines.

    Return list of lines read from a text file.

    Parameters
    ----------
    fname: str
        input file name

    Returns:
    list
        IDs

    """
    IDs = []
    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            IDs.append(line.rstrip())

    return IDs


def get_exp_shdu_missing(patches):
    """Get Exp Shdu Missing.

    Returns set of missing CCDs (single-exposure single-HDU IDs) from a list of patches.

    Parameters
    ----------
    patches: list
        input patches

    Returns
    -------
    set
       missing CCD IDs

    """
    exp_shdu_missing_all = set()
    
    for patch in patches:

        path_exp_shdu_missing = f"{patch}/summary/missing_job_32_all.txt"
        exp_shdu_missing = get_lines(path_exp_shdu_missing)

        print(f"Patch {patch}: Found {len(exp_shdu_missing)} missing ccds", end="; ")

        exp_shdu_missing_all.update(exp_shdu_missing)

        print(f"cumulative {len(exp_shdu_missing_all)} missing ccds")

    print()
    
    return exp_shdu_missing_all


def get_exp(patches):
    """Get Exp.

    Return set of exposures from a list of patches.

    Parameters
    ----------
    patches: list
        input patches

    Returns
    -------
    set
       exposure IDs

    """
    exp_all = set()
    
    for patch in patches:

        path_exp = f"{patch}/exp_numbers.txt"
        exp = get_lines(path_exp)

        print(f"Patch {patch}: Found {len(exp)} exposures", end="; ")

        exp_all.update(exp)

        print(f"cumulative {len(exp_all)} exposures")

    print()

    return exp_all


def get_ccds_with_psf(patches, n_CCD=40):
    """Get CCDs With PSF.

    Return set of CCDs from list of patches.

    Parameters
    ----------
    patches: list
        input patches

    Returns
    -------
    set
       CCD IDs

    """
    # Get missing CCDs
    print("=== get missing CCDs ===")
    exp_shdu_missing_all = get_exp_shdu_missing(patches)

    # Get all exposures used in tiles
    print("=== get exposures ===")
    exp_all = get_exp(patches)

    # Turn exposures into exposure-single-HDU names (CCDs)
    exp_shdu_all = summary.get_all_shdus(exp_all, n_CCD)

    print(f"Found {len(exp_shdu_all)} CCDs")

    return exp_shdu_all


def get_ccds_with_psf_method_2(patches, n_CCD=40):

    for patch in patches:
        directory = f"{patch}/exp_runs"

    for entry in os.scandir(directory):
        pass

def save(IDs, path):
    """Save.

    Save list of IDs to text file.

    Parameters
    ----------
    IDs: set
        input IDs

    path: str
        output file name

    """
    with open(path, "w") as f_out:
        for ID in IDs:
            print(ID, file=f_out)

def main(argv):
    """Main.

    Main program.

    """
    version = "v1.5"

    if version == "v1.4":
        n_patch = 7
    elif version == "v1.5":
        n_patch = 8
    else:
        raise ValueError(f"Invalid version {version}")
        
    patches = [f'P{x}' for x in np.arange(n_patch) + 1]

    print(f"=== get_ccds_with_psf for version {version}, patches {patches} ===")

    print("=== method 1: exp_list - missing === ")
    exp_shdu_all = get_ccds_with_psf(patches)

    save(exp_shdu_all, f"ccds_with_psf_{version}.txt")

    #print("=== method 2: star cats === ")
    #exp_shdu_all_method_2 = get_ccds_with_psf_method_2(patches)
    #save(exp_shdu_all_method_2, f"ccds_with_psf_{version}_method_2.txt")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
