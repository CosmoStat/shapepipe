#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script link_to_exp_for_tile.py

:Description: Link to exposure and PSF catalogue
 for a given tile.

:Author: Martin Kilbinger

"""

import os
import sys
import re
import copy

from optparse import OptionParser


class param:
    """General class to store (default) variables"""

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def print(self, **kwds):
        print(self.__dict__)

    def var_list(self, **kwds):
        return vars(self)


def params_default():
    """Params Default.

    Set default parameter values.

    Returns
    -------
    class param
        parameter values

    """
    p_def = param(
        tile_base_dir=".",
        exp_base_dir=".",
        sp_local=0,
    )

    return p_def


def parse_options(p_def):
    """Parse Options.

    Parse command line options.

    Parameters
    ----------
    p_def: class param
        parameter values

    Returns
    -------
    list
        command line options
        command line str

    """
    usage = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # IO
    parser.add_option(
        "-i",
        "--input_tile_dir",
        dest="tile_base_dir",
        type="string",
        default=p_def.tile_base_dir,
        help=f"input tile base directory, default='{p_def.tile_base_dir}'",
    )
    parser.add_option(
        "-t",
        "--tile_ID",
        dest="tile_ID",
        type="string",
        help=f"input tile ID",
    )
    parser.add_option(
        "-I",
        "--input_exp_dir",
        dest="exp_base_dir",
        type="string",
        default=p_def.exp_base_dir,
        help=f"input exposure base directory, default='{p_def.exp_base_dir}'",
    )
    parser.add_option(
        "-s",
        "--sp_local",
        dest="sp_local",
        type="int",
        default=p_def.sp_local,
        help=f"local runi of split_exposure_runner, default='{p_def.sp_local}'",
    )
    parser.add_option(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="verbose output",
    )

    options, args = parser.parse_args()

    return options, args


def check_options(options):
    """Check Options.

    Check command line options.

    Parameters
    ----------
    options: tuple
        Command line options

    Returns
    -------
    bool
        Result of option check. False if invalid option value.

    """
    return True


def update_param(p_def, options):
    """Update Param.

    Return default parameter, updated and complemented according to options.

    Parameters
    ----------
    p_def:  class param
        parameter values
    optiosn: tuple
        command line options

    Returns
    -------
    class param
        updated paramter values

    """
    param = copy.copy(p_def)

    # Update keys in param according to options values
    for key in vars(param):
        if key in vars(options):
            setattr(param, key, getattr(options, key))

    # Add remaining keys from options to param
    for key in vars(options):
        if not key in vars(param):
            setattr(param, key, getattr(options, key))

    # Do extra stuff if necessary

    return param


# TODO: move to cs_util
def matching_subdirs(base_dir, pattern):
    """Matching Subdirs.

    Return all subdirectories whose name start with a given string pattern.

    Parameters
    ----------
    base_dir: str
        base directory
    pattern: str
        pattern to match beginning of subdir name

    Returns
    -------
    list
        matched subdirectory names

    """
    subdirs = []
    if os.path.exists(base_dir):
        # Loop over directory entries
        for entry in os.listdir(base_dir):
            full_path = os.path.join(base_dir, entry)
            # Append if match
            if os.path.isdir(full_path) and entry.startswith(pattern):
                subdirs.append(full_path)
    else:
        print(f"Warning: {base_dir} does not exist, continuing...")

    # Sort according to creation date
    subdirs.sort(key=os.path.getctime)

    return subdirs


def get_tile_out_dir(tile_base_dir, tile_ID):
    """Get Tile Out Dir.

    Return output directory path for a given tile ID.

    Parameters
    ----------
    tile_base_dir: str
        base directory of tile runs
    tile_ID: str
        tile ID

    Returns
    -------
    str
        path

    """
    return f"{tile_base_dir}/{tile_ID}/output"


def get_exp_IDs(tile_base_dir, tile_ID, verbose=False):
    """Get Exp IDs.

    Return exposure IDs used for given tile.

    Parameters
    ----------
    tile_base_dir: str
        base directory of tile runs
    tile_ID: str
        tile ID
    verbose: bool, optional
        verbose output if ``True``; default is ``False``

    Returns
    -------
    list
        exposure IDs
        
    """
    # Get tile output path
    tile_out_dir = get_tile_out_dir(tile_base_dir, tile_ID)

    # Get subdirectories with runs "Fe" (find exposures)
    pattern = "run_sp_GitFeGie"
    subdirs = matching_subdirs(tile_out_dir, pattern)

    # Raise error if not exactly one run
    if len(subdirs) == 0:
        raise IOError(
            f"No matching directory '{pattern}' in {tile_out_dir} found"
        )
    if len(subdirs) != 1:
        raise IOError(
            f"Exactly one directory natching {pattern} in {tile_out_dir} "
            + f"expected, not {len(subdirs)}"
        )

    # Replace dot with dash in tile ID
    tile_ID_sp = re.sub(r"\.", "-", tile_ID)
    # Get output file of find_exposure_runner for this tile ID
    exp_ID_file = 
        f"{subdirs[0]}/find_exposures_runner/output/"
        + f"exp_numbers-{tile_ID_sp}.txt"
    )

    exp_IDs = []
    # Read file to get exposure IDs
    with open(exp_ID_file) as f_in:
        for line in f_in:
            name = line.strip()
            # Remove any letter
            ID = re.sub("[a-zA-Z]", "", name)
            exp_IDs.append(ID)

    if verbose:
        print("Exposures: ", exp_IDs)
    return exp_IDs


def get_exp_single_HDU_IDs(exp_IDs, n_CPU):
    """Get Exp Single HDU IDs.

    Return all single-HDU IDs for given exposures.

    Parameters
    ----------
    exp_IDs: list
        input exposure IDs, list of str
    n_CPU: int
        number of CPUs (=HDUs) of an exposure; n_CPU=40 for MegaCAM

    Returns
    -------
    list
        list of single-HDU exposure IDs

    """
    exp_shdu_IDs = []
    # Loop over input exposure IDs
    for exp_ID in exp_IDs:
        # Append int up to n_CPU to each
        for idx in range(n_CPU):
            ID = f"{exp_ID}-{idx}"
            exp_shdu_IDs.append(ID)

    return exp_shdu_IDs


def get_paths(exp_base_dir, exp_shdu_IDs, pattern):
    """Get Paths.

    Return (newest) subdirectory for each input single-exposure HDU ID that matches pattern.

    Parameters
    ----------
    exp_base_dir: str
        base directory for (single-HDU) exposure runs
    exp_shdu_IDs: list
        single-HDU exposure IDs; list of str
    pattern: str
        pattern to match beginning of subdir names
    Returns
    -------
    list
        matching paths; list of str, one entry for each input single-exp ID

    """
    paths = []
    # Loop over single-HDU exposure IDs
    for exp_shdu_ID in exp_shdu_IDs:

        # output path of runs
        name = f"{exp_base_dir}/{exp_shdu_ID}/output"
        path = os.path.abspath(name)

        # get matching subdirs
        subdirs = matching_subdirs(path, pattern)
        n_subdirs = len(subdirs)

        if n_subdirs != 1:
            msg = (
                f"Exactly one directory matching {pattern} in {path} expected,"
                + f"  not {n_subdirs}"
            )

            # If more than one found: sort by name = sort by date
            subdirs = sorted(subdirs)

            # No match
            if n_subdirs == 0:
                continue

        # Append matching subdir; if more than one append newest
        paths.append(f"{subdirs[-1]}")

    return paths


def create_links_paths(tile_base_dir, tile_ID, paths, verbose=False):
    """Create Links Paths.

    Create links to paths.

    Parameters
    ----------
    tile_base_dir: str
        base directory for tile runs
    tile_ID: str
        tile ID
    paths: list
        paths; list of str
    verbose: bool, optional
        verbose output if ``True``; default is ``False``
 
    """
    # Get tile output path
    tile_out_dir = get_tile_out_dir(tile_base_dir, tile_ID)

    # Loop over paths
    for path in paths:

        # Get destination = tile output dir + path tail (part of path after last slash)
        head, tail = os.path.split(path)
        src = path
        dst = f"{tile_out_dir}/{tail}"
        
        if os.path.exists(dst):
            
            src_existing = os.readlink(dst)
            if src_existing == src:
                # destination already points to source: skip
                if verbose:
                    print(
                        f"Warning: {src} <- {dst} already exists, no link created"
                    )
                continue
            else:
                # destination points to different source: create links with added index to distinguish
                # from existing one(s)
                idx = 1
                dst_orig = dst
                while True:
                    # Find destination name with lowest appended index that does not exist
                    dst = f"{dst_orig}_{idx}"
                    if os.path.exists(dst):
                        idx += 1
                    else:
                        # Found: create new link
                        if verbose:
                            print(f"link {src} <- {dst}")
                        os.symlink(src, dst)
                        break
        else:

            # destination does not exist: create link
            if verbose:
                print(f"link {src} <- {dst}")
            os.symlink(src, dst)


def main(argv=None):

    # Set default parameters
    p_def = params_default()

    # Command line options
    options, args = parse_options(p_def)

    if check_options(options) is False:
        return 1

    param = update_param(p_def, options)

    tile_base_dir = param.tile_base_dir
    exp_base_dir = param.exp_base_dir
    tile_ID = param.tile_ID
    n_CPU = 40
    verbose = param.verbose

    exp_IDs = get_exp_IDs(tile_base_dir, tile_ID, verbose=verbose)
    exp_shdu_IDs = get_exp_single_HDU_IDs(exp_IDs, n_CPU)

    # Note: psfex P3 is mostly run_sp_exp_SxSePsf
    patterns = ["run_sp_exp_SxSePsfPi"]  # , "run_sp_exp_Pi"]
    if param.sp_local == 1:
        patterns.append("run_sp_exp_Sp_shdu")
    for pattern in patterns:
        paths = get_paths(exp_base_dir, exp_shdu_IDs, pattern)

        create_links_paths(tile_base_dir, tile_ID, paths, verbose=verbose)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
