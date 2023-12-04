#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script link_to_esp_for_tile.py

:Author: Martin Kilbinger

"""

import os
import sys
import re
import copy

from optparse import OptionParser


class param:
    """General class to store (default) variables

    """
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
        tile_base_dir  = '.',
        exp_base_dir = '.',
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
    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    # IO
    parser.add_option(
        '-i',
        '--input_tile_dir',
        dest='tile_base_dir',
        type='string',
        default=p_def.tile_base_dir,
        help=f'input tile base directory, default=\'{p_def.tile_base_dir}\''
    )
    parser.add_option(
        '-t',
        '--tile_ID',
        dest='tile_ID',
        type='string',
        help=f"input tile ID",
    )
    parser.add_option(
        '-I',
        '--input_exp_dir',
        dest='exp_base_dir',
        type='string',
        default=p_def.exp_base_dir,
        help=f'input exposure base directory, default=\'{p_def.exp_base_dir}\''
    )
    parser.add_option(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        help='verbose output'
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

    # Find all matching subdirectories
    subdirs = []
    for entry in os.listdir(base_dir):
        full_path = os.path.join(base_dir, entry)
        if os.path.isdir(full_path) and entry.startswith(pattern):
            subdirs.append(full_path)

    # Sort according to creation date
    subdirs.sort(key=os.path.getctime)

    return subdirs


def get_tile_out_dir(tile_base_dir, tile_ID):

    tile_out_dir = f"{tile_base_dir}/{tile_ID}/output"

    return tile_out_dir


def get_exp_IDs(tile_base_dir, tile_ID, verbose=False):

    tile_out_dir = get_tile_out_dir(tile_base_dir, tile_ID)

    subdirs = matching_subdirs(tile_out_dir, "run_sp_GitFeGie")

    if len(subdirs) != 1:
        raise IOError(f"Exactly one matching directory in {tile_out_dir} expected, not {len(subdirs)}")

    # Replace dot with dash in tile ID
    tile_ID_sp = re.sub(r"\.", "-", tile_ID)
    exp_ID_file = f"{subdirs[0]}/find_exposures_runner/output/exp_numbers-{tile_ID_sp}.txt"

    exp_IDs = []
    with open(exp_ID_file) as f_in:
        for line in f_in:
            name = line.strip()
            # Remove any letter
            ID = re.sub("[a-zA-Z]", "", name) 
            exp_IDs.append(ID)

    return exp_IDs


def get_exp_single_HDU_IDs(exp_IDs, n_CPU):

    exp_shdu_IDs = []
    for exp_ID in exp_IDs:
        for idx in range(n_CPU):
            ID = f"{exp_ID}-{idx}"
            exp_shdu_IDs.append(ID)

    return exp_shdu_IDs

    
def get_paths(exp_base_dir, exp_shdu_IDs):

    number = {}
    paths = []
    for exp_shdu_ID in exp_shdu_IDs:
        name = f"{exp_base_dir}/{exp_shdu_ID}/output"
        path = os.path.abspath(name)
        subdirs = matching_subdirs(path, "run_sp_exp_SxSePsf")
        n_subdirs = len(subdirs)

        if n_subdirs not in number:
            number[n_subdirs] = 1
        else:
            number[n_subdirs] += 1

        if n_subdirs != 1:
            msg = f"Exactly one matching directory in {path} expected, not {n_subdirs}"
            print(msg)
            if n_subdirs == 0:
                continue
        paths.append(f"{subdirs[0]}")

    return paths, number


def create_links_paths(tile_base_dir, tile_ID, paths):

    tile_out_dir = get_tile_out_dir(tile_base_dir, tile_ID)

    for path in paths:

        head, tail = os.path.split(path)
        src = path
        dst = f"{tile_out_dir}/{tail}"
        if os.path.exists(dst):
            print("Warning: {dst} already exists, no link created")
        else:
            print(f"ln -s {src} {dst}")
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
    paths, number = get_paths(exp_base_dir, exp_shdu_IDs)
    print(number)

    create_links_paths(tile_base_dir, tile_ID, paths)


    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
