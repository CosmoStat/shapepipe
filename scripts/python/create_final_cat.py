#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script create_final_cat.py

Create and update hdf5 file of all final ShapePipe output FITS files, runs of
ShapePipe module ``make_catalogue_runner``. Supercedes `merge_final_cat.py`.

:Author: Martin Kilbinger

"""

import sys
import os
import re
import numpy as np
import tqdm
import glob
import h5py
from astropy.io import fits

from cs_util import args as cs_args
from cs_util import logging


def set_params_from_command_line(args):                               
    """Set Params From Command Line.                                        
                                                                                
    """                                                                     
    _params, _short_options, _types, _help_strings = (
        params_default()
    )
    
    # Read command line options                                             
    options = cs_args.parse_options(                                        
        _params,                                                       
        _short_options,                                                
        _types,                                                        
        _help_strings,                                                 
    )                                                                       

    # Update parameter values from options                                  
    for key in options:                                               
        _params[key] = options[key]                           
                                                                                
    del options                                                             
                                                                                
    # Save calling command                                                  
    logging.log_command(args)

    return _params


def params_default():                                                   
    """Params Default.                                                      
                                                                                
    Set default parameter values.                                           
                                                                                
    """                                                                     
    _params = {                                                        
        "input_root_dir": ".",
        "output_path": "final_cat.hdf5",
        "param_path": None,
        "hdu_num": 1,
    }
    _short_options = {                                                 
        "input_root_dir": "-i",                                           
        "output_path": "-o",                                                 
        "param_path": "-p",
    }                                                                       
    _types = {                                                         
        "hdu_num": "int",                                                   
    }                                                                       
    _help_strings = {                                                  
        "input_root_dir": "input root_dir for tile catalogues",     
        "output_path": "output path (hdf5 file)",
        "param_path": "parameter file path, if not given use all columns",
    }

    return _params, _short_options, _types, _help_strings


def read_param_file(path, verbose=False):
    """Read Param File.

    Return parameter list read from file.

    Parameters
    ----------
    path: str
        input file name
    verbose: bool, optional, default=False
        verbose output if True

    Returns
    -------
    list of str
        parameter names

    """
    param_list = []

    if path:
        if verbose:
            print(f"Reading parameter file {path}")

        with open(path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                entry = line.rstrip()
                if not entry or entry == "":
                    continue
                param_list.append(entry)

    if verbose:
        if len(param_list) > 0:
            print(f"Copying {len(param_list)} columns", end="")
        else:
            print("Copying all columns", end="")
        print(" into merged catalogue")

    param_list_unique = list(set(param_list))
    
    if verbose:
        n = len(param_list) - len(param_list_unique)
        if n > 1:
            print("Removed {n} duplicate entries")
    
    return param_list_unique


def check_params(params):
 
    if bool(params["ID"] is None) != bool(params["single_op"] is None):
        print("Both or none of the options 'ID' and 'single_op' need to be specified")
        return False

    return True


def remove_ID(merged_cat_path, ID, verbose=False):
    """Remove ID.

    Remove ID from merged catalogue file.

    Parameters
    ----------
    merged_cat_path : str
        catalogue hdf5 file path
    ID : str
        tile ID
    verbose : bool, optional
        verbose output if True

    """
    # First check ID
    patch_found = check_ID(merged_cat_path, ID, verbose=False)
    if not patch_found:
        raise KeyError(f"ID {ID} not found in file {merged_cat_path}")

    with h5py.File(merged_cat_path, "a") as hdf5_file:

        patch_group = get_patch_group(hdf5_file, patch_found, verbose=verbose)
        del patch_group[ID]
        if verbose:
            print(f"Removed {patch_found}/{ID}.")


def check_ID(merged_cat_path, ID, verbose=False):
    """Check ID.

    Check whether ID exists in merged catalogue file.

    Parameters
    ----------
    merged_cat_path : str
        catalogue hdf5 file path
    ID : str
        tile ID
    verbose : bool, optional
        verbose output if True

    Returns
    -------
    str
        patch of ID; "" if not found

    """
    with h5py.File(merged_cat_path, "r") as hdf5_file:
        for patch in hdf5_file["patches"]:
            for id in hdf5_file[f"patches/{patch}"]:
                if id == ID:
                    if verbose:
                        print(f"ID {ID} found under patch {patch}")
                    return patch
    if verbose:
        print(f"ID {ID} not found")
    return ""


def print_list(params):
    """Print List.

    Print number of tile IDs found in hdf5 file.

    Parameters
    ----------
    params: dict
        parameters

    """
    n_tiles = 0
    # Open hdf5 merge cat file
    with h5py.File(params["merged_cat_path"], "r") as hdf5_file:

        # Loop over patches contained in hdf5 file
        for patch in hdf5_file["patches"]:
            # Count tiles summed over all patches
            for id in hdf5_file[f"patches/{patch}"]:
                n_tiles += 1

    # Write number of tiles to output file
    with open(params["output_summary"], "w") as f_out:
        print(n_tiles, file=f_out)


def get_patch_group(hdf5_file, patch, verbose=False):
    """Get Patch group.

    Return group from hdf5 file for a given patch.

    Parameters
    ----------
    hdf5_file : class h5py.File
        input hdf5 file
    patch : str
        patch name
    verbose : bool, optional
        verbose output if ``True``; default is ``False``

    Returns
    -------
    h5py.Group
        group for given patch

    """
    if f"patches/{patch}" in hdf5_file:
        if verbose:
            print(f"Processing group {patch}")
        patch_group = hdf5_file[f"patches/{patch}"]
    else:
        # Create new group
        if verbose:
            print(f"Creating new group for patch {patch}")
        patch_group = hdf5_file.create_group(f"patches/{patch}")

    return patch_group


def read_data(fits_file, params):
    """Read Data.

    Read data from final merge catalogue FITS file.

    Parameters
    ----------
    fits_file: str
        intput FITS file path
    params: dict
        parameters

    Returns
    -------
    dict
        data
    list
        dtypes of data columns

    """
    # Open FITS file
    with fits.open(fits_file) as hdu_list:
        # Get data of indicated HDU number
        try:
            data = hdu_list[params["hdu_num"]].data
        except:
            print(f"Error with ID {id}, file{fits_file}")
            raise

    # If columns not given on input: read all column names
    if params["param_list"] is None:
        params["param_list"] = [col for col in data.keys()]

    # Copy data to output dict; copy dtypes
    try:
        extracted_data = {col: data[col] for col in params["param_list"]}
        dtype = data.dtype
    except:
        print(f"Error for ID {id}, path {fits_file}")
        for col in params["param_list"]:
            if col not in data:
                print(col, end=" ")
            print()
            continue

    return extracted_data, dtype


def copy_data(param_list, extracted_data, dtype):
    """Copy Data.

    """
    # Initialize new data structure
    structured_data = np.empty(
        len(extracted_data[param_list[0]]),
        dtype=dtype,
    )

    # Loop over parameters
    for col in param_list:
        if not col in extracted_data:
            print(f"Column {col} not in file with ID {id}")
        structured_data[col] = extracted_data[col]
    
    return True


def process(params):
    
    #patch_pattern = re.compile(r"P\d")  # Matches P followed by a single digit
    patch_pattern = re.compile(r"P3")  # Matches P followed by a single digit
    
    # Define the nested path pattern for locating `.fits` files
    file_pattern = "output/run_sp_Mc_*/make_cat_runner/output/*.fits"

    # Regex pattern for tile IDs
    id_pattern = re.compile(r"^\d+\.\d+$")
    
    # Open the HDF5 file (create it if it doesn't exist)
    if params["verbose"]:
        print(f"Initializing file {params['output_path']}")
    with h5py.File(params["output_path"], "a") as hdf5_file:

        # Iterate over patches
        for patch in os.listdir(params["input_root_dir"]):
            if not patch_pattern.fullmatch(patch):  # Skip non-matching entries
                continue

            patch_path = os.path.join(params["input_root_dir"], patch)

            if not os.path.isdir(patch_path):
                continue
        
            # Get or create a group for the patch
            if f"patches/{patch}" in hdf5_file:
                if params["verbose"]:
                    print(f"Extending group for patch {patch}")
                patch_group = hdf5_file[f"patches/{patch}"]
            else:
                if params["verbose"]:
                    print(f"Creating new group for patch {patch}")
                patch_group = hdf5_file.create_group(f"patches/{patch}")
        
            tile_runs_path = os.path.join(patch_path, "tile_runs")
            #if params["verbose"]:
                #print(f"Loop over tile IDs in {tile_runs_path}...")
            subdirs = os.listdir(tile_runs_path)
            #if params["verbose"]:
                #print(f"{len(subdirs)} IDs found")
                
            for id in tqdm.tqdm(subdirs, total=len(subdirs)):
                if id_pattern.match(id) and os.path.isdir(os.path.join(tile_runs_path, id)):
                    id_path = os.path.join(tile_runs_path, id)
                
                    output_path = os.path.join(id_path, "output", "run_sp_Mc_*", "make_cat_runner", "output", "*.fits")    
                    fits_files = glob.glob(output_path, recursive=True)
                    for fits_file in fits_files:
                    
                        # Skip if the dataset already exists
                        if id in patch_group:
                            if params["verbose"]:
                                print(f"Skipping {id} (already processed)")
                                continue
                    
                        #if params["verbose"]:
                            #print(f"Reading file {fits_file}...")
                        with fits.open(fits_file) as hdul:
                            data = hdul[params["hdu_num"]].data

                            # If columns not given on input: read all column names
                            if params["param_list"] is None:
                                params["param_list"] = [col for col in data.keys()]

                            try:
                                extracted_data = {col: data[col] for col in params["param_list"]}
                            except:
                                print("Error for ID {id}. Keys not found:")
                                for col in params["param_list"]:
                                    if col not in data:
                                        print(col, end=" ")
                                    print()
                                continue

                        # For columns with 2-entry quantities, ensure proper dtype
                        dtype = []
                        structured_data = np.empty(len(extracted_data[params["param_list"][0]]), dtype=dtype)
                        for col in params["param_list"]:
                            if isinstance(extracted_data[col][0], (np.ndarray, tuple, list)):
                                num_elements = len(extracted_data[col][0])
                                for idx in range(num_elements):
                                    new_col_name = f"{col}_{idx}"
                                    dtype.append((new_col_name, extracted_data[col].dtype))
                                    structured_data[new_col_name] = [x[i] for x in extracted_data[col]]  # Extract element `i`
                            else:
                                dtype.append((col, extracted_data[col].dtype))
                                structured_data[col] = extracted_data[col]
                        
                        # Create a new dataset
                        patch_group.create_dataset(
                            str(id),
                            data=structured_data,
                            dtype=dtype,
                        )


def main(argv=None):
    
    if argv is None:
        argv = sys.argv[0:]
    params = set_params_from_command_line(argv)

    check_params(params)

    params["param_list"] = read_param_file(params["param_path"], verbose=params["verbose"])

    process(params)

    return 0


if __name__ == "__main__":                                                      
    sys.exit(main(sys.argv))



