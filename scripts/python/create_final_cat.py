#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Script create_final_cat.py

Create and update hdf5 file of all final ShapePipe output FITS files, runs of
ShapePipe module ``make_catalogue_runner``. Supercedes `merge_final_cat.py`.

Usage: in parent dir of patches:
create_final_cat.py -p ~/shapepipe/example/cfis/final_cat.param -i . -P 7 -v -m final_cat_P7.hdf5
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
        "merged_cat_path": "final_cat.hdf5",
        "param_path": None,
        "hdu_num": 1,
        "patch": "3",
        "list_only": False,
    }
    _short_options = {                                                 
        "input_root_dir": "-i",                                           
        "merged_cat_path": "-m",
        "param_path": "-p",
        "patch": "-P",
        "list_only": "-l",

    }                                                                       
    _types = {                                                         
        "hdu_num": "int",                                                   
        "list_only": "bool",
    }                                                                       
    _help_strings = {                                                  
        "input_root_dir": "input root_dir for tile catalogues, default={}",
        "merged_cat_path": "merged catalogue path (hdf5 file), default={}",
        "param_path": "parameter file path, if not given use all columns, default={}",
        "patch": "patch number, default={}",
        "list_only": "print list of patches and IDs only, default={}",
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
    
    return True


def print_list(params):

    with h5py.File(params["merged_cat_path"], "a") as hdf5_file:

        for patch in hdf5_file["patches"]:
            print(patch)
            for id in hdf5_file[f"patches/{patch}"]:
                print(f" {id}")


def process(params):
    
    patch = rf"P{params['patch']}"
    patch_pattern = re.compile(patch)
    
    # Define the nested path pattern for locating `.fits` files
    file_pattern = "output/run_sp_Mc_*/make_cat_runner/output/*.fits"

    # Regex pattern for tile IDs
    id_pattern = re.compile(r"^\d+\.\d+$")
    
    # Open the HDF5 file (create it if it doesn't exist)
    if params["verbose"]:
        print(f"Initializing file {params['merged_cat_path']}")
    with h5py.File(params["merged_cat_path"], "a") as hdf5_file:

        # Iterate over patches
        for patch in os.listdir(params["input_root_dir"]):
            if not patch_pattern.fullmatch(patch):  # Skip non-matching entries
                #if params["verbose"]:
                    #print(f"Patch {patch} not found in {params['input_root_dir']}, skipping")
                continue

            patch_path = os.path.join(params["input_root_dir"], patch)

            if not os.path.isdir(patch_path):
                if params["verbose"]:
                    print(f"Path {patch_path} not found, skipping")
                continue
        
            # Get or create a group for the patch
            if f"patches/{patch}" in hdf5_file:
                if params["verbose"]:
                    print(f"Processing group {patch}")
                patch_group = hdf5_file[f"patches/{patch}"]
            else:
                if params["verbose"]:
                    print(f"Creating new group for patch {patch}")
                patch_group = hdf5_file.create_group(f"patches/{patch}")
        
            tile_runs_path = os.path.join(patch_path, "tile_runs")
            subdirs = os.listdir(tile_runs_path)
                
            for id in tqdm.tqdm(subdirs, total=len(subdirs)):

                # Skip if the dataset already exists
                if id in patch_group:
                    if params["verbose"]:
                        print(f"Skipping {id} (already processed)")
                        continue

                if id_pattern.match(id) and os.path.isdir(os.path.join(tile_runs_path, id)):
                    id_path = os.path.join(tile_runs_path, id)
                
                    output_path = os.path.join(id_path, "output", "run_sp_Mc_*", "make_cat_runner", "output", "*.fits")    
                    fits_files = glob.glob(output_path, recursive=True)
                    for fits_file in fits_files:
                    
                        # Read data
                        with fits.open(fits_file) as hdul:
                            data = hdul[params["hdu_num"]].data

                            # If columns not given on input: read all column names
                            if params["param_list"] is None:
                                params["param_list"] = [col for col in data.keys()]

                            try:
                                extracted_data = {col: data[col] for col in params["param_list"]}
                                dtype = data.dtype
                            except:
                                print("Error for ID {id}")
                                for col in params["param_list"]:
                                    if col not in data:
                                        print(col, end=" ")
                                    print()
                                continue

                        # For columns with 2-entry quantities, ensure proper dtype
                        #dtype = []
                        structured_data = np.empty(
                            len(extracted_data[params["param_list"][0]]),
                            dtype=dtype,
                        )

                        # Loop over parameters
                        for col in params["param_list"]:

                            #if isinstance(extracted_data[col][0], (np.ndarray, tuple, list)):
                            if False:

                                continue
                                # If multi-entry loop over entries
                                print(col)
                                num_elements = len(extracted_data[col][0])
                                for idx in range(num_elements):
                                    new_col_name = f"{col}_{idx}"
                                    dtype.append((new_col_name, extracted_data[col].dtype))
                                    try:
                                        structured_data[new_col_name] = [x[idx] for x in extracted_data[col]]  # Extract element `idx`
                                    except:
                                        print(f"Error for ID {id} for column {new_col_name}")
                                        raise
                            else:
                                # single entry
                                if not col in extracted_data:
                                    print(f"Column {col} not in file with ID {id}")
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

    if params["list_only"]:
        print_list(params)
    else:
        params["param_list"] = read_param_file(
            params["param_path"],
            verbose=params["verbose"],
        )
        process(params)

    return 0


if __name__ == "__main__":                                                      
    sys.exit(main(sys.argv))



