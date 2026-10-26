#!/usr/bin/env python

"""Script create_final_cat.py

Create or append per-tile final catalogues to an HDF5 file, one structured
dataset per tile under ``patches/<patch>``. Inputs are ``make_cat_runner``
outputs in HDF5 (one dataset per column) or FITS format.

Existing tile datasets are skipped, even if their sources change. For a
workflow output reconciled against campaign membership and source changes,
use ``workflow/scripts/merge_final_cat.py``.

Usage: in parent dir of patches:
create_final_cat.py -p ~/shapepipe/workflow/config/cfis/final_cat.param -i . -P 7 -v -m final_cat_P7.hdf5

:Author: Martin Kilbinger

:Date: 2025

"""

import sys
import os
import re
import numpy as np
import tqdm
import glob
import h5py

from cs_util import args as cs_args
from cs_util import logging

from shapepipe.utilities.final_cat import read_final_cat


def params_from_run_config(params, defaults):
    """Replace default-valued paths from an image-simulation workflow run config.

    Use ``workflow/scripts/run_config.py`` for config resolution and precedence.
    Derive image-simulation paths only; data patch naming is not handled here.

    Values equal to their defaults are replaced, including explicitly supplied
    flags with default values. Nondefault flag values win.
    """
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, os.path.join(repo, "workflow", "scripts"))
    import run_config as _rc

    cfg = _rc.load(os.path.join(repo, "workflow", "config.yaml"),
                   params["run_config"])
    if cfg.get("input_type") != "image_sims":
        raise ValueError(
            f"run config {params['run_config']} has input_type="
            f"{cfg.get('input_type', 'data')!r}; -c derives paths for "
            "image_sims only")

    outputs = cfg.get("outputs") or {}
    products = outputs.get("products_dir") or outputs.get("run_dir")
    if not products:
        raise ValueError(f"run config {params['run_config']} sets neither "
                         "outputs.products_dir nor outputs.run_dir")

    # The group and file are named for `run:`, as the workflow's
    # final_cat_merge names them. -P is also the directory the -I walk
    # matches under -i, so this needs the run template's layout,
    # products_dir = <root>/<run>/product, and -i is <root>.
    patch = cfg["run"]
    patch_dir = os.path.dirname(os.path.normpath(products))
    if os.path.basename(patch_dir) != patch:
        raise ValueError(
            f"run config {params['run_config']}: products_dir {products} is "
            f"not <root>/{patch}/product, so -c cannot locate the tiles; "
            "pass -i and -P explicitly")
    derived = {
        "image_sims": True,
        "input_root_dir": os.path.dirname(patch_dir),
        "patch": patch,
        "merged_cat_path": os.path.join(products, f"final_cat_{patch}.hdf5"),
        "param_path": os.path.join(repo, "workflow", "config",
                                   "cfis_image_sims", "final_cat.param"),
    }
    for key, value in derived.items():
        if params.get(key) == defaults.get(key):
            params[key] = value
    # The tile count is the file's n_tiles attribute; the text summary is
    # written only when -o asks for it.
    if params.get("output_summary") == defaults.get("output_summary"):
        params["output_summary"] = None
    return params


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

    # See params_from_run_config for default-value replacement semantics.
    if _params.get("run_config"):
        _defaults, _, _, _ = params_default()
        _params = params_from_run_config(_params, _defaults)
                                                                                
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
        "output_summary": "n_tiles_final.txt",
        "ID": None,
        "single_op": None,
        "image_sims": False,
        "run_config": None,
    }
    _short_options = {
        "input_root_dir": "-i",
        "merged_cat_path": "-m",
        "param_path": "-p",
        "patch": "-P",
        "list_only": "-l",
        "output_summary": "-o",
        "single_op": "-s",
        "image_sims": "-I",
        "run_config": "-c",
    }
    _types = {
        "hdu_num": "int",
        "list_only": "bool",
        "image_sims": "bool",
    }
    _help_strings = {
        "input_root_dir": "input root_dir for tile catalogues, default={}",
        "merged_cat_path": "merged catalogue path (hdf5 file), default={}",
        "param_path": "parameter file path, if not given use all columns, default={}",
        "patch": "patch number (data) or grid subdir (image_sims), default={}",
        "list_only": "print list of patches and IDs only, default={}",
        "output_summary": "output file for number of tiles (with -c, written"
                          " only if given), default={}",
        "ID": "ID for single-ID operation, default={}",
        "single_op": "single ID operation, allowed are 'check', 'add', 'remove'; default={}",
        "image_sims": "image simulations mode (different dir layout and run prefix), default={}",
        "run_config": "workflow run config (e.g. sp_1p2z_grid_1.yaml); fills"
                      " in the paths below that were not given explicitly,"
                      " default={}",
    }

    return _params, _short_options, _types, _help_strings


def read_param_file(path, verbose=False):
    """Return unique parameter names in file order.

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
            print(f"Read {len(param_list)} columns", end="")
        else:
            print("No parameters read", end="")
        print(" into merged catalogue")

    # Preserve parameter-file order while removing duplicates; column order
    # is part of the output's structured dtype.
    param_list_unique = list(dict.fromkeys(param_list))

    if verbose:
        n = len(param_list) - len(param_list_unique)
        if n > 0:
            print(f"Removed {n} duplicate entries")

    return param_list_unique


def check_params(params):
 
    if bool(params["ID"] is None) != bool(params["single_op"] is None):
        print("Both or none of the options 'ID' and 'single_op' need to be specified")
        return False

    if bool(params["single_op"] is not None):
        allowed = ("check", "add", "remove")
        if params["single_op"] not in allowed:
            print(
                f"Invalid single_op '{params['single_op']}', allowed are"
                + f" {', '.join(allowed)}"
            )
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

    verbose = params.get("verbose", False)
    n_tiles = 0

    if not os.path.exists(params["merged_cat_path"]):
        print(f"File {params['merged_cat_path']} not found")
        return

    with h5py.File(params["merged_cat_path"], "r") as hdf5_file:

        if "patches" not in hdf5_file:
            print("Warning: no 'patches' group in output file (0 tiles added?)")
        else:
            for patch in hdf5_file["patches"]:
                for id in hdf5_file[f"patches/{patch}"]:
                    n_tiles += 1
                    if verbose:
                        print(f"  {patch}/{id}")

    if verbose:
        print(f"Total: {n_tiles} tiles")

    if params["output_summary"]:
        with open(params["output_summary"], "w") as f_out:
            print(n_tiles, file=f_out)

    # Write n_tiles to HDF5 file header
    with h5py.File(params["merged_cat_path"], "a") as hdf5_file:
        hdf5_file.attrs["n_tiles"] = n_tiles


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


def read_data(cat_file, params):
    """Read the parameter list's columns out of one catalogue.

    The catalogue is HDF5 or FITS, decided per file
    (``shapepipe.utilities.final_cat.read_final_cat``); a FITS table is read
    at HDU ``params["hdu_num"]`` (default 1). Returns the requested columns
    and a structured dtype describing them.

    @sc [label:schema] read-data-raises-on-missing-column
    A requested column the catalogue lacks raises `KeyError` naming it; it is
    never skipped or filled. `copy_data` keeps only columns present in the
    source, so this raise is the one place a missing name stops a merge, and
    without it a tile short a per-epoch slot would land in the merged file
    silently narrower, with that slot's exposure identity gone. Enforced by
    tests/unit/test_final_cat_merge_invariants.py.
    """
    data = read_final_cat(
        cat_file, params["param_list"], hdu=params.get("hdu_num", 1)
    )

    # If columns not given on input: read all column names
    if params["param_list"] is None:
        params["param_list"] = list(data.dtype.names)

    extracted_data = {col: data[col] for col in params["param_list"]}

    return extracted_data, data.dtype


def copy_data(param_list, extracted_data, dtype):
    """Copy requested columns into a structured array with their stored dtypes.

    @sc [label:schema] final-copy-requested-column-order
    Allocate only requested columns present in the source, in parameter-file
    order. Unrequested fields must not enter the output as uninitialised memory,
    and source column order must not make tile dtypes incompatible at
    concatenation. ``read_data`` enforces missing-column errors before this call.
    """
    wanted = set(dtype.names or ())
    columns = [col for col in param_list if col in wanted]
    subset = np.dtype([(col, dtype[col]) for col in columns])

    # Initialize new data structure
    structured_data = np.empty(
        len(extracted_data[param_list[0]]),
        dtype=subset,
    )

    # Loop over parameters
    for col in columns:
        structured_data[col] = extracted_data[col]
    
    #if isinstance(extracted_data[col][0], (np.ndarray, tuple, list)):
    if False:

        # If multi-entry loop over entries
        print(col)
        num_elements = len(extracted_data[col][0])
        for idx in range(num_elements):
            new_col_name = f"{col}_{idx}"
            dtype.append((new_col_name, extracted_data[col].dtype))
            try:
                structured_data[new_col_name] = [x[idx] for x in extracted_data[col]]
            except:
                print(f"Error for ID {id} for column {new_col_name}")
                raise

    return structured_data


def collect_tile_ids_image_sims(patch_path):
    """Collect tile IDs from image-sims layout: tiles/<prefix>/<tile_id>/

    Parameters
    ----------
    patch_path : str
        path to the grid subdir (e.g. .../1p2z_grid_1)

    Returns
    -------
    list of (tile_id, tile_path) tuples
    """
    id_pattern = re.compile(r"^\d+\.\d+$")
    # Use the first existing tiles root: direct, product/, then products/.
    # See workflow/scripts/clean_tile.py for persistent catalogue ownership.
    result = []
    for sub in ("tiles", os.path.join("product", "tiles"),
                os.path.join("products", "tiles")):
        tiles_root = os.path.join(patch_path, sub)
        if os.path.isdir(tiles_root):
            break
    else:
        return result
    for prefix in os.listdir(tiles_root):
        prefix_path = os.path.join(tiles_root, prefix)
        if not os.path.isdir(prefix_path):
            continue
        for tile_id in os.listdir(prefix_path):
            if id_pattern.match(tile_id):
                result.append((tile_id, os.path.join(prefix_path, tile_id)))
    return result


# Prefer HDF5 when both supported formats exist.
FINAL_CAT_SUFFIXES = (".hdf5", ".fits")


def find_final_cat(id, id_path, run_prefix):
    """Path of this tile's final catalogue, or None.

    Flat layout first: the unified workflow copies one catalogue per tile to
    <products>/tiles/<shard>/<tile>/final_cat-<tile>.hdf5, in DOT form and
    with no run sub-tree. Then the run-dir layout, where the catalogue sits
    under the tile's own shapepipe run dir in DASH form, newest run wins.
    Either layout may hold an hdf5 or a FITS catalogue; hdf5 is preferred.
    """
    for suffix in FINAL_CAT_SUFFIXES:
        flat = os.path.join(id_path, f"final_cat-{id}{suffix}")
        if os.path.exists(flat):
            return flat

    base_pattern = os.path.join(id_path, "output", run_prefix)
    all_matches = [d for d in glob.glob(base_pattern) if os.path.isdir(d)]
    if not all_matches:
        return None
    newest_dir = max(all_matches, key=os.path.getmtime)
    id_dash = re.sub(r"\.", "-", id)
    for suffix in FINAL_CAT_SUFFIXES:
        in_run = (f"{newest_dir}/make_cat_runner/output/"
                  f"final_cat-{id_dash}{suffix}")
        if os.path.exists(in_run):
            return in_run
    return None


def process(params):

    if params["image_sims"]:
        patch_name = params["patch"]
        run_prefix = "run_sp_tile_Mc*"
    else:
        patch_name = rf"P{params['patch']}"
        run_prefix = "run_sp_tile_Mc_*"

    patch_pattern = re.compile(patch_name)

    # Regex pattern for tile IDs
    id_pattern = re.compile(r"^\d+\.\d+$")

    n_added = 0
    IDs_added = []

    # Open the HDF5 file (create it if it doesn't exist)
    if params["verbose"]:
        print(f"Initializing file {params['merged_cat_path']}")
    with h5py.File(params["merged_cat_path"], "a") as hdf5_file:

        # Iterate over patches
        for patch in os.listdir(params["input_root_dir"]):

            # Skip non-matching entries
            if not patch_pattern.fullmatch(patch):
                continue

            # Full path to patch
            patch_path = os.path.join(params["input_root_dir"], patch)
            if not os.path.isdir(patch_path):
                if params["verbose"]:
                    print(f"Path {patch_path} not found, skipping")
                continue

            # Get hdf5 group for this patch
            patch_group = get_patch_group(hdf5_file, patch, params["verbose"])

            # Collect (tile_id, tile_path) pairs depending on layout
            if params["image_sims"]:
                tile_items = collect_tile_ids_image_sims(patch_path)
            else:
                tile_runs_path = os.path.join(patch_path, "tile_runs")
                tile_items = [
                    (tid, os.path.join(tile_runs_path, tid))
                    for tid in os.listdir(tile_runs_path)
                    if id_pattern.match(tid)
                    and os.path.isdir(os.path.join(tile_runs_path, tid))
                ]

            for id, id_path in tqdm.tqdm(tile_items, total=len(tile_items)):

                # Skip if the patch/ID data already exists
                if id in patch_group:
                    if params["verbose"]:
                        print(f"Skipping {id} (already processed)")
                    continue

                cat_file = find_final_cat(id, id_path, run_prefix)
                if cat_file is None:
                    if params["verbose"]:
                        print(f"Final cat for {id} not found, continuing")
                    continue

                extracted_data, dtype = read_data(cat_file, params)

                structured_data = copy_data(params["param_list"], extracted_data, dtype)

                # Use copy_data's selected-column dtype, not the source dtype.
                try:
                    patch_group.create_dataset(
                        str(id),
                        data=structured_data,
                        dtype=structured_data.dtype,
                    )
                except:
                    print(f"Error for {id}: Could not create dataset in group {patch}")
                    raise

                n_added += 1
                IDs_added.append(id)

        if params["verbose"]:
            print(f"{n_added} tiles added ({' '.join(IDs_added)})")


def single_action(params):
    """Single Action.

    Check or remove one ID. ``add`` returns ``None`` and lets ``main`` run the
    ordinary catalogue walk; it does not restrict that walk to the requested ID.

    Parameters
    ----------
    params : dict
        parameter options)

    Returns
    -------
    int
        0 after check or remove; ``None`` when no single-ID action is performed.
        The check branch returns 0 even if the ID is absent.

    """
    if params["single_op"] == "check":
        # Check wheter ID is part of hdf5 file
        found = check_ID(params["merged_cat_path"], params["ID"], verbose=True)

        if found == "":
            # ID not found
            res = 1
        # This assignment also runs when the ID is absent.
        res = 0

    elif params["single_op"] == "remove":
        # Remove ID
        remove_ID(params["merged_cat_path"], params["ID"], verbose=True)
        res = 0

    else:

        # No action performed
        res = None

    return res


def main(argv=None):
    
    if argv is None:
        argv = sys.argv[0:]
    params = set_params_from_command_line(argv)

    if check_params(params) == False:
        return 1

    res = single_action(params)
    if res is not None:
        return res

    if params["list_only"] == False:
        params["param_list"] = read_param_file(
            params["param_path"],
            verbose=params["verbose"],
        )
        process(params)

    print_list(params)

    return 0


if __name__ == "__main__":                                                      
    sys.exit(main(sys.argv))
