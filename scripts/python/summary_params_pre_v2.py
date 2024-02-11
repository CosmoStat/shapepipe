# Parameters for summary run

import os
from shapepipe.utilities.summary import *

def init_par_runtime(list_tile_IDs):
    
    # Numbers updated at runtime 
    par_runtime = {}

    par_runtime["n_tile_IDs"] = len(list_tile_IDs)
    par_runtime["list_tile_IDs"] = list_tile_IDs

    return par_runtime


def update_par_runtime_after_find_exp(par_runtime, all_exposures):
    
    n_CCD = 40
    
    # Single-HDU single exposure images
    par_runtime["n_shdus"] = get_par_runtime(par_runtime, "exposures") * n_CCD
    par_runtime["list_shdus"] = get_all_shdus(all_exposures, n_CCD)

    ## For split_exposure_runner, the output is image, weight,flag per single-HDU image
    ## and a header per exposure.
    par_runtime["n_3*n_shdus+n_exposures"] = (
        3 * get_par_runtime(par_runtime, "shdus")
        + get_par_runtime(par_runtime, "exposures")
    )
    
    return par_runtime


def set_jobs_v2_pre_v2(patch, verbose):
    """ Return information about shapepipe jobs
    
    """
    print(f"Set job info for patch {patch}")

    # Main input and output directory
    path_main = f"{os.environ['HOME']}/cosmostat/v2/pre_v2/psfex/{patch}"

    # Logging
    path = f"{path_main}/summary"
    if not os.path.isdir(path):
        os.mkdir(path)
    log_file_name = f"{path}/summary_log.txt"
    handlers = [
        logging.FileHandler(log_file_name, mode="w"),
        logging.StreamHandler()
    ]
    logging.basicConfig(
        level=logging.INFO, format="%(message)s", handlers=handlers
    )

    logging.info(f"Checking main directory = {path_main}")

    # Tile IDs
    tile_ID_path = f"{path_main}/tile_numbers.txt"

    ## Tile IDs with dots
    list_tile_IDs_dot = get_IDs_from_file(tile_ID_path)

    jobs = {}

    # Set the first job (retrieve images)
    
    # With "CFIS_" only the linked images are counted. The original
    # ones do not match the IDdash pattern.
    # If images were downloaded in several runs:
    # - Only copy original images, then (re-)set links in SP numbering format
    # - get_images_runner_run_[12] consistent
    # - remove previous output dirs since only last is searched
    jobs["1"] = job_data(
        1,
         "run_sp_Git",
       [
            "get_images_runner_run_1",
            "find_exposures_runner",
            "get_images_runner_run_2",
        ],
        ["tile_IDs", "tile_IDs", "exposures"],
        pattern=["CFIS_", "", ""],
        n_mult=[2, 1, 3],
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    #        n_mult=[1, 1, 1],
    jobs["2"] = job_data(
        2,
        ["run_sp_Uz", "run_sp_exp_SpMh", "run_sp_exp_SpMh"],
        ["uncompress_fits_runner", "merge_headers_runner", "split_exp_runner"],
        ["tile_IDs", 0, "3*n_shdus+n_exposures"],
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    run_dir_mask_tiles = "run_sp_tile_Ma"
    run_dir_mask_exp = "run_sp_exp_Ma"
    mask_module_tiles = "mask_runner"
    mask_module_exp = "mask_runner"

    jobs["4"] = job_data(
        4,
        run_dir_mask_tiles,
        [mask_module_tiles],
        ["tile_IDs"],
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    jobs["8"] = job_data(
        8,
        run_dir_mask_exp,
        [mask_module_exp],
        ["shdus"],
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    jobs["16"] = job_data(
        16,
        "run_sp_tile_Sx",
        ["sextractor_runner"],
        ["tile_IDs"],
        n_mult=2,
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    # TODO 1 setools_runner output/rand_split
    # TODO 2 add back Pi
    jobs["32"] = job_data(
        32,
        [
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            #"run_sp_exp_Pi"
        ],
        [
            "sextractor_runner",
            "setools_runner",
            "psfex_runner",
            # "psfex_interp_runner"],
        ],
        "shdus",
        n_mult=[2, 2, 2],  # 1],
        path_main=path_main,
        path_left="exp_runs",
        output_subdirs="shdus",
        path_right="output",
        verbose=verbose,
    )

    # For P3
    #jobs["33"] = job_data(
    #    33,
    #    "run_sp_exp_Pi",
    #    ["psfex_interp_runner"],
    #    "shdus",
    #    path_main=path_main,
    #    path_left="exp_runs",
    #    output_subdirs="shdus",
    #    path_right="output",
    #    verbose=verbose,
    #)

    jobs["64"] = job_data(
        "64",
        "run_sp_tile_PsViSmVi",
        [
            "psfex_interp_runner",
            "vignetmaker_runner_run_1",
            "spread_model_runner",
            "vignetmaker_runner_run_2",
        ],
        "tile_IDs",
        n_mult=[1, 1, 1, 4],
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    n_sh = 8
    run_dirs = [f"run_sp_tile_ngmix_Ng{idx+1}u" for idx in range(n_sh)]
    output_path_missing_IDs = [
        f"{path_main}/summary/missing_job_128_ngmix_runner_{idx+1}.txt" for idx in range(n_sh)
    ]
    jobs["128"] = job_data(
        "128",
        run_dirs,
        ["ngmix_runner"] * 8,
        "tile_IDs",
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        output_path_missing_IDs=output_path_missing_IDs,
        verbose=verbose,
    )

    jobs["256"] = job_data(
        "256",
        ["run_sp_Ms", "run_sp_Mc"],
        ["merge_sep_cats_runner", "make_cat_runner"],
        "tile_IDs",
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    # Post-processing
    jobs["1024"] = job_data(
        "1024",
        "run_sp_combined_psf",
        ["psfex_interp_runner"],
        "shdus",
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )
    
    return jobs, list_tile_IDs_dot
