# Parameters for summary run

import os
from shapepipe.utilities.summary import *


def set_jobs_v2_pre_v2(patch, verbose):
    """Return information about shapepipe jobs"""
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
        logging.StreamHandler(),
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
        "run_sp_GitFeGie",
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

    if patch == "P8":
        jobs["2"] = job_data(
            2,
            ["run_sp_Uz", "run_sp_exp_Sp_shdu"],
            ["uncompress_fits_runner", "split_exp_runner"],
            ["tile_IDs", "shdus"],
            n_mult=[1, 3],
            path_main=path_main,
            path_left=["output", "exp_runs"],
            output_subdirs=[None, "shdus"],
            path_right=[None, "output"],
            verbose=verbose,
        )
        jobs["4096"] = job_data(
            4096,
            ["run_sp_exp_Sp_shdu"],
            ["split_exp_runner"],
            ["shdus"],
            n_mult=4,
            path_main=path_main,
            path_left="exp_runs",
            output_subdirs="shdus",
            path_right="output",
            verbose=verbose,
        )
    else:
        jobs["2"] = job_data(
            2,
            ["run_sp_Uz", "run_sp_exp_SpMh"],
            ["uncompress_fits_runner",  "split_exp_runner"],
            ["tile_IDs", "shdus"],
            n_mult=[1,  121],
            path_main=path_main,
            path_left="output",
            verbose=verbose,
        )

    jobs["4"] = job_data(
        4,
        ["run_sp_Ma_tile"],
        ["mask_runner"],
        ["tile_IDs"],
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    if patch != "P8":
        jobs["8"] = job_data(
            8,
            ["run_sp_Ma_exp"],
            ["mask_runner"],
            ["shdus"],
            path_main=path_main,
            path_left="output",
            verbose=verbose,
        )
    else:
        jobs["8"] = job_data(
            8,
            ["run_sp_exp_Ma"],
            ["mask_runner"],
            ["shdus"],
            n_mult=[1],
            path_main=path_main,
            path_left="exp_runs",
            output_subdirs= "shdus",
            path_right="output",
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

    jobs["32"] = job_data(
        32,
        [
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            #"run_sp_exp_Pi"
        ],
        [
            "sextractor_runner",
            "setools_runner",
            "setools_runner",
            "setools_runner",
            "psfex_runner",
            # "psfex_interp_runner"],
        ],
        "shdus",
        n_mult=[2, 2, 1, 1, 2],
        path_main=path_main,
        path_left="exp_runs",
        output_subdirs="shdus",
        path_right="output",
        path_output=[
            "output",
            "output/rand_split",
            "output/stat",
            "logs",
            "output",
        ],
        special=[False, False, True, True, False],
        verbose=verbose,
    )

    jobs["64"] = job_data(
        "64",
        "run_sp_tile_PsViSmVi",
        [
            "psfex_interp_runner",
            "psfex_interp_runner",
            "vignetmaker_runner_run_1",
            "spread_model_runner",
            "vignetmaker_runner_run_2",
        ],
        "tile_IDs",
        n_mult=[1, 1, 1, 1, 4],
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        path_output=["output", "logs", "output", "output", "output"],
        special=[False, True, False, False, False],
        verbose=verbose,
    )

    if patch in ("P2", "P5"):
        n_sh = 1
    else:
        n_sh = 8
    run_dirs = [f"run_sp_tile_ngmix_Ng{idx+1}u" for idx in range(n_sh)]
    output_path_missing_IDs = [
        f"{path_main}/summary/missing_job_128_ngmix_runner_{idx+1}.txt"
        for idx in range(n_sh)
    ]
    jobs["128"] = job_data(
        "128",
        run_dirs,
        ["ngmix_runner"] * n_sh,
        "tile_IDs",
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        output_path_missing_IDs=output_path_missing_IDs,
        verbose=verbose,
    )

    jobs["256"] = job_data(
        "256",
        ["run_sp_Ms"],
        ["merge_sep_cats_runner"],
        "tile_IDs",
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    jobs["512"] = job_data(
        "512",
        ["run_sp_Mc"],
        ["make_cat_runner"],
        "tile_IDs",
        path_main=path_main,
        path_left="tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    # Post-processing
    jobs["1024"] = job_data(
        "1024",
        ["run_sp_combined_final"],
        ["make_catalog_runner"],
        "tile_IDs",
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    jobs["2048"] = job_data(
        "2048",
        "run_sp_combined_psf",
        ["psfex_interp_runner"],
        "shdus",
        path_main=path_main,
        path_left="output",
        verbose=verbose,
    )

    return jobs, list_tile_IDs_dot


