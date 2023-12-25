#!/usr/bin/env python

import sys

from shapepipe.utilities.summary import *


def main(argv=None):


    # Setting
    verbose = True
    log_file_name = "summary_log.txt"
    handlers = [logging.FileHandler(log_file_name), logging.StreamHandler()]
    logging.basicConfig(
        level=logging.INFO, format="%(message)s", handlers=handlers
    )

    main_dir = "."
    retrieve = "vos"
    tile_ID_path = f"{main_dir}/tile_numbers.txt"

    # tile IDs with dots
    list_tile_IDs_dot = get_IDs_from_file(tile_ID_path)

    # tile IDs with dashes
    list_tile_IDs = job_data.replace_dot_dash(list_tile_IDs_dot)
    n_tile_IDs = len(list_tile_IDs)
    n_CCD = 40

    par_runtime = {}

    par_runtime["n_tile_IDs"] = n_tile_IDs
    par_runtime["list_tile_IDs"] = list_tile_IDs

    jobs = {}

    if retrieve == "vos":
        n_link = 2
    else:
        n_link = 1

    jobs["1"] = job_data(
        1,
        "run_sp_GitFeGie_",
        [
            "get_images_runner_run_1",
            "find_exposures_runner",
            "get_images_runner_run_2",
        ],
        ["tile_IDs", "tile_IDs", "exposures"],
        pattern=["CFIS_", "", ""],
        n_mult=[1 * n_link, 1, 3],
        path_left=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["2"] = job_data(
        2,
        ["run_sp_Uz", "run_sp_exp_SpMh", "run_sp_exp_SpMh_2023-12"],
        ["uncompress_fits_runner", "merge_headers_runner", "split_exp_runner"],
        ["tile_IDs", 0, "3*n_shdus+n_exposures"],
        n_mult=[1, 1, 1],
        path_left=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["4"] = job_data(
        4,
        "run_sp_combined_flag",
        ["mask_runner_run_1"],
        ["tile_IDs"],
        path_left=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["8"] = job_data(
        8,
        "run_sp_combined_flag",
        ["mask_runner_run_2"],
        ["shdus"],
        path_left=f"{main_dir}/output",
        verbose=verbose,
    )

    jobs["16"] = job_data(
        16,
        "run_sp_tile_Sx",
        ["sextractor_runner"],
        ["tile_IDs"],
        n_mult=2,
        path_left=f"{main_dir}/tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    # TODO setools_runner output/rand_split
    jobs["32"] = job_data(
        32,
        [
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
            "run_sp_exp_SxSePsf",
        ],  # "run_sp_exp_Pi"],
        [
            "sextractor_runner",
            "setools_runner",
            "psfex_runner",
        ],  # "psfex_interp_runner"],
        "shdus",
        n_mult=[2, 2, 2],  # 1],
        path_left=f"{main_dir}/exp_runs",
        output_subdirs="shdus",
        path_right="output",
        verbose=verbose,
    )

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
        path_left=f"{main_dir}/tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        verbose=verbose,
    )

    n_sh = 8
    run_dirs = [f"run_sp_tile_ngmix_Ng{idx+1}u" for idx in range(n_sh)]
    output_path_missing_IDs = [
        f"missing_job_128_ngmix_runner_{idx+1}.txt" for idx in range(n_sh)
    ]
    jobs["128"] = job_data(
        "128",
        run_dirs,
        ["ngmix_runner"] * 8,
        "tile_IDs",
        path_left=f"{main_dir}/tile_runs",
        output_subdirs=[f"{tile_ID}/output" for tile_ID in list_tile_IDs_dot],
        output_path_missing_IDs=output_path_missing_IDs,
        verbose=verbose,
    )


    jobs["1024"] = job_data(
        "1024",
        "run_sp_combined_psf",
        ["psfex_interp_runner"],
        "shdus",
        path_left=f"{main_dir}/output",
        verbose=verbose
    )

    job_data.print_stats_header()

    for key in "1":
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime, indices=[0, 1])

        all_exposures = get_all_exposures(job._paths_in_dir[1], verbose=True)
        par_runtime["n_exposures"] = len(all_exposures)

        job.check_numbers(par_runtime, indices=[2])

    # Update runtime parameter
    par_runtime["n_shdus"] = get_par_runtime(par_runtime, "exposures") * n_CCD
    par_runtime["n_3*n_shdus+n_exposures"] = 3 * get_par_runtime(
        par_runtime, "shdus"
    ) + get_par_runtime(par_runtime, "exposures")
    par_runtime["list_shdus"] = get_all_shdus(all_exposures, n_CCD)

    print_par_runtime(par_runtime, verbose=verbose)

    #for key in ["2", "4", "8", "16", "32", "64", "128"]:
    for key in ["1024"]:
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
