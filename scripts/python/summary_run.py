#!/usr/bin/env python

import sys
import os

from shapepipe.utilities.summary import *


def main(argv=None):

    patch = argv[1]

    retrieve = "vos"
    verbose = False
    
    import summary_params_pre_v2 as params

    jobs, list_tile_IDs_dot = params.set_jobs_v2_pre_v2(patch, retrieve, verbose)
    
    list_tile_IDs = job_data.replace_dot_dash(list_tile_IDs_dot)
    
    # Numbers updated at runtime
    par_runtime = params.init_par_runtime(list_tile_IDs)
    
    job_data.print_stats_header()

    for key in "1":
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime, indices=[0, 1])

        all_exposures = get_all_exposures(job._paths_in_dir[1], verbose=True)
        par_runtime["n_exposures"] = len(all_exposures)
        par_runtime["list_exposures"] = all_exposures

        job.check_numbers(par_runtime, indices=[2])

    par_runtime = params.update_par_runtime_after_find_exp(par_runtime, all_exposures)
    
    print_par_runtime(par_runtime, verbose=verbose)

    # Get all keys after "1"
    keys = sorted(jobs.keys(), key=int)
    _ = keys.pop(0)

    for key in keys:
        job = jobs[key]
        job.print_intro()
        job.check_numbers(par_runtime=par_runtime)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
