#!/usr/bin/env python

import argparse
import sys
import os

from shapepipe.utilities import summary

from shapepipe.utilities import summary_params_pre_v2 as summary_params


def run(patch, job_exclusive=None, verbose=False):

    jobs, list_tile_IDs_dot = summary_params.set_jobs_v2_pre_v2(patch, verbose)

    list_tile_IDs = summary.job_data.replace_dot_dash(list_tile_IDs_dot)

    # Numbers updated at runtime
    par_runtime = summary.init_par_runtime(list_tile_IDs)

    summary.job_data.print_stats_header()

    exp_IDs_path = "exp_numbers.txt"
    if os.path.exists(exp_IDs_path):
        # Read exposure ID list if file exists
        all_exposures = summary.get_IDs_from_file(exp_IDs_path)
        par_runtime = summary.update_par_runtime_after_find_exp(
            par_runtime, all_exposures
        )

    if (
        not os.path.exists(exp_IDs_path)
        or not job_exclusive
        or int(job_exclusive) & 1
    ):
        # Run job 1 if exposure ID list file does not exist or
        # job_exclusive is 1 or not set
        key = "1"
        jobs[key].print_intro()
        jobs[key].check_numbers(par_runtime=par_runtime, indices=[0, 1])

        all_exposures = summary.get_all_exposures(
            jobs[key]._paths_in_dir[1], verbose=True
        )
        par_runtime = summary.update_par_runtime_after_find_exp(
            par_runtime, all_exposures
        )

        jobs[key].write_IDs_to_file("exp_numbers.txt", all_exposures)

        jobs[key].check_numbers(par_runtime, indices=[2])

    summary.print_par_runtime(par_runtime, verbose=verbose)

    # Get all keys after "1"
    keys = sorted(jobs.keys(), key=int)
    _ = keys.pop(0)

    for key in keys:
        if job_exclusive and not int(key) & int(job_exclusive):
            continue
        jobs[key].print_intro()
        jobs[key].check_numbers(par_runtime=par_runtime)

    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            'Print summary of ShapePipe job status for a given UNIONS patch.'
        ),
    )
    parser.add_argument(
        'patch',
        help='Patch identifier (e.g. P3, P9).',
    )
    parser.add_argument(
        'job_exclusive',
        nargs='?',
        default=None,
        help=(
            'Bitmask selecting which jobs to run; '
            + 'if omitted, all jobs run.'
        ),
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help='Verbose output.',
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    run(args.patch, args.job_exclusive, args.verbose)


if __name__ == "__main__":
    sys.exit(main())
