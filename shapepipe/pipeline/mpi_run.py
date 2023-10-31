"""MPI RUN.

This module includes functions for handling the distribution of MPI jobs.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from shapepipe.pipeline.worker_handler import WorkerHandler


def split_mpi_jobs(jobs, batch_size):
    """Split MPI Jobs.

    Split the number of MPI jobs over the number of processes.

    Parameters
    ----------
    jobs : list
        List of MPI jobs
    batch_size : int
        Batch size

    Returns
    -------
    list
        Split list of jobs

    """
    return [jobs[_i::batch_size] for _i in range(batch_size)]


def submit_mpi_jobs(
    jobs,
    config,
    timeout,
    run_dirs,
    module_runner,
    module_config_sec,
    worker_log,
    verbose,
):
    """Submit MPI Jobs.

    This method distributes the jobs to the workers using MPI.

    """
    result = []

    for process in jobs:

        w_log_name = worker_log(module_runner.__name__, process[0])

        wh = WorkerHandler(verbose=verbose)
        print("MKDEBUG in submit_mpi_job")
        print("MKDEBUG process[1:] = ", process[1:])
        print("MKDEBUG process[0] = ", process[0])
        print("MKDEBUG w_log_name = ", w_log_name)
        print("MKDEBUG run_dirs = ", run_dirs)
        print("MKDEBUG config = ", config)
        print("MKDEBUG timeout = ", timeout)
        print("MKDEBUG module_runner = ", module_runner)
        print("MKDEBUG module_config_sec = ", module_config_sec)
        result.append(wh.worker(
            process[1:],
            process[0],
            w_log_name,
            run_dirs,
            config,
            module_config_sec,
            timeout,
            module_runner
        ))

    return result
