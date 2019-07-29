# -*- coding: utf-8 -*-

from shapepipe.pipeline.worker_handler import WorkerHandler


def split_mpi_jobs(jobs, batch_size):
    """ Split MPI Jobs

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


def submit_mpi_jobs(jobs, config, timeout, run_dirs, module_runner,
                    worker_log, verbose):
    """ Submit MPI Jobs

    This method distributes the jobs to the workers using MPI.

    """

    result = []

    for process in jobs:

        w_log_name = worker_log(module_runner.__name__, process[0])

        wh = WorkerHandler(verbose=verbose)
        result.append(wh.worker(process[1:], process[0], w_log_name,
                      run_dirs, config, timeout, module_runner))

    return result
