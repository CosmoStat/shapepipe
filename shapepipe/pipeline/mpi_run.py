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


def submit_mpi_jobs(jobs, filehd, config, timeout, module, verbose):
    """ Submit MPI Jobs

    This method distributes the jobs to the workers using MPI.

    """

    result = []

    for job in jobs:
        wh = WorkerHandler(verbose=verbose)
        result.append(wh.worker(job[0], job[1], filehd, config, timeout,
                      module))

    return result
