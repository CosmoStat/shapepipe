# -*- coding: utf-8 -*-

from shapepipe.pipeline.worker_handler import WorkerHandler


def mpi_run(job, filehd, config, timeout, module, verbose):
    """ Distribute MPI Jobs

    This method distributes the jobs to the workers using MPI.

    """

    wh = WorkerHandler(verbose=verbose)

    result = wh.worker(job[0], job[1], filehd, config, timeout, module)

    return result
