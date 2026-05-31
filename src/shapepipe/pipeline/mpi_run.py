"""MPI RUN.

This module includes functions for handling the distribution of MPI jobs.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os

from shapepipe.pipeline.worker_handler import WorkerHandler


def check_mpi_world(comm):
    """Check MPI World.

    Verify that the MPI world formed at the size the launcher requested, and
    abort loudly otherwise.

    This guards against the "N rank-0 singletons" failure: when the host MPI
    launcher and the container's MPI / PMIx stack are incompatible, each
    process initialises standalone (``COMM_WORLD`` size 1, rank 0). ShapePipe
    would then treat every process as master, hand each the full job list, and
    run that many uncoordinated copies of the pipeline into the same output
    directory -- silently, with exit code 0. Comparing the size that actually
    wired up against the size the launcher intended (``OMPI_COMM_WORLD_SIZE``,
    which is set per process even when the world fails to form) turns that
    silent corruption into an immediate, legible error.

    Parameters
    ----------
    comm : MPI.Comm
        MPI communicator instance (``MPI.COMM_WORLD``)

    Raises
    ------
    RuntimeError
        if the launcher requested more ranks than actually wired up

    """
    intended = os.environ.get("OMPI_COMM_WORLD_SIZE")
    actual = comm.Get_size()

    if intended is not None and int(intended) != actual:
        raise RuntimeError(
            f"MPI world mismatch: the launcher requested {intended} ranks but "
            f"only {actual} wired up (MPI_COMM_WORLD size {actual}). This is "
            f"the 'rank 0 of N singletons' failure -- the host MPI launcher and "
            f"the container's MPI/PMIx stack are incompatible. Aborting rather "
            f"than running {intended} uncoordinated copies of the pipeline into "
            f"the same output directory."
        )


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
    module_config_sec,
    timeout,
    run_dirs,
    module_runner,
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
        result.append(
            wh.worker(
                process[1:],
                process[0],
                w_log_name,
                run_dirs,
                config,
                module_config_sec,
                timeout,
                module_runner,
            )
        )

    return result
