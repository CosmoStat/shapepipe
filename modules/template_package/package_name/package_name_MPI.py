#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Package MPI Script

This module contain a class for executing the package using SMP.

:Authors: Samuel Farrens and Marc Gentile

:Date: 10/04/2018

"""

# -- Python imports
from sys import exit
from os.path import join
from shutil import copy
from mpi4py import MPI

# -- External Modules
from mpfx import mpfx_MPI
from shapepipe_base.args import PackageArgs
from shapepipe_base.helper import PackageHelper

# -- Module-specific imports
from info import __version__, __whoami__, __python_depend__, __system_depend__
from job import PackageJobProcessor


class PackageMasterMPI(mpfx_MPI.MpfxMasterMPI):

    """Package Master MPI

    This class contains methods for distributing jobs to workers and
    collecting and processing their results using MPI.

    Parameters
    ----------
    args : class
        PackageArgs instance
    comm : class
        MPI COMM_WORLD instance

    """

    def __init__(self, args, comm):

        try:

            # --- Master process
            mpfx_MPI.MasterMPI.__init__(self, args, comm)

            # --- Helper methods
            self._helper = PackageHelper(__version__, __whoami__,
                                         __python_depend__, __system_depend__)

            # --- Show config_summary
            self.helper.show_config_summary(self)

            # --- Job Processor
            self.job_processor = PackageJobProcessor(self)

            # --- Save a copy of the gfit configuration file to the log dir
            copy(join(args.options['-d'], args.options['-c']),
                 self.log_output_dir)

        except Exception as detail:

            print(detail)

    def create_job_processor(self):

        """Create Job Processor

        This method creats a Job Processor instance for managing the
        life-cycle of jobs.

        Returns
        -------
        class instance of the job processor object

        """

        return PackageJobProcessor(self)


def main():

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # --- Main process
    try:

        # --- Command-line arguments and options
        args = PackageArgs()

        # --- Check if Help requested
        if args.options["-h"] or args.options["--help"]:
            args.print_usage(args.helper.get_main_basename())
            exit(0)

        if rank == 0:
            # --- Master process
            master = PackageMasterMPI(args, comm)
        else:
            # --- Worker process
            master = PackageWorkerMPI(args, comm, rank)

        if isinstance(master, type(None)):
            exit(2)

        else:
            try:
                master.run()

            except Exception as detail:
                if rank == 0:
                    PackageHelper.print_error('The master process ended '
                                              'unexpectedly '
                                              '({0})'.format(detail))
                else:
                    PackageHelper.print_error('The worker process {0} ended '
                                              'unexpectedly ({1})'.format(
                                               worker.name, detail))

            finally:
                if rank == 0:
                    master.shutdown()
                exit(0)

    except Exception as detail:
        print(detail)


if __name__ == '__main__':
    main()
