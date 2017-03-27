#!/usr/bin/env python

# -- Python imports
import sys
import shutil
from mpi4py import MPI  # MPI interface

# -- External Modules
from mpfx.mpfx_MPI import *        # mpf: multiprocessing framework MPI mgt

# -- Module-specific imports
import ppe_args as arg             # command-line arguments and options
import ppe_job as jobp             # job processing
from ppe_help import *             # helper utility functions

"""!
   ppe_MPI.py - Multiprocessing with MPI - Parallel PSFEx
"""


# ----------------------------------------------------------------------------
class PpeMasterMPI(MpfxMasterMPI):

    """!
    Master calculator for MPI: distribute jobs to workers and collect/process
    their results

    """

    def __init__(self, args, comm):
        """!
        Master class constructor
        @param args list of command-line options
        @param comm communication channel

        """

        try:

            # --- Master process
            MasterMPI.__init__(self, args, comm)

            # --- Job Processor
            self.job_processor = jobp.PpeJobProcessor(self)

            # --- Helper methods
            self._helper = PpeHelper()

            # --- Show config_summary
            self.helper.show_config_summary(self)

            # --- Save a copy of the gfit configuration file to the log dir
            shutil.copy(os.path.join(args.options['-d'], args.options['-c']),
                        self.log_output_dir)

        except Exception as detail:

            print(detail)

    # ~~~~~~~~~~~
    # Properties
    # ~~~~~~~~~~~

    # ~~~~~~~~~~~~~~~
    # Public methods
    # ~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------
    def create_run_output_dir(self):
        """!
        Create run directory where output data will be stored
        @return time-stamped output directory tracing the run of the process

        """

        return MpfxMasterMPI.create_run_output_dir(self)

    # ------------------------------------------------------------------------
    def create_job_processor(self):
        """
        Factory method for creating a Job Processor for managing the
        life-cycle of jobs.
        @return instance of job processor object

        """

        return jobp.PpeJobProcessor(self)


# ----------------------------------------------------------------------------
class PpeWorkerMPI(WorkerMPI):

    def __init__(self, args, comm, rank):
        """! Worker class constructor
        @param args list of command-line options
        @param comm comunication channel
        @param rank rank of MPI process (> 0)

        """

        WorkerMPI.__init__(self, args, comm, rank)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:

        # --- Main process
        try:

            # --- Command-line arguments and options
            args = arg.PpeArgs()

            # --- Check if Help requested
            if args.options["-h"] or args.options["--help"]:
                args.print_usage(args.helper.get_main_basename())
                sys.exit(0)

            # --- Master process
            master = PpeMasterMPI(args, comm)
            if master is None:
                sys.exit(2)
            else:
                try:
                    master.run()
                except Exception as detail:
                    PpeHelper.print_error('The master process ended '
                                          'unexpectedly ({0})'.format(detail))
                finally:
                    master.shutdown()
                    sys.exit(0)

        except Exception as detail:
            print(detail)

    else:

        # --- Worker process
        try:

            # --- Command-line arguments and options
            args = arg.PpeArgs()

            # --- Check if Help requested
            if args.options["-h"] or args.options["--help"]:
                sys.exit(0)

            # --- Worker process
            worker = PpeWorkerMPI(args, comm, rank)
            if worker is None:
                sys.exit(2)
            else:
                try:
                    worker.run()
                except Exception as detail:
                    PpeHelper.print_error('The worker process {0} ended '
                                          'unexpectedly ({1})'.format(
                                           worker.name, detail))
                finally:
                    sys.exit(0)

        except Exception as detail:
            pass

# -- EOF ppe_MPI.py
