#!/usr/bin/env python

# -- Python imports
import sys
import shutil

# -- External Modules
from mpfx.mpfx_SMP import *

# -- Module-specific imports
import ppe_args as arg
import ppe_job as jobp
from ppe_help import PpeHelper

"""!
    ppe_SMP.py - Multiprocessing with SMP - Parallel PSFEx
"""


# ----------------------------------------------------------------------------
class PpeMasterSMP(MpfxMasterSMP):

    """!
    Master calculator for SMP: distribute jobs to workers and collect/process
    their results

    """

    def __init__(self):
        """!
        Master default class constructor

        """

        try:

            # --- Command-line arguments and options
            args = arg.PpeArgs()

            # --- Master process
            MasterSMP.__init__(self, args)

            # --- Check if Help requested
            if args.options['-h'] or args.options['--help']:
                args.print_usage(args.helper.get_main_basename())
                sys.exit(0)

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

    # ~~~~~~~~~~~~~~~
    # Public methods
    # ~~~~~~~~~~~~~~~

    # ------------------------------------------------------------------------
    def create_run_output_dir(self):
        """!
        Create run directory where output data will be stored
        @return time-stamped output directory tracing the run of the process

        """

        return MpfxMasterSMP.create_run_output_dir(self)

    # ------------------------------------------------------------------------
    def create_worker(self, arg_options):
        """
        Create a worker process (Override)

        """

        return PpeWorkerSMP(arg_options)

    # ------------------------------------------------------------------------
    def shutdown(self):
        """
        Shutdown master process: close files, stop processes, etc. (Override)

        """

        MasterSMP.shutdown(self)

    # ------------------------------------------------------------------------
    def create_job_processor(self):
        """
        Factory method for creating a Job Processor for managing the
        life-cycle of jobs.
        @return instance of job processor object

        """

        return jobp.PpeJobProcessor(self)


# ----------------------------------------------------------------------------
class PpeWorkerSMP(MpfxWorkerSMP):

    def __init__(self, arg_options):
        """
        Worker initialization

        """

        # --- Worker process
        MpfxWorkerSMP.__init__(self, arg_options)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    # --- Invoke the master calculator in module: mp_calc_SMP.py
    master = PpeMasterSMP()

    if master is None:
        sys.exit(2)
    else:
        try:
            print 'before' + '-' * 100
            master.run()
        except Exception:
            PpeHelper.print_error("The PPE master process ended unexpectedly")
        finally:
            master.shutdown()
            sys.exit(0)

# -- EOF ppe_SMP.py
