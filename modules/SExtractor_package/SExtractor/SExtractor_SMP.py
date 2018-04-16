#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Package SMP Script

This module contain a class for executing the package using SMP.

:Authors: Samuel Farrens and Marc Gentile

:Date: 10/04/2018

"""

# -- Python imports
from sys import exit
from os.path import join
from shutil import copy

# -- External Modules
from mpfx import mpfx_SMP
from shapepipe_base.args import PackageArgs
from shapepipe_base.helper import PackageHelper

# -- Module-specific imports
from info import __version__, __whoami__, __python_depend__, __system_depend__
from job import PackageJobProcessor


class PackageMasterSMP(mpfx_SMP.MpfxMasterSMP):

    """Package Master SMP

    This class contains methods for distributing jobs to workers and
    collecting and processing their results using SMP.

    """

    def __init__(self):

        try:

            # --- Command-line arguments and options
            args = PackageArgs()

            # --- Master process
            mpfx_SMP.MasterSMP.__init__(self, args)

            # --- Check if Help requested
            if args.options['-h'] or args.options['--help']:
                args.print_usage(args.helper.get_main_basename())
                exit(0)

            # --- Helper methods
            self._helper = PackageHelper(__version__, __whoami__,
                                         __python_depend__, __system_depend__)

            # --- Show config_summary
            self._helper.show_config_summary(self)

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

    # --- Invoke the master calculator in module: mp_calc_SMP.py
    master = PackageMasterSMP()

    if isinstance(master, type(None)):
        exit(2)

    else:
        try:
            master.run()

        except Exception:
            PackageHelper.print_error("The master process ended unexpectedly")
        finally:
            master.shutdown()
            exit(0)


if __name__ == '__main__':
    main()
