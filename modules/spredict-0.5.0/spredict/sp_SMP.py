#!/usr/bin/env python

# -- Python imports
import sys
from operator import itemgetter, attrgetter

# -- External Modules
from mpfg3.mpfg3_SMP import *       # mpf: multiprocessing framework, SMP mgt


# -- Module-specific imports
import sp_args as arg             # command-line arguments and options
import sp_job as jobp             # job processing
from sp_helper import *           # helper utility functions

"""!
   sp_SMP.py - spredict for GREAT3 - Multiprocessing with SMP
"""

# -------------------------------------------------------------------------------------------------
class SpredictMasterSMP(Mpfg3MasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor
      """

      # --- Command-line arguments and options
      args = arg.SpredictArgs()

      # --- Master process
      MasterSMP.__init__(self, args)

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(args.helper.get_main_basename())
         sys.exit(0)

      # --- Job Processor
      self.job_processor = jobp.SpredictJobProcessor(self)

      # --- Helper methods
      self._helper = SpredictHelper()

      # --- Show config_summary
      self.helper.show_config_summary(self)

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   @property
   def spredict_version(self):
      """! @return this version of the spredict code """
      return self._spredict_version

   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~

   def create_worker(self, arg_options):
      """ Create a worker process (Override) """
      return SpredictWorkerSMP(arg_options)


   def shutdown(self):
      """ Shutdown master process: close files, stop processes, etc. (Override) """
      MasterSMP.shutdown(self)


# -------------------------------------------------------------------------------------------------
class SpredictWorkerSMP(Mpfg3WorkerSMP):

   def __init__(self, arg_options):
      """ Worker initialization """

      # --- Worker process
      Mpfg3WorkerSMP.__init__(self, arg_options)


   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """ Main entry point of the spredict master process """

   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = SpredictMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         print_error("The spredict master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF sp_SMP.py
