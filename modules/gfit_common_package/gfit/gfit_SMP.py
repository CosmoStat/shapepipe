#!/usr/bin/env python
"""! 
   Multiprocessing with Symmetric Multiprocessing (SMP)
"""

# -- Python imports
import sys
from operator import itemgetter, attrgetter

# -- External Modules
from mpfx.mpfx_SMP import *         # multiprocessing framework, SMP mgt


# -- Module-specific imports     
import gfit_args as arg             # command-line arguments and options
import gfit_job as jobp             # job processing
from gfit_helper import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class GfitMasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      try:

         # --- Command-line arguments and options
         args = arg.GfitArgs()

         # --- Master process
         MasterSMP.__init__(self, args) 

         # --- Check if Help requested
         if args.options["-h"] or args.options["--help"]:
            args.print_usage(args.helper.get_main_basename())
            sys.exit(0)

         # --- Job Processor
         self.job_processor = jobp.GfitJobProcessor(self)

         # --- Helper methods
         self._helper = GfitHelper() 

         # --- Show config_summary
         self.helper.show_config_summary(self)

         # --- Save a copy of the gfit configuration file to the log directory
         shutil.copy(os.path.join(args.options["-d"], args.options["-c"]), self.log_output_dir)

      except Exception as detail:

         print(detail)

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def gfit_version(self):
      """! @return this version of the gfit code """
      return self._gfit_version

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_worker(self, arg_options):
      """ Create a worker process (Override) """
      return GfitWorkerSMP(arg_options)

   # -----------------------------------------------------------------------------------------------
   def shutdown(self):
      """ Shutdown master process: close files, stop processes, etc. (Override) """
      MasterSMP.shutdown(self)

   # -----------------------------------------------------------------------------------------------
   def create_job_processor(self):
      """
         Factory method for creating a Job Processor for managing the life-cycle of jobs. 
         @return instance of job processor object
      """
      return jobp.GfitJobProcessor(self)


# -------------------------------------------------------------------------------------------------
class GfitWorkerSMP(MpfxWorkerSMP):
   
   """! Worker calculator for SMP: process jobs supplied by the master calculator """

   def __init__(self, arg_options):
      """ Worker initialization """

      # --- Worker process
      MpfxWorkerSMP.__init__(self, arg_options)


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """ Main entry point of the gfit master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = GfitMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except Exception as detail:
         print_error("The gfit master process ended unexpectedly: {0}", detail)
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF gfit_SMP.py
