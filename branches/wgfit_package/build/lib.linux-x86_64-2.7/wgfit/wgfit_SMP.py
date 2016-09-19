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
import wgfit_args as arg             # command-line arguments and options
import wgfit_job as jobp             # job processing
from wgfit_helper import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class WGfitMasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      try:

         # --- Command-line arguments and options
         args = arg.WGfitArgs()

         # --- Master process
         MasterSMP.__init__(self, args) 

         # --- Check if Help requested
         if args.options["-h"] or args.options["--help"]:
            args.print_usage(args.helper.get_main_basename())
            sys.exit(0)

         # --- Job Processor
         self.job_processor = jobp.WGfitJobProcessor(self)

         # --- Helper methods
         self._helper = WGfitHelper() 

         # --- Show config_summary
         self.helper.show_config_summary(self)

         # --- Save a copy of the wgfit configuration file to the log directory
         shutil.copy(os.path.join(args.options["-d"], args.options["-c"]), self.log_output_dir)

      except Exception as detail:

         print(detail)

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def wgfit_version(self):
      """! @return this version of the wgfit code """
      return self._wgfit_version

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_worker(self, arg_options):
      """ Create a worker process (Override) """
      return WGfitWorkerSMP(arg_options)

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
      return jobp.WGfitJobProcessor(self)


# -------------------------------------------------------------------------------------------------
class WGfitWorkerSMP(MpfxWorkerSMP):
   
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
   """ Main entry point of the wgfit master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = WGfitMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except Exception as detail:
         print_error("The wgfit master process ended unexpectedly: {0}", detail)
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF wgfit_SMP.py
