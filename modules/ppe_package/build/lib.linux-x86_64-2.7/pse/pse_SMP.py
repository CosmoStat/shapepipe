#!/usr/bin/env python

# -- Python imports
import sys
import shutil
from operator import itemgetter, attrgetter

# -- External Modules
from mpfx.mpfx_SMP import *        # multiprocessing framework SMP mgt

# -- Module-specific imports     
import pse_args as arg             # command-line arguments and options
import pse_job as jobp             # job processing
from pse_help import *             # helper utility functions

"""! 
   pse_SMP.py - Multiprocessing with SMP - Parallel SExtractor
"""

## -------------------------------------------------------------------------------------------------
#def get_base_master_class():

#   parent_class = MasterSMP(arg.PseArgs())
#   base_type = parent_class.config.get_as_string("TYPE", "DATASET")
#   return eval(base_type + "MasterSMP")

# -------------------------------------------------------------------------------------------------
class PseMasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      try:
      
         # --- Command-line arguments and options
         args = arg.PseArgs()

         # --- Master process
         MasterSMP.__init__(self, args) 

         # --- Check if Help requested
         if args.options["-h"] or args.options["--help"]:
            args.print_usage(args.helper.get_main_basename())
            sys.exit(0)

         # --- Job Processor
         self.job_processor = jobp.PseJobProcessor(self)

         # --- Helper methods
         self._helper = PseHelper() 

         # --- Show config_summary
         self.helper.show_config_summary(self)

         # --- Save a copy of the gfit configuration file to the log directory
         shutil.copy(os.path.join(args.options["-d"], args.options["-c"]), self.log_output_dir)

      except Exception as detail:

         print(detail)


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_run_output_dir(self):
      """! 
         Create run directory where output data will be stored
         @return time-stamped output directory tracing the run of the process 
      """

      return  MpfxMasterSMP.create_run_output_dir(self)

   # -----------------------------------------------------------------------------------------------
   def create_worker(self, arg_options):
      """ Create a worker process (Override) """
      return PseWorkerSMP(arg_options)

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

      return jobp.PseJobProcessor(self)


## -------------------------------------------------------------------------------------------------
#def get_base_worker_class():

#   parent_class = MasterSMP(arg.PseArgs())
#   base_type = parent_class.config.get_as_string("TYPE", "DATASET")
#   return eval(base_type + "WorkerSMP")


# -------------------------------------------------------------------------------------------------
class PseWorkerSMP(MpfxWorkerSMP):
   
   def __init__(self, arg_options):
      """ Worker initialization """

      # --- Worker process
      MpfxWorkerSMP.__init__(self, arg_options)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """ Main entry point of the PSE master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = PseMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         PseHelper.print_error("The PSE master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF pse_SMP.py
