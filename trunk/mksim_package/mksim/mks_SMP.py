#!/usr/bin/env python

# -- Python imports
import sys
import shutil

# -- External Modules
from mpfx.mpfx_SMP import *        # multiprocessing framework SMP mgt

# -- Module-specific imports     
import mks_args as arg             # command-line arguments and options
import mks_job as jobp             # job processing
import mks_help as hlp             # helper utility functions

"""! 
   mks_SMP.py - Multiprocessing with SMP - MKSIM
"""

## -------------------------------------------------------------------------------------------------
#def get_base_master_class():

#   parent_class = MasterSMP(arg.MksArgs())
#   base_type = parent_class.config.get_as_string("TYPE", "DATASET")
#   return eval(base_type + "MasterSMP")

# -------------------------------------------------------------------------------------------------
class MksMasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      try:
      
         # --- Command-line arguments and options
         args = arg.MksArgs()

         # --- Master process
         MasterSMP.__init__(self, args) 

         # --- Check if Help requested
         if args.options["-h"] or args.options["--help"]:
            args.print_usage(args.helper.get_main_basename())
            sys.exit(0)

         # --- Helper methods
         self._helper = hlp.MksHelper()

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
      return MksWorkerSMP(arg_options)

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
      return jobp.MksJobProcessor(self)


## -------------------------------------------------------------------------------------------------
#def get_base_worker_class():

#   parent_class = MasterSMP(arg.MksArgs())
#   base_type = parent_class.config.get_as_string("TYPE", "DATASET")
#   return eval(base_type + "WorkerSMP")


# -------------------------------------------------------------------------------------------------
class MksWorkerSMP(MpfxWorkerSMP):
   
   def __init__(self, arg_options):
      """ Worker initialization """

      # --- Worker process
      MpfxWorkerSMP.__init__(self, arg_options)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """ Main entry point of the MKSIM master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = MksMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         hlp.MksHelper.print_error("The MKSIM master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF mks_SMP.py
