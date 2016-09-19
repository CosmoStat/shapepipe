#!/usr/bin/env python

# -- Python imports
import sys

# -- External Modules
#from mpf.mp_SMP import *          # mpf: multiprocessing framework, SMP mgt
from mpfg.mp_calc_SMP import *     # mpf: multiprocessing framework, SMP mgt
from mpfx.mpfx_SMP import *        # mpf: multiprocessing framework, SMP mgt

# -- Module-specific imports     
import mpfg3_args as arg           # command-line arguments and options
import mpfg3_job as jobp           # job processing
import mpfg.mp_helper              # helper functions


"""! 
   mpfg3_SMP.py - Multiprocessing with SMP 
"""

# -------------------------------------------------------------------------------------------------
class Mpfg3MasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      # --- Command-line arguments and options
      args = arg.Mpfg3Args()

      # --- Master process
      MasterSMP.__init__(self, args)    

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(self.helper.get_main_basename())
         sys.exit(0)


#      # --- Job Processor
#      self.job_processor = jobp.Mpfg3JobProcessor(self)

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
      """! 
         Create a worker process (Override)
         @param arg_options list of command-line options
      """
      return Mpfg3WorkerSMP(arg_options)

   # -----------------------------------------------------------------------------------------------
   def shutdown(self):
      """! Shutdown master process: close files, stop processes, etc. (Override) """
      MpfxMasterSMP.shutdown(self)

   # -----------------------------------------------------------------------------------------------
   def create_job_processor(self):
      """
         Factory method for creating a Job Processor for managing the life-cycle of jobs. 
         @return instance of job processor object
      """

      # --- Decide if should create a Mono- or Multi-Epoch job processor
      if self.config.has_section("JOB_MANAGEMENT") and \
         self.config.has_key("MULTIPLE_EPOCHS_PER_JOB", "JOB_MANAGEMENT"): 

         if self.config.get_as_boolean("MULTIPLE_EPOCHS_PER_JOB", "JOB_MANAGEMENT"):
            return jobp.Mpfg3MultiEpochJobProcessor(self)
         else:
            return jobp.Mpfg3JobProcessor(self)
      else:
         return jobp.Mpfg3JobProcessor(self)


# -------------------------------------------------------------------------------------------------
class Mpfg3WorkerSMP(MpfxWorkerSMP):
   
   def __init__(self, arg_options):
      """!
         Worker class constructor 
         @param arg_options list of command-line options
      """

      # --- Worker process
      MpfxWorkerSMP.__init__(self, arg_options)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point of the Mpfg3 master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = Mpfg3MasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         mpfg.mp_helper.Helper.print_error("The Mpfg3 master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF mpfg3_SMP.py
