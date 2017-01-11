#!/usr/bin/env python

# -- Python imports
import sys

# -- External Modules
from mpf.mp_calc_SMP import *       # mpf: multiprocessing framework, SMP mgt
from mpfx.mpfx_SMP import *       # mpf: multiprocessing framework, SMP mgt

# -- Module-specific imports     
import mpfcs82_args as arg             # command-line arguments and options
import mpfcs82_job as jobp             # job processing


"""! 
   mpfcs82_SMP.py - Multiprocessing with SMP 
"""

# -------------------------------------------------------------------------------------------------
class MpfCS82MasterSMP(MpfxMasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """

      # --- Command-line arguments and options
      args = arg.MpfCS82Args()

      # --- Master process
      MasterSMP.__init__(self, args)    

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(self.helper.get_main_basename())
         sys.exit(0)

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
      return MpfCS82WorkerSMP(arg_options)

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
            return jobp.MpfCS82MultiEpochJobProcessor(self)
         else:
            return jobp.MpfCS82JobProcessor(self)
      else:
         return jobp.MpfCS82JobProcessor(self)


# -------------------------------------------------------------------------------------------------
class MpfCS82WorkerSMP(MpfxWorkerSMP):
   
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
   """! Main entry point of the MpfCS82 master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = MpfCS82MasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         print_error("The MpfCS82 master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF mpfcs82_SMP.py
