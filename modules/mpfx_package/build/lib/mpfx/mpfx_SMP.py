#!/usr/bin/env python

# -- Python imports
import sys
from operator import itemgetter, attrgetter

# -- External Modules
from mpfg.mp_calc_SMP import *       # mpf: multiprocessing framework, SMP mgt

# -- Module-specific imports     
import mpfx_args as arg             # command-line arguments and options
import mpfx_job as jobp             # job processing
import mpfg.mp_helper               # helper functions


"""! 
   mpfx_SMP.py - Multiprocessing with SMP 
"""

# -------------------------------------------------------------------------------------------------
class MpfxMasterSMP(MasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """
      # --- Command-line arguments and options
      args = arg.MpfxArgs()

      # --- Master process
      MasterSMP.__init__(self, args) 

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(self.helper.get_main_basename())
         sys.exit(0)

#      # --- Job Processor
#      self.job_processor = jobp.MpfxJobProcessor(self)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_run_output_dir(self):
      """! 
         Create run directory where output data will be stored
         @return time-stamped output directory tracing the run of the process 
      """

      run_output_dir = self._build_run_dir_name(self)

      base_output_dir = self.config.get_as_string("BASE_OUTPUT_DIR", "DIR.OUTPUT")
      if not os.path.isdir(base_output_dir):
         self.helper.make_dir(base_output_dir)
      run_output_dir = os.path.join(base_output_dir, run_output_dir)
      self.helper.make_dir(run_output_dir)

      return run_output_dir

   # -----------------------------------------------------------------------------------------------
   def create_worker(self, arg_options):
      """! 
         Create a worker process (Override)
         @param arg_options list of command-line options
      """
      return MpfxWorkerSMP(arg_options)

   # -----------------------------------------------------------------------------------------------
   def shutdown(self):
      """! Shutdown master process: close files, stop processes, etc. (Override) """
      MasterSMP.shutdown(self)

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
            return jobp.MpfxMultiEpochJobProcessor(self)
         else:
            return jobp.MpfxJobProcessor(self)
      else:
         return jobp.MpfxJobProcessor(self)
         

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   def _build_run_dir_name(self, master):

      log_time = time.strftime("%d.%m.%y_%H.%M.%S", time.localtime())
      run_dir_name = 'run_{0}'.format(log_time)

      return run_dir_name

# -------------------------------------------------------------------------------------------------
class MpfxWorkerSMP(WorkerSMP):
   
   def __init__(self, arg_options):
      """!
         Worker class constructor 
         @param arg_options list of command-line options
      """

      # --- Worker process
      WorkerSMP.__init__(self, arg_options)


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point of the Mpfx master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = MpfxMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         mpfg.mp_helper.Helper.print_error("The Mpfx master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF mpfx_SMP.py
