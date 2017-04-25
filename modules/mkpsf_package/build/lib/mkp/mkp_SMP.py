#!/usr/bin/env python

# -- Python imports
import sys
from operator import itemgetter, attrgetter

# -- External Modules
from mpfg3.mpfg3_SMP import *       # mpfg3: multiprocessing framework SMP mgt

# -- Module-specific imports     
import mkp_args as arg             # command-line arguments and options
import mkp_job as jobp             # job processing
from mkp_help import *             # helper utility functions


"""! 
   mkp_SMP.py - Multiprocessing with SMP - Parallel SExtractor
"""

# -------------------------------------------------------------------------------------------------
class MkpMasterSMP(Mpfg3MasterSMP):

   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self):
      """!
         Master default class constructor 
      """
      # --- Command-line arguments and options
      args = arg.MkpArgs()

      # --- Master process
      MasterSMP.__init__(self, args) 

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(args.helper.get_main_basename())
         sys.exit(0)

      # --- Job Processor
      self.job_processor = jobp.MkpJobProcessor(self)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   def create_worker(self, arg_options):
      """ Create a worker process (Override) """
      return MkpWorkerSMP(arg_options)


   def shutdown(self):
      """ Shutdown master process: close files, stop processes, etc. (Override) """
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

         return jobp.Mpfg3JobProcessor

# -------------------------------------------------------------------------------------------------
class MkpWorkerSMP(Mpfg3WorkerSMP):
   
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
   """ Main entry point of the PSE master process """
 
   # --- Invoke the master calculator in module: mp_calc_SMP.py
   master = MkpMasterSMP()
   if master is None:
      sys.exit(2)
   else:
      try:
         master.run()
      except:
         MkpHelper.print_error("The PSE master process ended unexpectedly")
      finally:
         master.shutdown()
         sys.exit(0)

# -- EOF mkp_SMP.py
