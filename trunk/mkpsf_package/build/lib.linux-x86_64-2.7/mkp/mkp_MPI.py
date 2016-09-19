#!/usr/bin/env python

# -- Python imports
import sys
from mpi4py import MPI  # MPI interface

# -- External Modules
from mpfg3.mpfg3_MPI import *       # mpfg3: multiprocessing framework MPI mgt

# -- Module-specific imports
import mkp_args as arg             # command-line arguments and options
import mkp_job as jobp             # job processing
from mkp_help import *             # helper utility functions


"""! 
   mkp_MPI.py - Multiprocessing with MPI - Parallel SExtractor
"""


# -------------------------------------------------------------------------------------------------
class MkpMasterMPI(Mpfg3MasterMPI):
   
   """! Master calculator for MPI: distribute jobs to workers and collect/process their results """

   def __init__(self, args, comm):
      """! 
         Master class constructor 
         @param args list of command-line options
         @param comm communication channel   
      """

      # --- Master process
      Mpfg3MasterMPI.__init__(self, args, comm)

      # --- Job Processor
      self.job_processor = jobp.MkpJobProcessor(self)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_job_processor(self):
      """
         Factory method for creating a Job Processor for managing the life-cycle of jobs. 
         @return instance of job processor object
      """

      # --- Decide if should create a Mono- or Multi-Epoch job processor
      return jobp.MkpJobProcessor(self)
      

# -------------------------------------------------------------------------------------------------
class MkpWorkerMPI(Mpfg3WorkerMPI):      

   def __init__(self, args, comm, rank):
      """! Worker class constructor          
         @param args list of command-line options
         @param comm comunication channel
         @param rank rank of MPI process (> 0)
     """
      Mpfg3WorkerMPI.__init__(self, args, comm, rank)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point, either for the master or for a particular worker """

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()                 # rank: 0 for the main process, > 0 for workers

   # --- Command-line arguments and options
   args = arg.MkpArgs()

   if rank == 0:

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         args.print_usage(args.helper.get_main_basename())
         sys.exit(0)

      # --- Master process
      master = MkpMasterMPI(args, comm)
      if master is None:
         sys.exit(2)
      else:
         try:
            master.run()
         except:
            MkpHelper.print_error("The master process ended unexpectedly")
         finally:
            master.shutdown()
            sys.exit(0)

   else:

      # --- Check if Help requested
      if args.options["-h"] or args.options["--help"]:
         sys.exit(0)

      # --- Worker process
      worker = MkpWorkerMPI(args, comm, rank)  
      if worker is None:
         sys.exit(2)
      else:
         try:
            worker.run()
         except:
            MkpHelper.print_error("The worker process {0} ended unexpectedly".format(worker.name))
         finally:
            sys.exit(0)

# -- EOF mkp_MPI.py
