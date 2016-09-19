#!/usr/bin/env python

# -- Python imports
import sys
from mpi4py import MPI  # MPI interface
from operator import itemgetter, attrgetter

# -- External Modules
from mpfx.mpfx_MPI import *      # mpfx: multiprocessing framework MPI mgt

# -- Module-specific imports
import mpfcs82_args as arg             # command-line arguments and options
import mpfcs82_job as jobp             # job processing

"""! 
   mpfcs82_MPI.py - Multiprocessing with MPI
"""


# -------------------------------------------------------------------------------------------------
class MpfCS82MasterMPI(MpfxMasterMPI):
   
   """! Master calculator for MPI: distribute jobs to workers and collect/process their results """

   def __init__(self, args, comm):
      """! 
         Master class constructor 
         @param args list of command-line options
         @param comm communication channel   
      """

      # --- Master process
      MpfxMasterMPI.__init__(self, args, comm)


   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_run_output_dir(self):
      """! 
         Create run directory where output data will be stored
         @return time-stamped output directory tracing the run of the process 
      """

      return MpfxMasterMPI.create_run_output_dir(self)

   # -----------------------------------------------------------------------------------------------
   def create_job_processor(self):
      """
         Factory method for creating a Job Processor for managing the life-cycle of jobs. 
         @return instance of job processor object
      """

      return jobp.MpfCS82JobProcessor(self)

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
class MpfCS82WorkerMPI(MpfxWorkerMPI):      

   def __init__(self, args, comm, rank):
      """! Worker class constructor          
         @param args list of command-line options
         @param comm comunication channel
         @param rank rank of MPI process (> 0)
     """
      MpfxWorkerMPI.__init__(self, args, comm, rank)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point, either for the master or for a particular worker """

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()                 # rank: 0 for the main process, > 0 for workers

   if rank == 0:

      # --- Command-line arguments and options
      args = arg.MpfCS82Args()

      # --- Master process
      master = MpfCS82MasterMPI(args, comm)
      if master is None:
         sys.exit(2)
      else:
         try:
            master.run()
         except:
            print_error("The master process ended unexpectedly")
         finally:
            master.shutdown()
            sys.exit(0)

   else:

      # --- Command-line arguments and options
      args = arg.MpfCS82Args()

      # --- Worker process
      worker = MpfCS82WorkerMPI(args, comm, rank)  
      if worker is None:
         sys.exit(2)
      else:
         try:
            worker.run()
         except:
            print_error("The worker process {0} ended unexpectedly".format(worker.name))
         finally:
            sys.exit(0)

# -- EOF mpfcs82_MPI.py
