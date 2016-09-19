#!/usr/bin/env python

# -- Python imports
import sys
from mpi4py import MPI  # MPI interface
from operator import itemgetter, attrgetter

# -- External Modules
from mpfg.mp_calc_MPI import *       # mpf: multiprocessing framework MPI mgt

# -- Module-specific imports
import mpfx_args as arg             # command-line arguments and options
import mpfx_job as jobp             # job processing

"""! 
   mpfx_MPI.py - Multiprocessing with MPI
"""


# -------------------------------------------------------------------------------------------------
class MpfxMasterMPI(MasterMPI):
   
   """! Master calculator for MPI: distribute jobs to workers and collect/process their results """

   def __init__(self, args, comm):
      """! 
         Master class constructor 
         @param args list of command-line options
         @param comm communication channel   
      """

      # --- Master process
      MasterMPI.__init__(self, args, comm)

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

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
class MpfxWorkerMPI(WorkerMPI):      

   def __init__(self, args, comm, rank):
      """! Worker class constructor          
         @param args list of command-line options
         @param comm comunication channel
         @param rank rank of MPI process (> 0)
     """
      WorkerMPI.__init__(self, args, comm, rank)


# ~~~~~~~~~~~~~~~~~
# Main entry point
# ~~~~~~~~~~~~~~~~~

if __name__ == '__main__':
   """! Main entry point, either for the master or for a particular worker """

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()                 # rank: 0 for the main process, > 0 for workers

   if rank == 0:

      # --- Command-line arguments and options
      args = arg.MpfxArgs()

      # --- Master process
      master = MpfxMasterMPI(args, comm)
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
      args = arg.MpfxArgs()

      # --- Worker process
      worker = MpfxWorkerMPI(args, comm, rank)  
      if worker is None:
         sys.exit(2)
      else:
         try:
            worker.run()
         except:
            print_error("The worker process {0} ended unexpectedly".format(worker.name))
         finally:
            sys.exit(0)

# -- EOF mpfx_MPI.py
