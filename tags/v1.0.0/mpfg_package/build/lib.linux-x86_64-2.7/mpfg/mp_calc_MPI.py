"""! 
   @package mpfg.mp_calc_MPI  multiprocessing infrastructure for calculation with the MPI interface 
   @author Marc Gentile
   @file mp_calc_MPI.py

   mp_calc_MPI.py: nultiprocessing infrastructure for calculation with the MPI interface
""" 

# -- Python imports
import os, sys
import Queue
from mpi4py import MPI  # MPI interface

# --- Module-specific imports
import mp_args as arg            # command-line arguments handling
from mp_calc import *            # base calculator class


# Request types
TAG_WORK      = 1     # tag to assign a job to a worker
TAG_RESULT    = 2     # tag to signal a job result is available
TAG_TERMINATE = 0    # tag to end a worker process

# -------------------------------------------------------------------------------------------------
class MasterMPI(Master):
   
   """! Master calculator for MPI: distribute jobs to workers and collect/process their results """

   def __init__(self, arg_options, comm):
      """! 
         Master class constructor 
         @param arg_options list of command-line options
         @param comm communication channel   
      """
      Master.__init__(self, arg_options)

      self._comm = comm                            # communication channel
      self._rank = 0                               # rank of master (i.e. main) process
      self._nb_workers = self.comm.Get_size()-1    # nb. of worker processes              
      self._input_queue  = Queue.Queue()           # input queue of jobs to process by workers

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def comm(self):
      """! @return the communicaton channel """
      return self._comm

   @property
   def rank(self):
      """! @return the rank of the master (0) """
      return self._rank

   @property
   def nb_workers(self):
      """! @return the number of worker processes managed by the master """
      return self._nb_workers

   @property
   def nb_workers(self):
      """! @return the number of worker processes managed by the master """
      return self._nb_workers
   @property
   def input_queue(self):
      """! @return the input queue of jobs to process by workers """
      return self._input_queue

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def submit_jobs(self):
      """! 
         Populate the job input queue containing the list of jobs to process to be dispached among 
         workers
         @note this queue is private to the master process and not shared by workers.
      """ 

      # Add jobs to process to the input queue (reverse order)
      for job in self.job_list:
         self.input_queue.put(job)

   # -----------------------------------------------------------------------------------------------
   def run(self):
      """! Run the jobs and process their results. Overrides Master.run(). """

      try:

         # --- Inform user of work to be done
         if self.logging_enabled:
            self.logger.log_info_p(
                   'Main Process - Starting processing data (MPI)...')
#            self.logger.log_info_p(
#                   'Main Process - Configuration file: {0}\n'.format(self._get_config_filepath()))

         # --- Run
         Master.run(self)

      except:
         self.helper.print_error("{0}".format(sys.exc_info()))
         #self.helper.print_error("{0}".format(sys.exc_info()[1]))

   # -----------------------------------------------------------------------------------------------
   def process_jobs(self):
      """! Have workers process the list of jobs and process the corresponding job results. """

      # --- Broadcast shared dictionary so that workers are able to access the corresponding data
      self.comm.bcast(self.shared_dico, root=0)

      # --- Only keep a minimum number of worker processes to optimally process the sets
      nb_jobs = len(self.job_list)
      nb_workers = self.nb_workers

      if nb_jobs < self.nb_workers:    # if more workers than necessary
         for iproc in range(nb_jobs+1, self.nb_workers+1):
            self.comm.send(None, dest=iproc, tag=TAG_TERMINATE)
            if self.logging_enabled():
               self.logger.log_info_p(
                     '{0} - Terminating redundant worker process {1} ...'.format(self.name, iproc))
         nb_workers = nb_jobs    # effective nb. of workers

      # --- Send data sets to available workers for processing
      for iproc in range(1, nb_workers+1):
         job = self._get_next_job(self.input_queue)   # next set to process

         if job is not None:
            if self.logging_enabled():
               self.logger.log_info_p(
                  '{0} - Dispatching job {1} to worker {2} ...'.format(self.name, job, iproc))

            self.comm.send(job, dest=iproc, tag=TAG_WORK)

      nb_received_results = 0  

      if self.logging_enabled():
         self.logger.log_info_p(
               '\n{0} - Waiting for worker processes to terminate ...\n'.format(self.name))

      # --- Process results and feed workers
      while not self.input_queue.empty():

         # Wait for any worker's result
         [rank, job_result] = self.comm.recv(source=MPI.ANY_SOURCE, tag=TAG_RESULT)
#         if self.logging_enabled():
#            self.logger.log_info_p('Main process - Received results from worker {0}'.format(rank))

         if job_result is None:
            if self.logging_enabled():
               self.logger.log_error_p(
                       '{0} - Worker {1} failed to deliver results'.format(self.name, rank))
            
         nb_received_results += 1
 
         # Give worker a new set to process
         new_job = self._get_next_job(self.input_queue)

         if new_job is not None:
            if self.logging_enabled():
               self.logger.log_info_p(
                  '{0} - Dispatching job {1} to worker {2} ...'.format(
                                                                  self.name, new_job, rank))

#         if __debug__ and self.logging_enabled():
#            self.logger.log_info_p(
#               'Main process - Master sending to worker {0} the new job {1}'.format(
#                                                                              rank, new_job.name))
         self.comm.send(new_job, dest=rank, tag=TAG_WORK)

         if job_result is not None:   
            # Store the Job result
            self.job_processor.record_job_result(job_result, self)

            # Process the results received from last worker
            self.process_job_result(job_result)

      # --- No more sets to process, wait for the last workers to finish their work...  
      nb_remaining_workers = nb_jobs - nb_received_results

      for iproc in range(0, nb_remaining_workers):

         [rank, job_result] = self.comm.recv(source=MPI.ANY_SOURCE, tag=TAG_RESULT)

         if job_result is None:
            if self.logging_enabled():
               self.logger.log_error_p(
                            '{0} - Worker {1} failed to deliver results'.format(self.name, rank))

         nb_received_results += 1

         if job_result is not None:   
            # Store the Job result
            self.job_processor.record_job_result(job_result, self)

            # Process the results from the last workers
            self.process_job_result(job_result)

         # Terminate this worker since no more work left
         self.comm.send(None, dest=rank, tag=TAG_TERMINATE)

      # --- Once all the jobs have been processed, give the opportunity to do something with
      #     the entire job result list.
      self.job_processor.all_jobs_processed(self)

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _get_next_job(self, queue):
      """! 
         Extract a job from a Queue until none. 
         @param queue private queue to hold the job list   
      """
      try:
         return queue.get_nowait()
      except (Queue.Empty):
         return None


# -------------------------------------------------------------------------------------------------
class WorkerMPI(Worker):
   
   """! Worker calculator: process jobs supplied by master and send back the results. """

   def __init__(self, arg_options, comm, rank):
      """! Worker class constructor          
         @param arg_options list of command-line options
         @param comm comunication channel
         @param rank rank of MPI process (> 0)
     """
      
      # --- Command-line options
      self._arg_options = arg_options
      Worker.__init__(self, self._arg_options)

      self._comm = comm                            # communication channel
      self._rank = rank                            # rank of worker process
      self._name = "Worker-{0:d}".format(rank)   # name of worker process

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def comm(self):
      """! @return the communication channel. """
      return self._comm

   @property
   def rank(self):
      """! @return the rank of the master (0). """
      return self._rank
      
   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def run(self):
      """! run the Worker process. """

      try:   

         # --- Access and set the shared data broadcastef by the master
         self.shared_dico = self.comm.bcast(None, root=0)

         # logging scope: none, per_worker or per job
         logging_scope = self.config.get_as_string("CREATE_WORKER_LOGGERS", "LOGGING").upper()

         # --- Logging
         if self.shared_dico["logging_enabled"] and logging_scope.upper() == "PER_WORKER":
            self.logger = self.create_logger()
            self.logger.log_info_p('{0} - Waiting for jobs to process...'.format(self.name))
            self.logger.flush()

         # --- Job processing
         self._job_processor = self.shared_dico["job_processor"]
         

         # --- Process jobs 
         while(True):

            try:   

               # --- Got a job to process from the master process
               job = self.comm.recv(source=0, tag=MPI.ANY_TAG)

               # --- Worker running...
               self._start_time = time.time()  # record start time

               # --- Check if worker requested to terminate
               if job is None:
                  if self.logging_enabled():
                     self.logger.log_info_p(
                        '{0} - "Terminating, processor time: {1:.3f} secs. ({2:.3f} hours)'.format(
                                          self.name, self.run_duration, self.run_duration/3600.0))
                  break;

               # --- Create log file if required
               if self.shared_dico["logging_enabled"] and logging_scope.upper() == "PER_JOB":
                  self.logger = self.create_logger(job.name)

               # --- pre-Process the job
               job_result = self.preprocess_job(job)

               # --- Process the job
               job_result = self.process_job(job)

               # --- pre-Process the job
               self.postprocess_job(job_result)

               # --- Send the job results to the master process
               self.comm.send([self.rank, job_result], dest=0, tag=TAG_RESULT)

               # --- Close log per job
               if self.shared_dico["logging_enabled"] and logging_scope.upper() == "PER_JOB":
                  hours = self.run_duration/3600.0
                  if hours > 1.0:
                     self.logger.log_info_p(
                        '{0} - Work complete - Processor time: {1:.3f} secs. ({2:.3f} hours)'.format(
                                    self.name, self.run_duration, hours))
                  else:
                     self.logger.log_info_p(
                        '{0} - Work complete - Processor time: {1:.3f} secs. ({2:.3f} minutes)'.format(
                                    self.name, self.run_duration, hours * 60))
                  self._close_logger(self.logger)        

            except Exception as detail:

               self.helper.print_error(
                   "{0} - Exception thrown while processing a job - {1}".format(
                                                                        self.name, detail))     
 
               # --- Notify the master of the exception
               self.comm.send([self.rank, None], dest=0, tag=TAG_RESULT)

         # --- end while

         # --- Inform user of work completion
         if self.logging_enabled() and logging_scope.upper() == "PER_WORKER":
            hours = self.run_duration/3600.0
            if hours > 1.0:
               self.logger.log_info_p(
                  '{0} - Work complete - Processor time: {1:.3f} secs. ({2:.3f} hours)'.format(
                                    self.name, self.run_duration, hours))
            else:
               self.logger.log_info_p(
                  '{0} - Work complete - Processor time: {1:.3f} secs. ({2:.3f} minutes)'.format(
                                    self.name, self.run_duration, hours * 60))
      except Exception as detail:
         self.helper.print_error(
               "{0} A fatal exception occurred: {1} - Shutting down Master".format(
                                                                     self.name, detail))

      finally:
         self._close_logger(self.logger)        


# -- EOF mp_calc_MPI.py_
