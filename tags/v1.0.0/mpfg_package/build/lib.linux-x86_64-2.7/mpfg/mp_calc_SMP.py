"""! 
   @package mpfg.mp_calc_SMP multiprocessing infrastructure for calculation on a SMP architecture  
   @author Marc Gentile
   @file mp_calc_SMP.py

   mp_calc_SMP.py: multiprocessing infrastructure for calculation on a SMP architecture
""" 

# -- Python imports
import os, sys
import Queue
import multiprocessing as multi
import time #FSC MODIF

# --- Module-specific imports
from mp_calc import *            # base calculator class
from mp_helper import *          # miscellaneous utility functions

# -------------------------------------------------------------------------------------------------
class MasterSMP(Master):
   
   """! Master calculator for SMP: distribute jobs to workers and collect/process their results """

   def __init__(self, arg_options):
      """!
         Master class constructor 
         @param arg_options list of command-line options
      """

      Master.__init__(self, arg_options)

      self._input_queue  = multi.Queue()     # input queue of jobs to process by workers
      self._output_queue = multi.Queue()     # output queue of jobs having results from workers

      self._start_event = multi.Event()      # event to have all the worker processes start together

      self._workers = {} # dictionary of started worker processes ([process index, process object])

      # Nb. of worker among which the jobs will be distributed
      self._nb_workers = self.config.get_as_int('WORKER_PROCESSES_COUNT', 'MULTIPROC_SMP') 
      if self._nb_workers == 0:   
         self._nb_workers = multi.cpu_count()-1  # auto-detected

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def workers(self):
      """! @return a dictionary of workers: [process index, process object] managed by the Master """
      return self._workers

   @property
   def nb_workers(self):
      """! @return the number of worker processes managed by the Master """
      return self._nb_workers

   @property
   def input_queue(self):
      """! @return the input queue of jobs to process by Worker objects """
      return self._input_queue

   @property
   def output_queue(self):
      """! @return the output queue of jobs having results from Worker objects """
      return self._output_queue

   @property
   def start_event(self):
      """! Event to have all the processes start together """
      return self._start_event

   # --- Setters

   @workers.setter 
   def workers(self, workers):
      """! 
         Set the list of Worker object managed by the Master 
         @param workers list of Worker object instances
      """
      self._workers = workers 

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def run(self):
      """! 
         Run the jobs and process their results.
         @note overrides Master.run() """

      try:

         # --- Inform user of work to be done
         if self.logging_enabled():
            self.logger.log_info_p(
                   'Main Process - Starting processing data (SMP)...')
            self.logger.log_info_p(
                   'Main Process - Configuration file: {0}\n'.format(self._get_config_filepath()))

         # --- Create worker processes
	 #self.start_event.clear() #TO BE CONFIRMED
         self.create_workers()

         # --- Run workers
         #FSC Added: problem, in code, jobs are populated in subsequent call Master.run(self)
	 #i.e. AFTER start signal launched in run_workers(). Cause sometimes problems. MOVED TO end of submit_jobs 
         #self.run_workers()  

         Master.run(self)

      except Exception as detail:
         Helper.print_error("{0}".format(detail))

   # -----------------------------------------------------------------------------------------------
   def create_workers(self):

      """! Create the pool of workers that will process the jobs """

      # Create the worker processes
      for iproc in range(0, self.nb_workers):

         # Arguments to pass to the worker process  
         worker_args = (self.input_queue, self.output_queue, self.start_event, self.shared_dico)
         
         # Create & Start worker process
         worker = self.create_worker(self.arg_options)
         self.workers[iproc] = multi.Process(target=worker.run, args=worker_args)

         if self.logging_enabled():
            self.logger.log_info_p('Main Process - Starting worker {0} ...'.format(
                                                                        self.workers[iproc].name))

   # -----------------------------------------------------------------------------------------------
   def create_worker(self, arg_options):
      """! 
         Create a worker process, method is to be overridden pas subclasses
         @param arg_options list of command-line options
      """
      return WorkerSMP(arg_options)
      
   # -----------------------------------------------------------------------------------------------
   def submit_jobs(self):
      """! Submit the jobs for processing """
      
      # Add jobs to process to the input queue (reverse order)
      for job in self.job_list:
         self.input_queue.put(job)
      self.run_workers()

   # -----------------------------------------------------------------------------------------------
   def run_workers(self):
      """! Start worker processes """

      # Start the worker processes
      for iproc in range(0, self.nb_workers):
         self.workers[iproc].start()

   # -----------------------------------------------------------------------------------------------
   def process_jobs(self):
      """! Have workers process the list of jobs and process the corresponding job results """

      # Read results of processed jobs from the output queue
      nb_jobs_processed = 0

      # Signal the waiting processes that they can start running
      self.start_event.set()

      if self.logging_enabled():
         self.logger.log_info_p('Main Process - Waiting for worker processes to terminate ...\n')

      while nb_jobs_processed < len(self.job_list):

         # --- Wait until a result gets available
         job_result = self.output_queue.get()    

         if job_result is None:
            if self.logging_enabled():
               self.logger.log_error_p(
                            '{0} - A worker failed to deliver results'.format(self.name))
         else:

            # --- Store Job result
            self.job_processor.record_job_result(job_result, self)

            # --- Process job result
            self.process_job_result(job_result)

         nb_jobs_processed += 1

      # --- Once all the jobs have been processed, give the opportunity to do something with
      #     the entire job result list.
      self.job_processor.all_jobs_processed(self)

   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed.
         @param master Master process instance 
      """
      pass


# -------------------------------------------------------------------------------------------------
class WorkerSMP(Worker):
   
   """! Worker calculator: process jobs supplied by master and send back the results """

   def __init__(self, arg_options):
      """
         Worker class constructor 
         @param arg_options list of command-line options
      """
      Worker.__init__(self, arg_options)

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def run(self, input_queue, output_queue, start_event, shared_dico):
      """! 
         Run Worker process 
         @param input_queue queue for input containing the Job objects to process by the worker
         @param output_queue queue for output containing the JobResult objects of processed jobs
         @param start_event synchronization event between the Master and the Worker objects
         @param shared_dico dictionary for fast communication between Master and Worker objects
      """

      try:

         self.name = multi.current_process().name  # set worker process name
         self.shared_dico = shared_dico            # common dictionary between master and workers

         # Wait until the queue gets populated by the main thread
         start_event.wait()
         # logging scope: none, per_worker or per job
         logging_scope = self.config.get_as_string("CREATE_WORKER_LOGGERS", "LOGGING").upper()

         # --- Logging per worker
         if self.shared_dico["logging_enabled"] and logging_scope.upper() == "PER_WORKER":
            self.logger = self.create_logger()

         # --- Job processing
         self._job_processor = self.shared_dico["job_processor"]

         # --- Process jobs 
         #print(self.name," EMPTY QUEUE ?",input_queue.empty(),input_queue.qsize())
	 time.sleep(20)
         # Worker processing loop
         #while not input_queue.empty(): 
         while True: 
            try:

               # --- Retrieve available job
               job = input_queue.get(True, 20)
               print(self.name," PROCESS ",job.name)

               # Worker running...
               self._start_time = time.time()  # record actual start time

               # --- Logging per job
               if self.shared_dico["logging_enabled"] and logging_scope.upper() == "PER_JOB":
                  self.logger = self.create_logger(job.name)

               # --- pre-Process the job
               job_result = self.preprocess_job(job)

               # --- Process the job
               job_result = self.process_job(job)

               # --- pre-Process the job
               self.postprocess_job(job_result)

               # --- Send result to master
               output_queue.put(job_result)

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

	    except Queue.Empty as detail:
	       print(self.name," queue empty EXCEPTION")

            except Exception as detail:
               print(self.name," EXCEPTION")
               if not input_queue.empty():
                  self.helper.print_error(
                      "{0} - Exception thrown while processing a job - {1}".format(
                                                                        self.name, detail))     
                  # --- Notify the master of the exception
                  output_queue.put(None)
	    if input_queue.empty() :
	        break
         print(self.name," EMPTY QUEUE")

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

      except(Queue.Empty):
         pass     # queue empty, ignoring

      except Exception as detail:
         self.helper.print_error(
               "{0} A fatal exception occurred: {1} - Shutting down Master".format(
                                                                     self.name, detail))
      finally:
         self._close_logger(self.logger)        


# -- EOF mp_calc_SMP.py_
