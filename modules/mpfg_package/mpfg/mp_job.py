"""!
   @package mpfg.mp_job job management for SMP or MPI interface
   @author Marc Gentile
   @file mp_job.py
   Job management for SMP or MPI interface
"""

# -- Python imports
import os, sys
import time
from operator import itemgetter, attrgetter

# -- External imports

# --- Module-specific imports
from mp_data import *               # data access
from mp_helper import *             # utility functions

# -------------------------------------------------------------------------------------------------
class JobProcessor(object):

   """! Submit jobs and process associated job results. """

   def __init__(self, master):
      """! Job Processor default constructor """

      self._job_list = []           # list of jobs (Job objects) to process
      self._job_result_list = []    # list of JobResult objects for each processed job
      self._helper = Helper()       # helper utility functions

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper singleton. """
      return self._helper

   @property
   def job_list(self):
      """! @return the jobs to process as a list. """
      return self._job_list

   @property
   def job_result_list(self):
      """! @return the results of the processed jobs as a list. """
      return self._job_result_list

   @property
   def job_count(self):
      """! @return the number of jobs to process. """
      return len(self._job_list)


   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """! Invoked by worker to perform some optional pre-processing on the job

         @param job an object of class Job to pre-process
         @param worker instance of the worker process

         @note this method is called <b>before</b> process_job_result()
         @note  this method is to be overridden by subclasses to perform
                any processing required <b>before</b> process_job() is called
         @see process_job, postprocess_job, Job
      """
      pass

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """!
         Process a job and return the corresponding results to the Master (here sample data).
         Processing the results may consist e.g. in storing files to the output result directory
         tree, draw plots, compute statistics...

         @param job an object of class Job to process
         @param worker instance of the worker process
         @return an object of class JobResult containing the data of the processed job
         @note This method should be overridden in a subclass of JobProcessor
         @see preprocess_job, postprocess_job, Job, JobResult
      """

      if worker.logging_enabled():
         worker.logger.log_info_p("{0} - Processing job {1}...".format(worker.name, job))

      sample_data = {}  # sample result data

      return JobResult(worker, job, sample_data)  # create a JobResult containing the result data

   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """!
         Invoked by worker to perform some optional post-processing on the job.
         @param job_result object of class JobResult with processed data
         @param worker instance of the worker process

         @note this method is called <b>after</b> process_job_result()
         @note this method is to be overridden by subclasses to perform
               any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult
      """
      pass

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """!
         Process the result associated with a processed job.

         @param job_result object of class JobResult with processed data
         @param master Master object instance
         @note This method should be overridden in a subclass of JobProcessor
         @see Job, JobResult
      """
      if master.logging_enabled():
         master.logger.log_info_p("Main process - Processing result from job {0}...".format(
                                                                             job_result.job))

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """!
         Locate all objects to process and create the corresponding jobs.
         @param master Master object instance
         @return the total number of created jobs
      """

      # --- Test Dataset
      dataset_name = master.config.get_as_string("NAME", "PRIMARY_DATASET")
      dataset_base_dir = master.config.get_as_string("BASE_DIR", "PRIMARY_DATASET")

      dataset = self.create_dataset(master, dataset_name, dataset_base_dir)

      # --- Search for matching files to process
      dataset_pattern_list = master.config.get_as_list("FILE_PATTERNS", "PRIMARY_DATASET")

      self._job_list = [FileJob(master, filename, dataset) for filename in dataset.query(
                                                      master, dataset_pattern_list, recurse=True)]

      if master.logging_enabled():
         master.logger.log_info_p(
           "Main process - Looking for files to process in dataset: {0} with query: {1}...".format(
                           dataset.name, dataset_pattern_list))

      return len(self._job_list)   # return number of jobs

   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, name, base_dir):
      """!
         Create a Dataset object that represent the data source for images and catalogs.
         @param master Master object
         @param name dataset name
         @param base_dir dataset base directory
         @return the Dataset instance
      """
      return Dataset(master, name, base_dir)

   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """!
         Invoked by the Master to add a job result to the list of job results.
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job
         @see Job, JobResult
      """

      self.job_result_list.append(job_result)

   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """!
         This method is called once all the jobs have been processed.
         @param master instance of the Master
         @note Job results can be obtained with by calling job_result_list
         @see Job, JobResult
      """
      pass


# -------------------------------------------------------------------------------------------------
class Job(object):

   """!
       Represents a job to process by a calculator. Can be extended to include application-specific
       information.
   """

   def __init__(self, master, name, dataset):
      """!
         Class constructor
         @param master Master object instance
         @param name the name of the Job object to create
         @param dataset the data source used to create the job
      """

      self._name = name          # job name (assumed unique)
      self._dataset = dataset

   def __str__(self):
      """!
          Formatting for display
          @return string representation of the object
      """
      return "<{0}>".format(self.name)


   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def name(self):
      """! @return the name of the job. """
      return self._name

   @property
   def dataset(self):
      """! @return data source used to create the job. """
      return self._dataset

   # --- Setters

   @name.setter
   def name(self, name):
      """!
         Set the job name
         @param name name of he job
      """
      self._name = name


# -------------------------------------------------------------------------------------------------
class FileJob(Job):

   """!
      Represents a job to process in the form of a file path. Can be extended to include
      application-specific information.
      @see Job
   """

   def __init__(self, master, filepath, dataset):

      """!
         Class constructor
         @param master Master object instance
         @param filepath filepath of the job ro process
         @param dataset the data source used to create the job
      """

      # Directory + file name
      self._filepath = filepath
      self._directory, self._filename = os.path.split(filepath)

      Job.__init__(self, master, self._filename, dataset)

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def filepath(self):
      """! @return associated full file path. """
      return self._filepath

   @property
   def filename(self):
      """! @return associated file name. """
      return self._filename

   @property
   def directory(self):
      """! @return associated directory. """
      return self._directory

   # --- Setters

   @filepath.setter
   def filepath(self, filepath):
      """!
         Set the file path associated with the Job.
         @param filepath file path associated with the Job
      """
      self._filepath = filepath


# -------------------------------------------------------------------------------------------------
class JobResult(object):

   """!
       Represents the result of a Job object processed by a calculator. Can be extended to include
       application-specific information.
   """

   def __init__(self, worker, job, result):
      """!
         Class constructor
         param job instance of a Job class
         param job instance of a JobResult class with the results from the processed job
      """

      self._job = job
      self._result = result

   def __str__(self):
      """!
          Formatting for display
          @return string representation of the object
      """
      return "Result for <{0}>: {1}".format(self.job.name, self.result)

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def job(self):
      """! @return associated job object. """
      return self._job

   @property
   def result(self):
      """! @return the result associated with the corresponding job. """
      return self._result


   # --- Setters

   @result.setter
   def result(self, result):
      """!
         Set the JobResult object containing the results from the processed job
         @param result JobResult object with the result associated with the corresponding job
      """
      self._result = result


# -- EOF mp_job.py
