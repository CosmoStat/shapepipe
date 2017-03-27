"""! 
   @package mpfg.mp_calc multiprocessing infrastructure for calculation - Base module
   @author Marc Gentile
   @file mp_calc.py
   Multiprocessing infrastructure for calculation - Base module
""" 

# -- Python imports
import os, sys
import time
from operator import itemgetter, attrgetter

# -- External imports
from slogger import FileLogger   # logging

# --- Module-specific imports
from mp_job import *                # job processing
from mp_helper import *             # utility functions

# -------------------------------------------------------------------------------------------------
class Master(object):
   
   """! Master calculator: distribute jobs to workers and collect/process their results """

   def __init__(self, arg_options):
      """
         ! Master class constructor 
         @param arg_options list of command-line options
      """

      self._start_time = time.time()  # start time

      self._helper = Helper()

      # --- Command-line options
      self._arg_options = arg_options 

      # --- External configuration
      self._default_config_dir = self.get_default_config_dir()
      self._default_config_filename = self.get_default_config_filename()
      self._config = self.helper.read_config(self.arg_options.options)

      if self._config is None:
         sys.exit(2)  

      self._name = "Main Process"  # name of the master process

      # --- Shared dictionary between master and workers (NOT shared memory)
      self._shared_dico = {}

      # --- Setup directory tree
      self._run_output_dir = self.create_run_output_dir() 
      if not self._create_directory_tree():
         sys.exit(2)  
      self._record_dynamic_paths() # record dynamic paths of tree structure

      # --- Setup logging
      self._logger = self.create_logger()
      self.shared_dico["logging_enabled"] = not self._logger is None
      if self.shared_dico["logging_enabled"]:
         # Dump current configuration data
         self.helper.dump_config(self._config, self._logger) 
         self._logger.flush()

      # --- Instantiate a Job Processor
      self._job_processor = self.create_job_processor() 
      self.shared_dico["job_processor"] = self._job_processor  

   def __str__(self):
      """!
          Formatting for display 
          @return string representation of the object  
      """      
      return self.name

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def arg_options(self):
      """! @return command-line arguments and options """
      return self._arg_options

   @property
   def config(self):
      """! @return configuration data in the form of a SConfig object """
      return self._config

   @property
   def default_config_dir(self):
      """! @return default configuration directory (e.g. './config') """
      return self._default_config_dir

   @property
   def default_config_filename(self):
      """! @return default configuration filename (e.g. 'mpfg.cfg') """
      return self._default_config_filename

   @property
   def config_dir(self):
      """! @return actual configuration directory """
      return self._arg_options.options["-d"]

   @property
   def config_filename(self):
      """! @return actual configuration filename """
      return self._arg_options.options["-c"]

   @property
   def logger(self):
      """! @return the associated SLogger file logger object """
      return self._logger

   @property
   def name(self):
      """! @return the name of the worker """
      return self._name

   @property
   def base_input_dir(self):
      """! @return the base input directory """
      return self._base_input_dir

   @property
   def base_output_dir(self):
      """! @return the base output directory """
      return self._base_output_dir

   @property
   def run_output_dir(self):
      """! @return the directory were output is stored """
      return self._run_output_dir

   @property
   def log_output_dir(self):
      """! @return the directory were logs are created are stored """
      return self._log_output_dir

   @property
   def plot_output_dir(self):
      """! @return the directory were plots are stored """
      return self._plot_output_dir

   @property
   def map_output_dir(self):
      """! @return the directory were maps are stored """
      return self._map_output_dir

   @property
   def stat_output_dir(self):
      """! @return the directory were stattistics are recorded """
      return self._stat_output_dir

   @property
   def result_output_dir(self):
      """! @return the directory were results are stored """
      return self._result_output_dir

   @property
   def error_output_dir(self):
      """! @return the directory were errors are recorded """
      return self._error_output_dir

   @property
   def job_processor(self):
      """! @return the job processor """
      return self._job_processor

   @property
   def job_list(self):
      """! @return the jobs to process as a list """
      return self._job_processor.job_list

   @property
   def shared_dico(self):
      """! @return shared dictionary between master and workers """
      return self._shared_dico

   @property
   def run_duration(self):
      """! @return the total run duration of the master in seconds """
      return time.time() - self._start_time

   # --- Setters

   @name.setter 
   def name(self, name):
      """! 
         set the name of the Master process 
         @param name name of Master process 
      """
      self._name = name 

   @job_processor.setter 
   def job_processor(self, job_processor):
      """!
         set the job processor
         @param job_processor instance of a object of class JobProcessor 
      """
      self._job_processor = job_processor
      self.shared_dico["job_processor"] = self._job_processor 

   @logger.setter 
   def logger(self, logger):
      """! 
         Set the file logger 
         @param logger a FileLogger instance
      """
      self._logger = logger 

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def logging_enabled(self):
      """!
         Tell of a FileLogger has been created
         @retval True if created
         @retval False otherwise   
      """
      try:
         return self.logger is not None
      except:
         return False

   # -----------------------------------------------------------------------------------------------
   def run(self):
      """! Run the jobs and process their results """

      try:

         # --- Locate objects (e.g. images or catalogs) to process
         self.job_processor.create_jobs(self)

         if self.job_processor.job_count == 0:
            if self.logging_enabled():
               self.logger.log_info_p('Main Process - No job found to process in specified dataset.')
         else:
            if self.logging_enabled():
               if self.job_processor.job_count < 50:
                  self.logger.log_info_p('Main Process - Found {0} jobs to process: {1}'.format( 
                                         self.job_processor.job_count, 
                                         ["{0}".format(j) for j in self.job_processor.job_list]))
               else:
                  self.logger.log_info_p('Main Process - Found {0} jobs to process.'.format( 
                                         self.job_processor.job_count))
               self.logger.flush()
               
         # --- Submit jobs
         self.submit_jobs()

         if self.logging_enabled():
            self.logger.flush()

         # --- Process the jobs
         self.process_jobs()

         # --- Inform user of work completion
         if self.logging_enabled():
            elapsed_time = self.run_duration

            hours = elapsed_time/3600.0
            if hours > 1.0:
               self.logger.log_info_p(\
                'Main Process - Terminating, processor time: {0:.3f} secs. ({1:.3f} hours)'.format(
                                                      elapsed_time, hours))
            else:
               self.logger.log_info_p(\
               'Main Process - Terminating, processor time: {0:.3f} secs. ({1:.3f} minutes)'.format(
                                                      elapsed_time, hours * 60))

            self.logger.log_info_p('Main Process - Output directory: {0}'.format(
                                                               self.shared_dico["run_output_dir"])) 
            if self.logging_enabled():
               self.logger.flush()

      except Exception as detail:
         self.helper.print_error("{0}".format(detail))
         #self.helper.print_error("{0}".format(sys.exc_info()[1]))

   # -----------------------------------------------------------------------------------------------
   def shutdown(self):
      """! Shutdown master process: close files, stop processes, etc. """

      if self.logging_enabled():
         self.logger.log_info_p("Main Process - Stopping...")

         # --- Logging   
         self._close_logger(self.logger)

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result):
      """! 
         Process the result associated with a processed job.

         @param job_result object of class JobResult
         @note this method is to be overridden by subclasses to exploit the
               content of the job result object, like making plots of 
               creating catalogs
      """

      self.job_processor.process_job_result(job_result, self)

   # -----------------------------------------------------------------------------------------------
   def get_default_config_dir(self):
      """! @return the default configuration file directory (to be overridden by subclasses) """
      return self.arg_options.default_config_dir

   # -----------------------------------------------------------------------------------------------
   def get_default_config_filename(self):
      """! @return the default configuration filename (to be overridden by subclasses) """
      return self.arg_options.default_config_filename

   # -----------------------------------------------------------------------------------------------
   def create_logger(self):
      """! 
          Create master File Logger  
          @return the file logger if success, None otherwise
       """

      if self.arg_options.options['-q']:   
         # Force configuration to quiet mode
         return None
      
      if self.helper.is_logging(self.config):
         # Create file logger
         log_time = time.strftime("%d.%m.%y_%H.%M.%S", time.localtime())
         log_file_name = 'run_{0}.log'.format(log_time)
         log_output_dir = self.shared_dico["log_output_dir"]
         logger = FileLogger(logger='', directory=log_output_dir, filename=log_file_name)       
         if logger is None:
            self.helper.print_error("Could not create file logger {0}".format(log_file_name))  
         return logger
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def create_run_output_dir(self):
      """! 
         Create run directory where output data will be stored
         @return time-stamped output directory tracing the run of the process 
      """

      output_date = time.strftime("%d.%m.%y_%H.%M.%S", time.localtime())
      run_output_dir = ('run_{0}'.format(output_date))
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
         @return instance of JobProcessor
      """

      return JobProcessor(self) 



   # ~~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _create_directory_tree(self):
      """! 
         setup tree structure (input, output, ...) 
         @retval True if success
         @retval False otherwise
      """

      # --- Input
      self._base_input_dir = self.config.get_as_string("BASE_INPUT_DIR", "DIR.INPUT")
      if not os.path.isdir(self._base_input_dir):
         self.helper.print_info("Creating directory {0}...".format(self._base_input_dir))
         self.helper.make_dir(self._base_input_dir)

      # --- Output

      # Base output dir
      self._base_output_dir = self.config.get_as_string("BASE_OUTPUT_DIR", "DIR.OUTPUT")
      if not os.path.isdir(self._base_output_dir):
         self.helper.print_info("Creating directory {0}...".format(self._base_output_dir))
         self.helper.make_dir(self._base_output_dir)

      # Create sub-directories for logs, plots, results
      self._log_output_dir    = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_LOG_DIR', "DIR.OUTPUT"))
      self._plot_output_dir   = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_PLOT_DIR', "DIR.OUTPUT"))
      self._map_output_dir    = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_MAP_DIR', "DIR.OUTPUT"))
      self._stat_output_dir   = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_STAT_DIR', "DIR.OUTPUT"))
      self._result_output_dir = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_RESULT_DIR', "DIR.OUTPUT"))
      self._error_output_dir  = os.path.join(self.run_output_dir, 
                                       self.config.get_as_string('OUTPUT_ERROR_DIR', "DIR.OUTPUT"))

      self.helper.make_dir(self.log_output_dir)
      self.helper.make_dir(self.plot_output_dir)
      self.helper.make_dir(self.map_output_dir)
      self.helper.make_dir(self.stat_output_dir)
      self.helper.make_dir(self.result_output_dir)
      self.helper.make_dir(self.error_output_dir)

      return True

   # -----------------------------------------------------------------------------------------------
   def _record_dynamic_paths(self):
      """! Record dynamic path for IPC between master and workers """

      self.shared_dico["base_input_dir"]    = self.base_input_dir
      self.shared_dico["base_output_dir"]   = self.base_output_dir
      self.shared_dico["run_output_dir"]    = self.run_output_dir
      self.shared_dico["log_output_dir"]    = self.log_output_dir
      self.shared_dico["plot_output_dir"]   = self.plot_output_dir
      self.shared_dico["map_output_dir"]    = self.map_output_dir
      self.shared_dico["stat_output_dir"]   = self.stat_output_dir
      self.shared_dico["result_output_dir"] = self.result_output_dir
      self.shared_dico["error_output_dir"]  = self.error_output_dir         

   # -----------------------------------------------------------------------------------------------
   def _close_logger(self, logger):
      """!
         Close the file logger associated with the Master process
         @param logger instance of a fileLogger object
      """
      if not logger is None:
         logger.close()
         self.logger = None     

   # -----------------------------------------------------------------------------------------------
   def _get_config_filepath(self):
      """! Return the configuration file full path """
      return os.path.join(self.arg_options.options["-d"], self.arg_options.options["-c"])


# -------------------------------------------------------------------------------------------------
class Worker(object):      

   """! Worker calculator: process jobs supplied by master and send back the results """

   def __init__(self, arg_options):
      """! 
         Worker constructor
         @param arg_options list of command-line options
      """

      self._start_time = time.time()  # start time

      # --- Helper utility class
      self._helper = Helper()

      # --- Command-line options
      self._arg_options = arg_options     # command-line arguments and options      

      # --- External configuration
      self._config = self.helper.read_config(self.arg_options.options)
      if self._config is None:
         sys.exit(2)    

      self._default_config_dir = self.get_default_config_dir()

      # --- Common dictionary between master and workers
      self._shared_dico = None

      # --- Logger
      self._logger = None


   def __str__(self):
      """!
          Formatting for display 
          @return string representation of the object  
      """
      return self.name

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the Helper class """
      return self._helper

   @property
   def arg_options(self):
      """! @return command-line arguments and options """
      return self._arg_options

   @property
   def config(self):
      """! @return configuration data in the form of a SConfig object """
      return self._config

   @property
   def default_config_dir(self):
      """! @return default configuratio directory (e.g. './config') """
      return self._default_config_dir

   @property
   def logger(self):
      """! @return the associated SLogger file logger object """
      return self._logger

   @property
   def name(self):
      """! @return the name of the worker """
      return self._name

   @property
   def shared_dico(self):
      """! @return the common dictionary between master and workers """
      return self._shared_dico

   @property
   def job_processor(self):
      """! @return the job processor """
      return self._job_processor

   @property
   def run_duration(self):
      """! @return the total run duration of the master in seconds """
      return time.time() - self._start_time

   @property
   def base_input_dir(self):
      """! @return the base input directory were the data is taken from """
      return self.shared_dico["base_input_dir"]

   @property
   def base_output_dir(self):
      """! @return the base output directory were output is stored """
      return self.shared_dico["base_output_dir"]

   @property
   def run_output_dir(self):
      """! @return the directory were output is stored """
      return self.shared_dico["run_output_dir"]

   @property
   def log_output_dir(self):
      """! @return the directory were logs are created are stored """
      return self.shared_dico["log_output_dir"]

   @property
   def plot_output_dir(self):
      """! @return the directory were plots are stored """
      return self.shared_dico["plot_output_dir"]

   @property
   def map_output_dir(self):
      """! @return the directory were maps are stored """
      return self.shared_dico["map_output_dir"]

   @property
   def stat_output_dir(self):
      """! @return the directory were stattistics are recorded """
      return self.shared_dico["stat_output_dir"]

   @property
   def result_output_dir(self):
      """! @return the directory were results are stored """
      return self.shared_dico["result_output_dir"]

   @property
   def error_output_dir(self):
      """! @return the directory were errors are recorded """
      return self.shared_dico["error_output_dir"]

   # --- Setters

   @logger.setter 
   def logger(self, logger):
      """! 
         Set the file logger 
         @param logger a FileLogger instance
      """
      self._logger = logger

   @name.setter 
   def name(self, name):
      """! 
         Set the name of the worker 
         @param name worker name
      """
      self._name = name

   @shared_dico.setter 
   def shared_dico(self, shared_dico):
      """! 
         Set the common dictionary between master and workers 
         @param shared_dico instance of a dictionary containing shared data
                between the Master and Workers processes 
      """
      self._shared_dico = shared_dico

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_default_config_dir(self):
      """! @return the default configuration file directory (to be overridden by subclasses) """
      return self.arg_options.default_config_dir

   # -----------------------------------------------------------------------------------------------
   def get_default_config_filename(self):
      """! @return the default configuration file name (to be overridden by subclasses) """
      return self.arg_options.default_config_filename

   # -----------------------------------------------------------------------------------------------
   def logging_enabled(self):
      """!
         Tell of a FileLogger has been created
         @retval True if created
         @retval False otherwise   
      """
      try:
         return self.logger is not None
      except:
         return False

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job):
      """! 
         Pre-Process a job, delegating the actual processing to the JobProcessor 
         @param job an object of class Job to pre-process
         @see process_job, postprocess_job
      """
      return self.job_processor.preprocess_job(job, self)

   def process_job(self, job):
      """! 
         Process a job and return the corresponding results. Delegates
         the processing to the JobProcessor 
         @param job object of class Job to process
         @see preprocess_job, postprocess_job
      """
      
      return self.job_processor.process_job(job, self)

   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result):
      """! 
         Post-Process a job, delegating the actual processing to the JobProcessor
         @param job_result object of class JobResult with processed data
         @see process_job, preprocess_job
      """
      return self.job_processor.postprocess_job(job_result, self)

   # -----------------------------------------------------------------------------------------------
   def create_logger(self, extra=None):
      """! Create worker File Logger  """

      if self.config.get_as_string("CREATE_WORKER_LOGGERS", "LOGGING").upper() != "NONE":
         log_time = time.strftime("%d.%m.%y_%H.%M.%S", time.localtime())
         if extra is None:   
            log_file_name = 'run_{0}_{1}.log'.format(self.name, log_time)
         else:
            log_file_name = 'run_{0}_{1}.log'.format(extra, log_time)

         log_output_dir = self.shared_dico['log_output_dir']
         logger = FileLogger(logger='', directory=log_output_dir, filename=log_file_name)       
         if logger is None:
            self.helper.print_error("Could not create worker file logger {0}".format(log_file_name))  
         return logger
      else:
         return None      

   # ~~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _close_logger(self, logger):
      if not logger is None:
         logger.close()
         self.logger = None     


# -- EOF mp_calc.py
