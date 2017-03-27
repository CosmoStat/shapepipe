"""! 
   @package mks.mks_job Job Management
   @author Marc Gentile
   @file mks_job.py
   Job Management
""" 

# -- Python imports
import os, sys
import time
import numpy
from operator import itemgetter, attrgetter

# -- External imports
from mpfx.mpfx_job import *            # base job processing

# --- Module-specific imports
from mks_sim import MksSimulator       # simulation management
from mks_help import *                 # helper utility functions


# -------------------------------------------------------------------------------------------------
class MksJobProcessor(MpfxJobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpf.JobProcessor.   
   """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, master):

      """! 
         Job Processor constructor
         @param master master process instance
      """
      
      MpfxJobProcessor.__init__(self, master)            # job processor

      self._simulator = MksSimulator(master)             # simulation management      
      self._helper    = master.helper                    # helper utility functions

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def simulator(self):
      """! @return the MksSimulator instance. """
      return self._simulator

   @property
   def helper(self):
      """! @return the MksHelper instance. """
      return self._helper

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def match_file(self, master, directory, filename, pattern):
      """! 
          File matching predicate method. Must return True in case of matching, False otherwise.   
          May be overriden by subclasses to set additional criteria. 
          @param master master object instance
          @param directory directory of filename
          @param filename file name
          @param pattern Unix-like file pattern 
          @return True of a match is found, False otherwise
      """

      return glob.fnmatch.fnmatch(filename, pattern)

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """! 
         Locate all objects to process and create the corresponding jobs. 
         @param master Master object instance
         @return the list of created jobs 
      """

      # --- Determine how many jobs to create      
      nb_fields =  master.config.get_as_int("NB_FIELDS", "SIMULATION")
      nb_jobs = nb_fields

      # --- Create MKSIM output dataset
      dataset_base_dir = master.config.get_as_string("BASE_OUTPUT_DIR", "DIR.OUTPUT")
      dataset_name = "mksim"
      dataset_type = "Mpfx"  
      dataset_dir_list = []
      dataset_recurse_dirs = True   
      dataset = self.create_dataset(master, dataset_name, dataset_type, 
                                            dataset_base_dir, dataset_dir_list, 
                                            dataset_recurse_dirs) 

      # --- Create and populate the image path dictionary
      out_gal_image_pattern   =  master.config.get_as_string("OUTPUT_FILE_NAME", 
                                                   "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
      out_psf_image_pattern   =  master.config.get_as_string("OUTPUT_FILE_NAME", 
                                                   "PRIMARY_DATASET.IMAGE_PROPERTIES.STAR")
      out_gal_catalog_pattern =  master.config.get_as_string("OUTPUT_FILE_NAME", 
                                                   "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY")
      out_psf_catalog_pattern =  master.config.get_as_string("OUTPUT_FILE_NAME", 
                                                   "PRIMARY_DATASET.CATALOG_PROPERTIES.STAR")
      
      (_, gal_image_ext) = os.path.splitext(out_gal_image_pattern)
      (_, psf_image_ext) = os.path.splitext(out_psf_image_pattern)
      (_, gal_catalog_ext) = os.path.splitext(out_gal_catalog_pattern)
      (_, psf_catalog_ext) = os.path.splitext(out_psf_catalog_pattern)
      gal_image_type   = out_gal_image_pattern.split('-')[0] + gal_image_ext
      psf_image_type   = out_psf_image_pattern.split('-')[0] + psf_image_ext
      gal_catalog_type = out_gal_catalog_pattern.split('-')[0] + gal_catalog_ext
      psf_catalog_type = out_psf_catalog_pattern.split('-')[0] + psf_catalog_ext

      image_list = []
      image_range = []
      if master.config.has_key("IMAGE_LIST", "PRIMARY_DATASET.IMAGE_PROPERTIES"):
         image_list = master.config.get_as_list("IMAGE_LIST", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      if master.config.has_key("IMAGE_RANGE", "PRIMARY_DATASET.IMAGE_PROPERTIES"):
         image_range = master.config.get_as_list("IMAGE_RANGE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
         
      file_matcher = MpfxFileMatcher(master, dataset, out_gal_image_pattern, dataset_dir_list, 
                                             image_list, image_range, 
                                             sort=True, recurse=dataset_recurse_dirs)
      gal_image_patterns   =  master.config.get_as_list("FILENAME_PATTERNS", 
                                                       "PRIMARY_DATASET.IMAGE_PROPERTIES.GALAXY")
      psf_image_patterns   =  master.config.get_as_list("FILENAME_PATTERNS", 
                                                       "PRIMARY_DATASET.IMAGE_PROPERTIES.STAR")
      
      file_path_dico = {}

      for ijob in xrange(0, nb_jobs):
         job_key  = (".", ijob, 0)
         gal_image_name   = out_gal_image_pattern.format(ijob, 0)
         psf_image_name   = out_psf_image_pattern.format(ijob, 0)
         gal_catalog_name = out_gal_catalog_pattern.format(ijob, 0)
         psf_catalog_name = out_psf_catalog_pattern.format(ijob, 0)
         
         gal_image_filepath   = os.path.abspath(os.path.join(master.result_output_dir,
                                                             gal_image_name))
         psf_image_filepath   = os.path.abspath(os.path.join(master.result_output_dir, 
                                                             psf_image_name))
         gal_catalog_filepath = os.path.abspath(os.path.join(master.result_output_dir, 
                                                             gal_catalog_name))
         psf_catalog_filepath = os.path.abspath(os.path.join(master.result_output_dir, 
                                                             psf_catalog_name))
         
         # --- Make sure the to create jobs for matching files only
         matched = False
         for pattern in gal_image_patterns:
            matched = file_matcher.match_file(master, gal_image_filepath, gal_image_name, pattern)
            if matched:
               break

         if matched:
            for pattern in psf_image_patterns:
               matched = file_matcher.match_file(
                                             master, psf_image_filepath, psf_image_name, pattern)
               if matched:
                  break

         if not matched:
            continue   
         
         file_path_dico[job_key] = {}
         file_path_dico[job_key][gal_image_type]   = gal_image_filepath 
         file_path_dico[job_key][psf_image_type]   = psf_image_filepath
         file_path_dico[job_key][gal_catalog_type] = gal_catalog_filepath
         file_path_dico[job_key][psf_catalog_type] = psf_catalog_filepath

#       print "file_path_dico:", file_path_dico

      # --- Setup the values of the 'sequence' section keys
      sequence_dico = self.simulator.setup_sequence_config_dico(master)

      # --- Create the jobs
      self._job_list = [self.create_job(master, dataset, 
                                        file_path_dico, sequence_dico, *key) \
                                                         for key in sorted(file_path_dico.keys())]

      # --- Record the list of branches to which jobs belongs
      self._job_branches = list(set([j.branch for j in self._job_list]))      

      # --- Create output job directories
      output_dir_dico = master.config.get_section_data("DIR.OUTPUT")
      for (dir_key, dir_value) in output_dir_dico.items():
         if dir_key.startswith("OUTPUT_"):
            if dir_key.find("LOG") == -1:  # no tree for log dirs
               base_dir = os.path.join(master.run_output_dir, dir_value)
               self._create_dir_tree(base_dir, self._job_branches)

      return self._job_list

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):
      """! 
         Factory method for creating a Job object.
         Each subclass can define a job class derived from mpf.Job and instantiate it here.

         @param master the master process
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                data necessary to process the job
         @return a new MksJob job object
         @see Job, MksJob
      """

      return MksJob(master, dataset, self.helper, *args)

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """!
         Invoked by the Worker to perform some optional pre-processing on the job. 
         It may be for instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class MksJob to process
         @param worker instance of the worker process

         @return an object of class MksJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>before</b> process_job() is called
         @see preprocess_job, postprocess_job, Job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in mks

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """! 
         Process a job of class MksJob and return the corresponding results in the form of a 
         MksJobResult object to the Master.       

         @param job object of class MksJob with processed data
         @param worker Worker object instance

         @return an object of class MksJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.process_job() 
         @see Job, MpfxJob, JobResult, MpfxJobResult 
      """

      object_per_type_dico = {}  # will contain the job results per image type (images)

      # --- Iterate over all jobs to create the galaxy and PSF mosaics with their catalogs
      try:

#          # --- Check data integrity 
#          if not self._valid_file_types(job, worker):
#             if worker.logging_enabled():
#                worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
#                   "Check definition of galaxy and PSF image and catalog properties".format(
#                   worker.name, job.img_no, job.epoch))
            
         # --- Process each job...
         for file_type in job.get_file_types(): 

            filepath = job.get_file_path(file_type)

            if job.dataset.is_galaxy_image(worker, filepath):

               if not file_type in object_per_type_dico:
                  object_per_type_dico[file_type] = {}

               if worker.logging_enabled():
                  worker.logger.log_info_p("{0} - "\
                           "Creating simulation files {1}...".format(worker.name, job.file_names))
                  worker.logger.flush()

               # Create simulation image and catalogs
               success = self.simulator.run_simulations(file_type, job, worker)      
                  
               object_per_type_dico[file_type] = {}
               
      except:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                        "{0} - An error occurred while processing job: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                sys.exc_info()[1]))

      # --- Create a MksJobResult object with the results from all image types
      return MksJobResult(object_per_type_dico, job, worker)


   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.
         @param job_result object of class MksJobResult with processed data
         @param worker instance of the worker process

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in mks

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job.       

         @param job_result object of class MksJobResult with processed data
         @param master Master object instance
         @note overrides MpfxJobProcessor.process_job_result()
         @see Job, JobResult, MpfxJobResult
      """
  
      # TODO: add whatever is needed here: plots, stats... 
  
#       print("*** job result: {0}".format(job_result))
  
#       # --- Iterate through the results of all file types in the job...  
#       object_per_type_dico = job_result.result
#       for file_type in object_per_type_dico.keys():
#          se_xform_dico = object_per_type_dico[file_type]["se_xform_dico"]
# 
#          if len(se_xform_dico) > 0:
# 
#             # --- Log processing statistics
#             if master.logging_enabled():
#                msg = "{0} - /{1}/img {2:03}-{3:1d} - {4} - Catalog generated - "\
#                      "Object count: {5} -> {6} - Processor time: {7:.2f} sec"
#                master.logger.log_info_p(msg.format(
#                   master.name, job_result.job.get_branch_tree(), job_result.job.img_no, 
#                   job_result.job.epoch,  file_type, 
#                   se_xform_dico["initial_object_count"], se_xform_dico["final_object_count"],
#                   se_xform_dico["elamksd_time"]))
#                master.logger.flush()
# 
#             # --- Produce plots based on the job results   
#             if master.config.get_as_boolean("CREATE_SE_PLOTS", "DEBUGGING"):
#                output_plot_dir = os.path.join(master.plot_output_dir, 
#                                               job_result.job.get_branch_tree())
#                self.plotter.make_plots(job_result, output_plot_dir, master)     
# 
#             # --- Produce statistics based on the job results   
#             if master.config.get_as_boolean("CREATE_SE_STATS", "DEBUGGING"):
#                output_stat_dir = os.path.join(master.stat_output_dir, 
#                                               job_result.job.get_branch_tree())
#                self.helper.make_stats(job_result, output_stat_dir, master)     
      pass
      
   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed. 
         @param master instance of the Master
         @note Job results can be selectively obtained by calling the get_job_results() methods,
               specifying the relevant query criteria 
         @note the entire list of Job results can aldo be obtained with by calling 
               JobProcessor.job_result_list()
         @see Job, JobResult, get_job_result()
         @note MksJobResult objects can also be directly queried from the 
               mks.mks_job.MksJobProcessor.job_result_dico dictionary
      """

      pass


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _valid_file_types(self, job, worker):
 
      nb_star_image_types = len(job.get_star_image_file_types(worker))
      nb_gal_image_types  = len(job.get_galaxy_image_file_types(worker))
      nb_star_catalog_types = len(job.get_star_catalog_file_types(worker))
      nb_gal_catalog_types  = len(job.get_galaxy_catalog_file_types(worker))
 
      is_valid =  (nb_star_image_types   == 1 and nb_gal_image_types   == 1)
      is_valid &= (nb_star_catalog_types == 1 and nb_gal_catalog_types == 1) 

      if not is_valid: 
         if worker.logging_enabled():
            if nb_star_image_types != 1:
               worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Invalid or missing galaxy image: {3} - Check 'OUTPUT_FILE_NAME' key value; "\
                  "make sure 'OUTPUT_FILE_NAME' matches 'FILENAME_PATTERNS' key value".format(
                  worker.name, job.img_no, job.epoch,  job.get_star_image_file_types(worker)))
            if nb_gal_image_types != 1:
               worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Invalid or missing galaxy image: {3} - Check 'OUTPUT_FILE_NAME' key value; "\
                  "make sure 'OUTPUT_FILE_NAME' matches 'FILENAME_PATTERNS' key value".format(
                  worker.name, job.img_no, job.epoch, job.get_galaxy_image_file_types(worker)))
            if nb_star_catalog_types != 1:
               worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Invalid or missing star catalog: {3} - Check 'OUTPUT_FILE_NAME' key value; "\
                  "make sure 'OUTPUT_FILE_NAME' matches 'FILENAME_PATTERNS' key value".format(
                  worker.name, job.img_no, job.epoch, job.get_star_catalog_file_types(worker)))
            if nb_gal_catalog_types != 1:
               worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Invalid or missing galaxy catalog: {3} - Check 'OUTPUT_FILE_NAME' key value; "\
                  "make sure 'OUTPUT_FILE_NAME' matches 'FILENAME_PATTERNS' key value".format(
                  worker.name, job.img_no, job.epoch, job.get_galaxy_catalog_file_types(worker)))
                     
      return is_valid

   # -----------------------------------------------------------------------------------------------
   def _create_dir_tree(self, base_dir, branches):
      """! 
         Create the directory tree for storing the results.
         @param base_dir base directory
         @param branches list of ordered directory names, making a path tree
         @note example of directory layout: control/ground/constant
         @see configuration file
      """      

      for branch in branches:
         dir_names = branch.split("/")
         branch_path = base_dir
         for dir_name in dir_names:
            branch_path = os.path.join(branch_path, dir_name)
            self.helper.make_dir(branch_path)


# -------------------------------------------------------------------------------------------------
class MksJob(MpfxJob):      

   """! 
       Represents a Mks job. Based on mpf.Job and mpfcs82.MpfxJob.  
   """

   def __init__(self, master, dataset, helper, img_path_dico, sequence_dico, *args):
      """! 
         Construct a new job for identifiying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:

         @code self.img_path_dico[type]@endcode with one of the aforementioned types
         @param master instance of the Master
         @param dataset data source used to create the job
         @param helper helper object for this class
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param sequence_dico doctionary parameters to initialize as a sequence
         @param args a list of extra parameters

         @return an initialized MksJob object 
      """

      # --- Construct a new Job

      MpfxJob.__init__(self, master, dataset, img_path_dico, *args)

      # --- Postage stamp size of the source image
      self._stamp_size = master.config.get_as_int("GAL_STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      self._helper = helper
      
      # --- Compute the actual values of sequences defined in the configuration
      self._sequence_value_dico = {}   
      for key in sequence_dico:
         cycle = sequence_dico[key]
         self._sequence_value_dico[key] = cycle.next()

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the class' helper object """
      return self._helper

   @property
   def stamp_size(self):
      """! @return the postage stamp size """
      return self._stamp_size

   @property
   def sequence_value_dico(self):
      """! @return the dictionary of pre-computed sequence values """
      return self._sequence_value_dico

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_stamp_size(self):
      """! @return the postage stamp size """

      if self.stamp_size != -1:
         return self.stamp_size
      else:
         return MpfxJob.get_stamp_size(self)    # default stamp size of the job

   # -----------------------------------------------------------------------------------------------
   def get_stamp_center(self):
      """ ! @return the postage stamp geometrical center pixel no, indexed from zero """      
      
      return self.helper.get_stamp_center(self.get_stamp_size())




# -------------------------------------------------------------------------------------------------
class MksJobResult(MpfxJobResult):      

   """! 
       Represents the result of a GREAT3 job. See mpfcs82.MpfxJobResult and mpf.JobResult.
   """
   
   # -----------------------------------------------------------------------------------------------
   def __init__(self, result, job, worker):   

      """! Construct a MksJobResult job result object """
      MpfxJobResult.__init__(self, worker, job, result)

      # --- Cached Job result data and associated stats
      self._helper = MksHelper()  # helper utility functions
#       self._data_dico  = self.helper.collect_catalog_data(result, job, worker)
#       self._stats_dico = self.helper.compute_stats(result, self._data_dico, job, worker)

      print("{0}".format(self))

   # -----------------------------------------------------------------------------------------------
   def __str__(self):
      """!
          Formatting for display 
          @return string representation of the object  
      """
      return MpfxJobResult.__str__(self)
   
   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the MksHelper instance. """
      return self._helper

#    @property
#    def data_dico(self):
#       """! @return the data dictionary """
#       return self._data_dico
# 
#    @property
#    def stats_dico(self):
#       """! @return the statistics dictionary """
#       return self._stats_dico


# -- EOF mks_job.py
