"""! 
   @package pse.pse_job Job Management
   @author Marc Gentile
   @file pse_job.py
   Job Management
""" 

# -- Python imports
import os, sys
import time
import numpy
from operator import itemgetter, attrgetter

# -- External imports
from mpfx.mpfx_job import *    # base job processing

# --- Module-specific imports
from pse_sex import *            # SExtractor catalog management
from pse_plot import *           # plotter 
from pse_help import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class PseJobProcessor(MpfxJobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpf.JobProcessor.   
   """

   def __init__(self, master):

      """! Job Processor constructor """
      MpfxJobProcessor.__init__(self, master)
      
      self._helper = PseHelper()                         # helper utility functions
      self._plotter = PsePlotter()                       # plotter

      self._se_runner = SExtractorRunner(self)           # SExtractor execution
      self._se_processor = SExtractorProcessor(self)     # SExtractor catalog processor

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def se_runner(self):
      """! @return the SExtractor runner instance. """
      return self._se_runner

   @property
   def se_processor(self):
      """! @return the SExtractor processor instance. """
      return self._se_processor

   @property
   def helper(self):
      """! @return the PseHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the Plotter instance. """
      return self._plotter

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def create_dataset(self, master, dataset_name, dataset_type, 
                                    dataset_base_dir, dataset_dir_list, dataset_recurse):
      """! 
         Create a primary Dataset object that represent the data source for images and catalogs.
         @param master master process instance
         @param dataset_name dataset name
         @param dataset_type prefix of a Dataset class, assumed to be of the form 
                @code <dataset_type>Dataset @endcode, like @c MpfxDataset
         @param dataset_base_dir dataset base directory
         @param dataset_dir_list [optional] a list of specific dirs to search under base directory
         @param dataset_recurse [optional] tell whether to walk down directories (default @c True)
         @return the Dataset instance 
      """

      dataset_class = self.get_dataset_module(master, dataset_name, dataset_type, dataset_base_dir)
      if dataset_class is not None:
         return dataset_class(master, dataset_name, 
                              dataset_base_dir, dataset_dir_list, dataset_recurse)   
      else:
         return None

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """! 
         Locate all objects to process and create the corresponding jobs. 
         @param master Master object instance
         @return the list of created jobs 
      """

      job_list = MpfxJobProcessor.create_jobs(self, master)

      # --- Create output directories for job outputs:
      output_dir_dico = master.config.get_section_data("DIR.OUTPUT")
      for (dir_key, dir_value) in output_dir_dico.items():
         if dir_key.startswith("OUTPUT_"):
            if dir_key.find("LOG") == -1:  # no tree for log dirs
               base_dir = os.path.join(master.run_output_dir, dir_value)
               self._create_dir_tree(base_dir, self._job_branches)

      return job_list

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):
      """! 
         Factory method for creating a Job object.
         Each subclass can define a job class derived from mpf.Job and instanciate it here.
         @param master the master process
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                data necessary to process the job
         @return a new PseJob job object
         @see Job, PseJob
      """

      return PseJob(master, dataset, self.helper, *args)

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """!
         Invoked by the Worker to perform some optional pre-processing on the job. 
         It may be for instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class PseJob to process
         @param worker instance of the worker process

         @return an object of class PseJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>before</b> process_job() is called
         @see preprocess_job, postprocess_job, Job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in pse

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """! 
         Process a job of class PseJob and return the corresponding results in the form of a 
         PseJobResult object to the Master.       

         @param job object of class PseJob with processed data
         @param worker Worker object instance

         @return an object of class PseJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.process_job() 
         @see Job, MpfxJob, JobResult, MpfxJobResult 
      """

      # --- This is where we add the code to process the job

      # --- Depending on the object selection criteria, we may have to process either catalogs
      #     matching "*catalog-???.fits" or images matching "*image-???-??.fits"
      #     e.g. "star_catalog-002-0.fits" or "starfield_image-002-0.fits" for regular stars.
      #     The class attribute: job.img_path_dico contains all these paths as a dictionary
      #     (see description of class: PseJob below).  

      # --- Here we only want to process images and catalogs. We don't any other files accesible
      #     from this job, like deep_images or other dither files. 

      object_per_type_dico = {}  # will contain the job results per image type (images)

      # --- Iterate over the path dictionary content to generate and process SEXtractor catalogs
      try:

         for file_type in job.get_file_types():    # check the type of file (looking for images) 

            filepath = job.get_file_path(file_type)

            if job.dataset.is_image(worker, filepath):

               if not file_type in object_per_type_dico:
                  object_per_type_dico[file_type] = {}

               if worker.logging_enabled():
                  worker.logger.log_info_p("{0} - Processing image {1}...".format(
                                                                       worker.name, filepath))
                  worker.logger.flush()

               # --- Run SExtractor to analyse the image and produce its catalog
               se_run_dico = self.se_runner.run_SExtractor(file_type, job, worker)
               se_output_cat_filepath   = se_run_dico["se_output_cat_filepath"]
               se_output_check_filepath = se_run_dico["se_output_check_filepath"]
               object_per_type_dico[file_type]["se_run_dico"] = se_run_dico

               # --- Produce initial plots and statistics before the initial SE catalog is processed
               initial_results = PseJobResult(object_per_type_dico, job, worker)
               if worker.config.get_as_boolean("CREATE_INITIAL_SE_PLOTS", "DEBUGGING"):
                  plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
                  self.plotter.make_plots(initial_results, plot_output_dir, worker,
                                                                            plot_prefix="initial_")
               
               if worker.config.get_as_boolean("CREATE_INITIAL_SE_STATS", "DEBUGGING"):
                  stat_output_dir = os.path.join(worker.stat_output_dir, job.get_branch_tree())
                  self.helper.make_stats(initial_results, stat_output_dir, worker,
                                                                           stat_prefix="initial_")
               # --- Process the generated SExtractor catalog
               object_per_type_dico[file_type]["se_xform_dico"] =\
                                            self.se_processor.process_catalog(
                                                            se_output_cat_filepath, 
                                                            se_output_check_filepath,
                                                            file_type, job, worker)
      except:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                        "{0} - Some error occurred while processing job: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                sys.exc_info()[1]))

      # --- Create a PseJobResult object with the results from all image types
      return PseJobResult(object_per_type_dico, job, worker)


   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.
         @param job_result object of class PseJobResult with processed data
         @param worker instance of the worker process

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in pse

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job.       

         @param job_result object of class PseJobResult with processed data
         @param master Master object instance
         @note overrides MpfxJobProcessor.process_job_result()
         @see Job, JobResult, MpfxJobResult
      """
  
      # --- Iterate through the results of all file types in the job...  
      object_per_type_dico = job_result.result
      for file_type in object_per_type_dico.keys():
         se_xform_dico = object_per_type_dico[file_type]["se_xform_dico"]

         if len(se_xform_dico) > 0:

            # --- Log processing statistics
            if master.logging_enabled():
               msg = "{0} - /{1}/img {2:03}-{3:1d} - {4} - Catalog generated - "\
                     "Object count: {5} -> {6} - Processor time: {7:.2f} sec"
               master.logger.log_info_p(msg.format(
                  master.name, job_result.job.get_branch_tree(), job_result.job.img_no, 
                  job_result.job.epoch,  file_type, 
                  se_xform_dico["initial_object_count"], se_xform_dico["final_object_count"],
                  se_xform_dico["elapsed_time"]))
               master.logger.flush()

            # --- Produce plots based on the job results   
            if master.config.get_as_boolean("CREATE_SE_PLOTS", "DEBUGGING"):
               output_plot_dir = os.path.join(master.plot_output_dir, 
                                              job_result.job.get_branch_tree())
               self.plotter.make_plots(job_result, output_plot_dir, master)     

            # --- Produce statistics based on the job results   
            if master.config.get_as_boolean("CREATE_SE_STATS", "DEBUGGING"):
               output_stat_dir = os.path.join(master.stat_output_dir, 
                                              job_result.job.get_branch_tree())
               self.helper.make_stats(job_result, output_stat_dir, master)     


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
         @note PseJobResult objects can also be directly queried from the 
               pse.pse_job.PseJobProcessor.job_result_dico dictionary
      """

      pass
      
#      # --- Check that the actual number of objects in the generated SExtractor catalog is correct
#      all_checked_objects = []
#      wrong_job_results = []
#      all_jobs_results = self.get_job_results()
#      if len(all_jobs_results) > 0:
#         for job_result in all_jobs_results:
#            object_per_type_dico = job_result.result

#            if len(object_per_type_dico[t]["se_xform_dico"]) > 0:

#               actual_object_counts = [object_per_type_dico[t]["se_xform_dico"]["actual_object_count"] \
#                                       for t in object_per_type_dico.keys()]
#               final_object_counts  = [object_per_type_dico[t]["se_xform_dico"]["final_object_count"] \
#                                       for t in object_per_type_dico.keys()]

#               checked_objects = numpy.equal(actual_object_counts, final_object_counts)
#               if not numpy.all(checked_objects):
#                  wrong_job_results.append(job_result)

#               all_checked_objects.extend(checked_objects)

#         all_checked_objects = numpy.asarray(all_checked_objects)
#         if not numpy.all(all_checked_objects):

#            # --- Found a discrepancy between the initial (expected) and final object counts
#            wrong_counts = numpy.where(all_checked_objects == False)[0]
#            wrong_objects =\
#             ["/{0}:{1:03d}-{2:1d}".format(
#               r.job.get_branch_tree(), r.job.img_no, r.job.epoch) for r in wrong_job_results]

#            if master.logging_enabled():
#               master.logger.log_warning_p(
#                "{0} - Found {1} catalogs with wrong object counts: {2}".format(
#                                                                  master.name, 
#                                                                  len(wrong_counts), wrong_objects))
#         else:
#            if master.logging_enabled():
#               master.logger.log_info_p(
#                  "{0} - *** ALL catalogs have the correct object count ***".format(master.name) )          
#               master.logger.flush()

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

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
class PseJob(MpfxJob):      

   """! 
       Represents a Pse job. Based on mpf.Job and mpfcs82.MpfxJob.  
   """

   def __init__(self, master, dataset, helper, img_path_dico, *args):
      """! 
         Construct a new job for identiofiying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:

         @code self.img_path_dico[type]@endcode with one of the aforementioned types
         @param master instance of the Master
         @param dataset data source used to create the job
         @param helper helper object for this class
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param args a list of extra parameters

         @return an initialized PseJob object 
      """

      # --- Construct a new Job

      MpfxJob.__init__(self, master, dataset, img_path_dico, *args)

      # --- Postage stamp size of the source image
      self._stamp_size = master.config.get_as_int("STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      self._helper = helper

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

#      stamp_size = self.get_stamp_size()
#      if stamp_size % 2 == 0:
#         # Even size
#         center_pixel = int(stamp_size / 2.0 - 1.0)   # e.g. 48x48 -> 23x23
#      else:
#         # Odd size
#         center_pixel = int(stamp_size/ 2.0)          # e.g. 47x47 -> 23x23

#      return center_pixel

# -------------------------------------------------------------------------------------------------
class PseJobResult(MpfxJobResult):      

   """! 
       Represents the result of a GREAT3 job. See mpfcs82.MpfxJobResult and mpf.JobResult.
   """

   def __init__(self, result, job, worker):   


      """! Construct a PseJobResult job result object """
      MpfxJobResult.__init__(self, worker, job, result)

      # --- Cached Job result data and associated stats
      self._helper = PseHelper()  # helper utility functions
      self._data_dico  = self.helper.collect_catalog_data(result, job, worker)
      self._stats_dico = self.helper.compute_stats(result, self._data_dico, job, worker)

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the PseHelper instance. """
      return self._helper

   @property
   def data_dico(self):
      """! @return the data dictionary """
      return self._data_dico

   @property
   def stats_dico(self):
      """! @return the statistics dictionary """
      return self._stats_dico


# -- EOF pse_job.py
