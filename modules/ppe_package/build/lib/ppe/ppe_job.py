"""!
   @package ppe.ppe_job Job Management
   @author Marc Gentile
   @file ppe_job.py
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
from ppe_psfex import *            # PSFEx catalog management
from ppe_plot import *           # plotter
from ppe_help import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class PpeJobProcessor(MpfxJobProcessor):

   """!
      Job processor: submit jobs and process associated job results. Based on mpf.JobProcessor.
   """

   def __init__(self, master):

      """! Job Processor constructor """
      MpfxJobProcessor.__init__(self, master)

      self._helper = PpeHelper()                         # helper utility functions
      self._plotter = PpePlotter()                       # plotter

      self._pe_runner = PSFExRunner(self)           # PSFEx execution
      self._pe_processor = PSFExProcessor(self)     # PSFEx catalog processor

   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~

   @property
   def pe_runner(self):
      """! @return the PSFEx runner instance. """
      return self._pe_runner

   @property
   def pe_processor(self):
      """! @return the PSFEx processor instance. """
      return self._pe_processor

   @property
   def helper(self):
      """! @return the PpeHelper instance. """
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

      # SF NOTE: this method is called by mpfx_job.py MpfxJobProcessor.create_jobs()
      # No apparent connection to specific package, could be part of template.

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

      # SF NOTE: This method is called by mp_calc.py Master.run() and the
      # JobProcessor is defined in the ..._SMP.py module.
      # No apparent connection to specific package, could be part of template.

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
         @return a new PpeJob job object
         @see Job, PpeJob
      """

      # SF NOTE: This method is called by mpfx_job.py
      # MpfxJobProcessor.create_jobs() and Simply returns an instance of
      # PpeJob().Maybe package specific

      return PpeJob(master, dataset, self.helper, *args)

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """!
         Process a job of class PpeJob and return the corresponding results in the form of a
         PpeJobResult object to the Master.

         @param job object of class PpeJob with processed data
         @param worker Worker object instance

         @return an object of class PpeJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.process_job()
         @see Job, MpfxJob, JobResult, MpfxJobResult
      """

      # SF NOTE: This method is called by mp_calc_SMP.py
      # MasterSMP.process_jobs() and is package specific.


      # --- This is where we add the code to process the job

      # --- Depending on the object selection criteria, we may have to process either catalogs
      #     matching "*catalog-???.fits" or images matching "*image-???-??.fits"
      #     e.g. "star_catalog-002-0.fits" or "starfield_image-002-0.fits" for regular stars.
      #     The class attribute: job.img_path_dico contains all these paths as a dictionary
      #     (see description of class: PpeJob below).

      # --- Here we only want to process images and catalogs. We don't any other files accesible
      #     from this job, like deep_images or other dither files.

      object_per_type_dico = {}  # will contain the job results per image type (images)

      try:

          for file_type in job.get_file_types():    # check the type of file (looking for images)
              filepath = job.get_file_path(file_type)

              if not file_type in object_per_type_dico:
                 object_per_type_dico[file_type] = {}

              if worker.logging_enabled():
                 worker.logger.log_info_p("{0} - Processing SExtractor Output File {1}...".format(worker.name, filepath))
                 worker.logger.flush()

              pe_run_dico = self.pe_runner.run_PSFEx(file_type, job, worker) # THIS LINE CALLS PSFEX
              pe_output_cat_filepath   = pe_run_dico["pe_output_cat_filepath"]
              pe_output_check_filepath = pe_run_dico["pe_output_check_filepath"]
              object_per_type_dico[file_type]["pe_run_dico"] = pe_run_dico

              initial_results = PpeJobResult(object_per_type_dico, job, worker)

              if worker.config.get_as_boolean("CREATE_INITIAL_PE_PLOTS", "DEBUGGING"):
                 plot_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
                 self.plotter.make_plots(initial_results, plot_output_dir, worker,
                                                                           plot_prefix="initial_")

              if worker.config.get_as_boolean("CREATE_INITIAL_PE_STATS", "DEBUGGING"):
                 stat_output_dir = os.path.join(worker.stat_output_dir, job.get_branch_tree())
                 self.helper.make_stats(initial_results, stat_output_dir, worker,
                                                                          stat_prefix="initial_")


              # SF NOTE: NEED TO FIX THIS
              object_per_type_dico[file_type]["pe_xform_dico"] = {'initial_object_count': 1, 'final_object_count': 1, 'elapsed_time': 1}

            #                                                               # --- Process the generated PSFEx catalog
            #   object_per_type_dico[file_type]["pe_xform_dico"] =\
            #                                self.pe_processor.process_catalog(
            #                                                pe_output_cat_filepath,
            #                                                pe_output_check_filepath,
            #                                                file_type, job, worker)

              # --- Job result to return
              job_result = PpeJobResult(object_per_type_dico, job, worker)

      except:

          if worker.logging_enabled():
              worker.logger.log_error_p("{0} - Some error occurred while processing job: {1} ({2})".format(worker.name, job, sys.exc_info()[1]))

      # --- Create a PpeJobResult object with the results from all image types
      print "PSF JOB PROCESSED"
      return PpeJobResult(object_per_type_dico, job, worker)

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """!
         Process the result associated with a processed job.

         @param job_result object of class PpeJobResult with processed data
         @param master Master object instance
         @note overrides MpfxJobProcessor.process_job_result()
         @see Job, JobResult, MpfxJobResult
      """

      # SF NOTE: This method is called by mp_calc_SMP.py
      # MasterSMP.process_jobs() and is package specific

      # --- Iterate through the results of all file types in the job...
      object_per_type_dico = job_result.result
      for file_type in object_per_type_dico.keys():

         pe_xform_dico = object_per_type_dico[file_type]["pe_xform_dico"]

         if len(pe_xform_dico) > 0:

            # --- Log processing statistics
            if master.logging_enabled():
               msg = "{0} - /{1}/img {2:03}-{3:1d} - {4} - Catalog generated - "\
                     "Object count: {5} -> {6} - Processor time: {7:.2f} sec"
               master.logger.log_info_p(msg.format(
                  master.name, job_result.job.get_branch_tree(), job_result.job.img_no,
                  job_result.job.epoch,  file_type,
                  pe_xform_dico["initial_object_count"], pe_xform_dico["final_object_count"],
                  pe_xform_dico["elapsed_time"]))
               master.logger.flush()

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
         @note PpeJobResult objects can also be directly queried from the
               ppe.ppe_job.PpeJobProcessor.job_result_dico dictionary
      """

      # SF NOTE: This method simply overrides the mpfx_job.py
      # MpfxJobProcessor.all_jobs_processed() method. Could certainly be
      # part of a template

      pass

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

      # SF NOTE: this method does not appear to package specific, could be
      # part of a template


      for branch in branches:
         dir_names = branch.split("/")
         branch_path = base_dir
         for dir_name in dir_names:
            branch_path = os.path.join(branch_path, dir_name)
            self.helper.make_dir(branch_path)


# -------------------------------------------------------------------------------------------------
class PpeJob(MpfxJob):

   """!
       Represents a Ppe job. Based on mpf.Job and mpfcs82.MpfxJob.
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

         @return an initialized PpeJob object
      """

      # SF NOTE: Adds methods to the MpfxJob() class.
      # Package specific


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

# -------------------------------------------------------------------------------------------------
class PpeJobResult(MpfxJobResult):

   """!
       Represents the result of a GREAT3 job. See mpfcs82.MpfxJobResult and mpf.JobResult.
   """

   # SF NOTE: Not sure if this class is package specific

   def __init__(self, result, job, worker):

      """! Construct a PpeJobResult job result object """
      MpfxJobResult.__init__(self, worker, job, result)

      # --- Cached Job result data and associated stats
    #   self._helper = PpeHelper()  # helper utility functions
    #   self._data_dico  = self.helper.collect_catalog_data(result, job, worker)
    #   self._stats_dico = self.helper.compute_stats(result, self._data_dico, job, worker)

   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~

   # @property
   # def helper(self):
   #    """! @return the PpeHelper instance. """
   #    return self._helper
   #
   # @property
   # def data_dico(self):
   #    """! @return the data dictionary """
   #    return self._data_dico
   #
   # @property
   # def stats_dico(self):
   #    """! @return the statistics dictionary """
   #    return self._stats_dico


# -- EOF ppe_job.py
