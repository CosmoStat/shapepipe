"""! 
   Job Processing
"""

# -- Python imports
import os, sys
import time
import numpy
from operator import itemgetter, attrgetter
from copy import deepcopy

# -- External imports
from mpfx.mpfx_job import *      # base job processing
from mpfx.mpfx_data import *     # data access

# --- Module-specific imports
from gfit_shape import *         # galaxy shape estimation
from gfit_plot import *          # plotter 
from gfit_helper import *        # helper utility functions

# -------------------------------------------------------------------------------------------------
class GfitJobProcessor(MpfxJobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpfx.JobProcessor.   
   """

   def __init__(self, master):
      """! 
         Job Processor constructor
         @param master master process instance
      """

      MpfxJobProcessor.__init__(self, master)
      
      self._shape_estimator = GfitShapeEstimator(master) 
      self._plotter = GfitPlotter()    # plotter
      self._helper = GfitHelper()      # helper

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def shape_estimator(self):
      """! @return the shape estimator instance. """
      return self._shape_estimator

   @property
   def plotter(self):
      """! @return the GfitPlotter instance. """
      return self._plotter

   @property
   def helper(self):
      """! @return the GfitHelper instance. """
      return self._helper


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
         @param dataset_dir_list [optional] a list of specific directories to search under
                base directory
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
         @param master Master process instance
         @return the list of created jobs 
      """

      # --- Check that the jobs contains all necessary data 
      if master.logging_enabled():
         master.logger.log_info_p("{0} - Creating jobs...".format(master))
         master.logger.flush()
         
      is_valid = True
      job_list = MpfxJobProcessor.create_jobs(self, master)
      if len(job_list) > 0:
         # --- Check first job
         job = job_list[0]
         is_valid = len(job.galaxy_se_filepath_dico) > 0 and \
                    len(job.star_se_filepath_dico) > 0 and \
                    len(job.psf_filepath_dico) > 0

         if master.logging_enabled() and not is_valid:
            if len(job.galaxy_se_filepath_dico) == 0:
               file_patterns = master.config.get_as_list("FILE_PATTERNS",
                                                    "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES")
               master.logger.log_error_p("{0} - could not find to find input "\
                                            "SExtractor Galaxy catalogs {1}".format(
                                            master, file_patterns))
            if len(job.star_se_filepath_dico) == 0:
               file_patterns = master.config.get_as_list("FILE_PATTERNS",
                                                      "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES")
               master.logger.log_error_p("{0} - could not find input "\
                                         "SExtractor Star catalogs {1}".format(
                                            master, file_patterns))
            if len(job.psf_filepath_dico) == 0:
               file_patterns = master.config.get_as_list("FILE_PATTERNS",
                                                         "PSF_DATASET.IMAGE_PROPERTIES")
               master.logger.log_error_p("{0} - could not find input PSF images {1}".format(
                                             master, file_patterns))

      if not is_valid:
         self._job_list = []
         self._job_branches = []
         return self._job_list

      # --- Create output directories for job outputs:
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
         Factory method for creating a GfitJob object. 
         @param master the master process
         @param dataset source dataset object
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                     data necessary to process the job
         @return a new GfitJob job object
         @note overrides MpfxJobProcessor.create_job() 
         @see Job, MpfxJob
      """

      # --- create a MpfxJob instance, based on parent class mpf.Job
      job = GfitJob(master, dataset, *args)  

      # --- Attach to the job the datasets to be used
      
      # Source PSF datasets and SExtractor datasets
      se_galaxy_dataset = self._create_galaxy_SE_dataset(master)
      se_psf_dataset    = self._create_psf_SE_dataset(master)
      psf_image_dataset = self._create_psf_image_dataset(master)
      
      # Background sky noise specified as external images
      gal_sky_dataset   = self._create_galaxy_sky_dataset(master)
      psf_sky_dataset   = self._create_psf_sky_dataset(master)

      # Galaxy and PSF object noise specified as external rms maps
      gal_noise_dataset   = self._create_galaxy_noise_dataset(master)
      psf_noise_dataset   = self._create_psf_noise_dataset(master)

      job.psf_filepath_dico = self._locate_psf_image_files(psf_image_dataset, job, master)

      # --- Attach to the job the file paths of the relevant SEXtractor catalogs and PSF images
      job.galaxy_se_filepath_dico = self._locate_galaxy_se_files(se_galaxy_dataset, job, master)
      job.star_se_filepath_dico   = self._locate_psf_se_files(se_psf_dataset, job, master)

      # --- Attach to the job the file paths of sky the images for galaxies and PSFs if specified
      if gal_sky_dataset is not None:
         gal_sky_section_name = "SKY_MODEL_GALAXY.EXT_IMAGE"         
         job.galaxy_sky_filepath_dico = self._locate_sky_image_files(
                                                            gal_sky_dataset, gal_sky_section_name,
                                                            job, master)
      if psf_sky_dataset is not None:
         psf_sky_section_name = "SKY_MODEL_PSF.EXT_IMAGE"         
         job.psf_sky_filepath_dico = self._locate_sky_image_files(
                                                             psf_sky_dataset, psf_sky_section_name,
                                                             job, master)

      # --- Attach to the job the file paths of galaxies and PSFs rms maps if specified
      if gal_noise_dataset is not None:
         gal_noise_section_name = "GALAXY_NOISE.EXT_NOISE_MAP"
         job.galaxy_noise_map_filepath_dico = self._locate_noise_map_files(
                                                         gal_noise_dataset, gal_noise_section_name, 
                                                         job, master)
      if psf_noise_dataset is not None:
         psf_noise_section_name = "PSF_NOISE.EXT_NOISE_MAP"
         job.psf_noise_map_filepath_dico = self._locate_noise_map_files(
                                                         psf_noise_dataset, psf_noise_section_name, 
                                                         job, master)

      print job.img_path_dico

      print "job.galaxy_se_filepath_dico:", job.galaxy_se_filepath_dico
      print "job.star_se_filepath_dico:", job.star_se_filepath_dico
      print "job.psf_image_filepath_dico:", job.psf_filepath_dico
      print "job.galaxy_sky_filepath_dico:", job.galaxy_sky_filepath_dico
      print "job.psf_sky_filepath_dico:", job.psf_sky_filepath_dico
      print "job.galaxy_noise_map_filepath_dico", job.galaxy_noise_map_filepath_dico
      print "job.psf_noise_map_filepath_dico:", job.psf_noise_map_filepath_dico

      #print "job.img_path_dico:", job.img_path_dico

      self.dataset_dico = {dataset.name : dataset,\
                           se_galaxy_dataset.name : se_galaxy_dataset,\
                           se_psf_dataset.name : se_psf_dataset,\
                           psf_image_dataset.name : psf_image_dataset,\
                          }

      return job
      
   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """!
         Invoked by the Worker to perform some optional pre-processing on the job. 
         It may be for instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class GfitJob to process
         @param worker instance of the worker process

         @return an object of class GfitJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>before</b> process_job() is called
         @see preprocess_job, postprocess_job, Job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in gfit

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """! 
         Process a job of class GfitJob and return the corresponding results in the form of a 
         GfitJobResult object to the Master.       

         @param job object of class GfitJob with processed data
         @param worker Worker object instance

         @return an object of class GfitJobResult containing the data of the processed job

         @note overrides MpfxJobProcessor.process_job() 
         @see Job, MpfxJob, JobResult, MpfxJobResult 
      """

      # --- Depending on the object selection criteria, we may have to process either catalogs
      #     matching "*catalog-???.fits" or images matching "*image-???-??.fits"
      #     e.g. "star_catalog-002-0.fits" or "starfield_image-002-0.fits" for regular stars.
      #     The class attribute: job.img_path_dico contains all these paths as a dictionary
      #     (see description of class: GfitJob below).  

      # --- Here we only want to process images and catalogs. We don't any other files accesible
      #     from this job, like deep_images or other dither files. 

      object_per_type_dico = {}  # will contain all results (e1, e2, etc.) per image type   


      # --- Iterate over the path dictionary content and compute the ellipticities
      try:

         # --- For each file type referenced by the job...
         file_types = job.get_file_types()
         for file_type in file_types: 

            filepath = job.get_file_path(file_type)

            # --- A galaxy image to process
            if filepath is not None and job.dataset.is_galaxy_image(worker, filepath): #FCS ADDED FIRST CONDITION
               if worker.logging_enabled():
                  worker.logger.log_info_p(
                         "{0} - Processing {1}/image-{2:03d}-{3:1d} - {4}...".format(
                            worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))

               try:

                  shape_request = ShapeMeasurementRequest(file_type, job, worker)
                  if not shape_request is None:

                     # --- Process the ShapeMeasurementRequest object and produce shape measurements
                     object_per_type_dico[file_type] = self.shape_estimator.measure_galaxy_shapes(
                                                                                 shape_request,
                                                                                 job, worker)
                     
                     result_dico = object_per_type_dico[file_type]["result"]   # results
                     error_dico  = object_per_type_dico[file_type]["error"]    # measurement error
                     filter_dico = object_per_type_dico[file_type]["filtered"] # filtered objects
                     info_dico   = object_per_type_dico[file_type]["info"]     # info. measurement 
                     model_dico  = object_per_type_dico[file_type]["model"]    # info. model   
                     
                     if worker.logging_enabled():
                        worker.logger.log_info_p(
                          "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Generating catalogs...".format(
                          worker, job.get_branch_tree(), job.img_no, job.epoch,  file_type))

                     # --- Save the shape measurement results for this job
                     self._create_result_catalog(result_dico, error_dico, filter_dico, info_dico,
                                                job, worker)
            
                     # --- Create catalogs with all collected data (successful, filtered, in error)
                     if worker.config.get_as_boolean("DUMP_COLLECTED_DATA", "DATA_COLLECTION"):
                        self._create_data_catalog(
                                 result_dico, error_dico, filter_dico, info_dico, model_dico,
                                 file_type, job, worker)
 
                         # --- Summary plots at job level
                     if worker.config.get_as_boolean("PLOT_COLLECTED_DATA", "DATA_COLLECTION"):
                        self.plotter.create_image_plots(object_per_type_dico, file_type, 
                                                        job, worker)         

               except GfitShapeEstimator.ShapeMeasurementError as detail:
                  if worker.logging_enabled():
                     worker.logger.log_info_p(detail)

                  continue    # nevertheless, proceed with the next job
      

      except Exception as detail:
         if worker.logging_enabled():
            worker.logger.log_error_p(
                        "{0} - Some error occurred while processing job: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                detail))

      # --- Create a GfitJobResult object with the results from all image types
      return GfitJobResult(job, object_per_type_dico, worker)

   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.
         @param job_result object of class GfitJobResult with processed data
         @param worker instance of the worker process

         @note overrides MpfxJobProcessor.preprocess_job() 
               to perform any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult, MpfxJobResult
      """
      pass  # we no nothing here in gfit

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job.       

         @param job_result object of class GfitJobResult with processed data
         @param master Master object instance
         @note overrides MpfxJobProcessor.process_job_result()
         @see Job, JobResult, MpfxJobResult
      """

      # --- Iterate through the results of all file types in the job...  
      job = job_result.job
      object_per_type_dico = job_result.result

      for file_type in object_per_type_dico.keys():

         object_dico = object_per_type_dico[file_type]
         result_dico = object_per_type_dico[file_type]["result"]   # shape measurement results
         error_dico  = object_per_type_dico[file_type]["error"]    # shape measurement error
         filter_dico = object_per_type_dico[file_type]["filtered"] # object was filtered
         info_dico   = object_per_type_dico[file_type]["info"]     # info. about the measurement 
         model_dico  = object_per_type_dico[file_type]["model"]    # info. about the model

         # --- Gives a status of what has been processed at job level
         if master.logging_enabled():
            elapsed_time = info_dico["elapsed_time"]
            hours = elapsed_time/3600.0

            if info_dico["success_count"] > 0:
               # --- There are successsful results...
               bad_entry_marker = master.config.get_as_float("FAILED_ELLIPTICITY_VALUE", 
                                                              "SHAPE_MEASUREMENT")
               [mean_e1, mean_e2, mean_e] = self.shape_estimator.compute_average_ellipticity(
                                                                                 result_dico,
                                                                                 bad_entry_marker)

               processed_count = info_dico["success_count"] +\
                                  info_dico["failure_count"] + info_dico["filtered_count"]
#                processed_count = info_dico["success_count"] + info_dico["warning_count"] +\
#                                   info_dico["failure_count"] + info_dico["filtered_count"]

               msg = "{0} - {1}/img {2:03}-{3:1d} - {4} - Galaxies successfully processed: "\
                     "{5}/{6} ({7}-{8}-{9}-{10}) - <e1>={11:+.9f} <e2>={12:+.9f} <e>={13:+.9f} - "\
                     "Total Processing time: {13:.3f} secs. ({15:.3f} hours)"
               master.logger.log_info_p(msg.format(
                  master.name, job.get_branch_tree(), job.img_no, 
                  job.epoch,  file_type, 
                  info_dico["success_count"], processed_count, # info_dico["total_count"], 
                  info_dico["success_count"], info_dico["warning_count"], 
                  info_dico["failure_count"], info_dico["filtered_count"], 
                  mean_e1, mean_e2, mean_e, elapsed_time, hours))

               avg_time_per_gal = elapsed_time / processed_count

               master.logger.log_info_p("{0} - {1}/img {2:03}-{3:1d} - {4} - "\
                                        "Average time per galaxy: {5:.2f} secs.".format(
                            master.name, job.get_branch_tree(), job.img_no, 
                            job.epoch,  file_type, avg_time_per_gal))
            else:
               msg = "{0} - {1}/img {2:03}-{3:1d} - {4} - "\
                     "None of the galaxies could be succcessfully processed - "\
                     "Total Processing time: {5:.3f} secs. ({6:.3f} hours)"
               master.logger.log_info_p(msg.format(
                  master.name, job.get_branch_tree(), job.img_no, 
                  job.epoch,  file_type, elapsed_time, hours))

               master.logger.flush()

  
   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed. 
         @param master instance of the Master

         @code  
         Job results can be selectively obtained from: self.job_result_dico, e.g.: 
         print(self.job_result_dico.items()) # => all (key,value) pairs in the dico
         print(self.job_result_dico["control"].values()) # => all values in branch "/control"
         print(self.job_result_dico["control"]["ground"].values())  # => all values in "/control/ground"
         print(self.job_result_dico["control"]["ground"]["constant"].values()) # => all values in "/control/ground/constant"
         print(self.job_result_dico["control"]["ground"]["constant"][0]) # => job results for img_0 in "/control/ground/constant"g_0  
         @endcode
         @note Job results can be obtained with by accessing propertyies job_result_list and 
               job_result_dico
         @see Job, JobResult, MpfxJobResult
      """

      # --- Compute and generate statistics per branch

      selected_branches = self._get_selected_branches(master)
      if len(selected_branches) == 0:
         selected_branches = [None]

      for branch_tree in selected_branches:
         branch_job_results = self.get_job_results(branch=branch_tree)
         if len(branch_job_results) == 0:
            continue   

         if not branch_tree is None:
            branch_id = branch_tree.replace("/", "_")
         else:
            branch_id = ""

         if master.logging_enabled():
            if not branch_tree is None:
               master.logger.log_info_p(
                  "{0} - {1} - Generating statistics...".format(master, branch_id))
            else:
               master.logger.log_info_p(
                  "{0} - Generating statistics...".format(master))

         self._generate_result_statistics("result", 
                                          branch_id, branch_tree, branch_job_results, master)
         self._generate_result_statistics("error", 
                                          branch_id, branch_tree, branch_job_results, master)
         self._generate_result_statistics("filtered", 
                                          branch_id, branch_tree, branch_job_results, master)


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

   # -----------------------------------------------------------------------------------------------
   def _create_result_catalog(self, result_dico, error_dico, filter_dico, info_dico, job, master):
      """! Create a catalog with the columns: ID, e1, e2 in SE text format """            

      # --- Create a result catalog, marking failed ellipticity values
      
      failed_marking_value = master.config.get_as_float("FAILED_ELLIPTICITY_VALUE", 
                                                        "SHAPE_MEASUREMENT")

      is_ascii = True
      is_sextractor = False 
      if master.config.has_section("OUTPUT_CATALOGS.RESULT_CATALOGS"):
         is_ascii = master.config.get_as_boolean("IS_ASCII_FORMAT", 
                                                 "OUTPUT_CATALOGS.RESULT_CATALOGS")
         is_sextractor = master.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                      "OUTPUT_CATALOGS.RESULT_CATALOGS")

      if is_ascii:
         if job.branch_dirs != ".":
            output_filename = "results_{0:03d}-{1:1d}.txt".format(job.img_no, job.epoch)
         else:   
            output_filename = "results_{0}-{1:03d}-{2:1d}.txt".format(
                                                "_".join(job.branch_dirs), job.img_no, job.epoch)
      else:   
         if job.branch_dirs != ".":
            output_filename = "results_{0:03d}-{1:1d}.fits".format(job.img_no, job.epoch)
         else:   
            output_filename = "results_{0}-{1:03d}-{2:1d}.fits".format(
                                                "_".join(job.branch_dirs), job.img_no, job.epoch)
            
      output_directory = os.path.join(master.result_output_dir, job.get_branch_tree())

      col_list = ["Gal_id", "e1", "e2", "weight", "GAL_Xc", "GAL_Yc", "GAL_RA", "GAL_DEC"]
      col_key_map = {"GAL_id":0, "e1":1, "e2":2, "weight":3,\
                     "GAL_Xc":4, "GAL_Yc":5, "GAL_RA":6, "GAL_DEC":7}
      col_fmt_map = {"GAL_id":"%d", "e1":"%15.15f","e2":"%15.15f", "weight":"%4.3f",\
                     "GAL_Xc":"%9.6f", "GAL_Yc":"%9.6f", "GAL_RA":"%15.15f", "GAL_DEC":"%15.15f"}

      sub_result_dico = deepcopy(result_dico)   # Note: a deep copy is required here, otherwise  
                                                #       result_dico gets updated as well
            
      nb_unprocessed_ids = len(info_dico["unprocessed_ids"])

      sub_result_dico["GAL_id"].extend(info_dico["unprocessed_ids"])
      sub_result_dico["e1"].extend(numpy.repeat(failed_marking_value, nb_unprocessed_ids))      
      sub_result_dico["e2"].extend(numpy.repeat(failed_marking_value, nb_unprocessed_ids))      
      sub_result_dico["weight"].extend(numpy.repeat(0.0, nb_unprocessed_ids))      

      if "GAL_Xc" in sub_result_dico and "GAL_Yc" in sub_result_dico: 
         sub_result_dico["GAL_Xc"].extend(info_dico["unprocessed_GAL_Xc"])
         sub_result_dico["GAL_Yc"].extend(info_dico["unprocessed_GAL_Yc"])
      if "GAL_RA" in sub_result_dico and"GAL_DEC" in sub_result_dico: 
         sub_result_dico["GAL_RA"].extend(info_dico["unprocessed_GAL_RA"])
         sub_result_dico["GAL_DEC"].extend(info_dico["unprocessed_GAL_DEC"])

#       self.helper.save_from_list_dico(sub_result_dico, output_directory, output_filename, col_list,
#                                       key_index_map=col_key_map, key_fmt_map=col_fmt_map)

      self.helper.create_from_list_dico(sub_result_dico,output_directory, output_filename, \
                                       job, master, col_list, 
                                       key_index_map=col_key_map, key_fmt_map=col_fmt_map, 
                                       default_fmt="%.9e",
                                       is_ascii=is_ascii, is_sextractor=is_sextractor, hdu_no=1)

      # --- Create a unfiltered catalog if required
      if len(filter_dico["GAL_id"]) > 0:

         # --- Also Create a result catalog, without "failed" values but keeping all "filtered" values
         output_filename = "results_unfiltered_{0}-{1:03}.txt".format("_".join(job.branch_dirs), job.img_no)

          # --- A deep copy is required here, otherwise result_dico gets updated as well
         result_unfiltered_dico = deepcopy(sub_result_dico)

         if failed_marking_value in sub_result_dico["e1"] or failed_marking_value in sub_result_dico["e2"]:

            for (ID, e1) in zip(filter_dico["GAL_id"], filter_dico["e1"]):
               result_unfiltered_dico["e1"][sub_result_dico["GAL_id"].index(ID)] = e1       
      
            for (ID, e2) in zip(filter_dico["GAL_id"], filter_dico["e2"]):
               result_unfiltered_dico["e2"][sub_result_dico["GAL_id"].index(ID)] = e2 

         # --- Create the result "unfiltered" catalog
         self.helper.save_from_list_dico(result_unfiltered_dico, 
                                         output_directory, output_filename, col_list,
                                         key_index_map=col_key_map, key_fmt_map=col_fmt_map)

   # -----------------------------------------------------------------------------------------------
   def _create_data_catalog(self, result_dico, error_dico, filter_dico, info_dico, model_dico,
                                  file_type, job, worker):
      """! Create a catalog with the columns: ID, e1, e2 in SE text format """            

      if worker.config.get_as_boolean("DUMP_COLLECTED_DATA", "DATA_COLLECTION"):

         # --- Create a catalog that includes all collected data for this job
         output_result_data_directory = os.path.join(worker.result_output_dir, 
                                                     job.get_branch_tree())
         output_error_data_directory  = os.path.join(worker.error_output_dir, 
                                                     job.get_branch_tree())

         data_file_prefix_main, data_file_prefix_main_ext = os.path.splitext(
                                                                        file_type)
         
         is_ascii = True
         is_sextractor = False 
         if worker.config.has_section("OUTPUT_CATALOGS.DATA_CATALOGS"):
            is_ascii = worker.config.get_as_boolean("IS_ASCII_FORMAT", 
                                                   "OUTPUT_CATALOGS.DATA_CATALOGS")
            is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                         "OUTPUT_CATALOGS.DATA_CATALOGS")

         if is_ascii:
            file_ext = ".txt"
         else:      
            file_ext = ".fits"
         
         output_result_data_filename = "data_{0}-{1:03}-{2:1d}{3}".format(
                                       data_file_prefix_main, job.img_no, job.epoch, file_ext)
         output_error_data_filename = "errors_{0}_{1:03}-{2:1d}{3}".format(
                                  "_".join(job.branch_dirs), job.img_no, job.epoch, file_ext)
         output_filtered_data_filename = "filtered_{0}-{1:03}-{2:1d}{3}".format(
                                  "_".join(job.branch_dirs), job.img_no, job.epoch, file_ext)
         
         self._dump_all_collected_data(result_dico, 
                                       error_dico, filter_dico, model_dico, 
                                       output_result_data_directory, 
                                       output_result_data_filename,  
                                       is_ascii, is_sextractor,
                                       job, worker)

         if len(error_dico["GAL_id"]) > 0:
            self._dump_all_collected_data(error_dico, error_dico, 
                                          filter_dico, model_dico,
                                          output_error_data_directory, 
                                          output_error_data_filename,
                                          is_ascii, is_sextractor,
                                          job, worker)

         if len(filter_dico["GAL_id"]) > 0:
            self._dump_all_collected_data(filter_dico, 
                                          error_dico, filter_dico, model_dico,
                                          output_error_data_directory, 
                                          output_filtered_data_filename,
                                          is_ascii, is_sextractor,
                                          job, worker)

   # -----------------------------------------------------------------------------------------------
   def _generate_result_statistics(self, data_tag, branch_id, branch_tree, branch_results, master):

      # --- Where to store the statistics file
      if not branch_tree is None:
         stats_output_dir = os.path.join(master.stat_output_dir, branch_tree)
      else:
         stats_output_dir = master.stat_output_dir
      
      if len(branch_id) > 0:
         output_filename = "stats_{0}_{1}.txt".format(data_tag, branch_id)
      else:             
         output_filename = "stats_{0}.txt".format(data_tag)
         
      output_filepath = os.path.join(stats_output_dir, output_filename) 

      # --- Variables we do not want statistics on 
      except_list = ["GAL_id", "flag", "weight",\
                     #"xc", "yc",\
                     "gal_Xc", "gal_Yc", "psf_Xc", "psf_Yc",\
                     "GAL_Xc", "GAL_Yc",\
                     "GAL_RA","GAL_DEC"
                    ]

      with open(output_filepath, "a") as fd:

         title = "SUMMARY STATISTICS - Branch: {0}\n".format(branch_tree)
         fd.write(title)
         fd.write("=" * len(title) + "\n")         

         # --- Global stats   
         global_data_dico = {}
         global_data_dico[data_tag] = {}

         var_list = []

         # --- Compute statistics for the branch and generate the corresponding files 
         job_results = sorted(branch_results, key=attrgetter('job.img_no'))
         for job_result in job_results:

            job = job_result.job
            object_per_type_dico = job_result.result

            for file_type in object_per_type_dico.keys():

               if not file_type in object_per_type_dico[file_type]:
                  object_per_type_dico[file_type]["stats"] = {}

               # --- Variables we do want statistics on: model params, SE_*, extra fitting info
               data_dico  = object_per_type_dico[file_type][data_tag]

               # --- Renmoved failed fits    
               data_dico = deepcopy(data_dico)
               data_dico = self.helper.remove_failed_fits_from_result_dico(data_dico, master)
               data_dico = self.helper.remove_nans_from_result_dico(data_dico, master)    

               model_dico = object_per_type_dico[file_type]["model"]
               var_list = model_dico["param_names"]

               image_label = '\n--------------- {0} - {1:03d} ---------------\n\n'.format(
                                                                            file_type, job.img_no)

               object_per_type_dico[file_type]["stats"][data_tag], var_list = \
                          self._compute_and_dump_stats(fd, 
                                                       var_list,  
                                                       data_dico,
                                                       data_tag, image_label, 
                                                       except_list, 
                                                       job, master)

               # --- Compute and dump global stats 
               if len(object_per_type_dico[file_type]["stats"][data_tag]) > 0:

                  for var in var_list:
                     if var in data_dico:
                        if not var in global_data_dico[data_tag]:
                           global_data_dico[data_tag][var] = []
                        global_data_dico[data_tag][var].extend(data_dico[var])
                  #print "###global_data_dico:", global_data_dico

         # --- Dump global stats
         image_label = '\n---------------------------- All ----------------------------\n\n'

         self._compute_and_dump_stats(fd, 
                                       var_list,  
                                       global_data_dico[data_tag],
                                       data_tag, image_label,
                                       except_list, 
                                       job, master)

         fd.write("\n=== EOF ===")
         
   
      return object_per_type_dico



   # -----------------------------------------------------------------------------------------------
   def _compute_and_dump_stats(self, fd, var_list, data_dico, 
                                     data_tag, image_label, except_list, job, master):

      stats_dico = {}
      stats_dico[data_tag] = {}

      # --- If there are data available...
      if "iobj" in data_dico and len(data_dico["iobj"]) > 0:

         # --- Which variables to compute statistics on
         se_var_list = [var for var in data_dico.keys() \
                        if not var in var_list \
                        and not var in except_list and var.startswith("SE_")]
         var_list.extend(se_var_list)
         var_list.extend([var for var in data_dico.keys() \
                          if not var in var_list \
                          and not var in except_list and not var in se_var_list])

         fd.write(image_label)

         # --- Compute summary stats, populate stats_dico
         stats_dico = self._compute_summary_statistics(data_dico, var_list, job, master)

         # --- Dump stats values
         for var in var_list:
            fd.write("{0:<20}: ".format(var))
            for qoper in sorted(stats_dico[var].keys()):
               fd.write("{0} = {1:+.6e}  ".format(qoper.replace("numpy.", ""), 
                                                  stats_dico[var][qoper]))
            fd.write("\n") 

      return stats_dico, var_list


   # -----------------------------------------------------------------------------------------------
   def _compute_summary_statistics(self, result_dico, var_list, job, master):
      """! 
         Compute summary statistics on input data 
         @param result_dico dictionary with data arrays on which to compute statistics 
         @param car_list list of keys in the data dictionary that represents the stat. variables 
         @param job an object of class GfitJob to process
         @param worker instance of the master or a worker process      
         @return dictionary with summary statistics
      """

      # --- Compute stats and update the input data dictionary
      stats_dico = {}

      # --- Statistical operations to apply 
      qoper_tuples = [("numpy.size", numpy.size), ("numpy.sum", numpy.sum), \
                      ("numpy.min", numpy.min),   ("numpy.max", numpy.max), \
                      ("numpy.mean", numpy.mean), ("numpy.median", numpy.median), \
                      ("numpy.std", numpy.std)]

      # --- Compute statistics and create a dictionary with keys variable and sub-keys operation
      for qvar in var_list:
         stats_dico[qvar] = {}
         for (qlabel, qoper) in qoper_tuples:
            stats_dico[qvar][qlabel] = eval("qoper(result_dico[qvar])")

      return stats_dico

   # -----------------------------------------------------------------------------------------------
   def _dump_all_collected_data(self, result_dico, error_dico, filter_dico, model_dico, 
                                      output_directory, output_filename,
                                      is_ascii, is_sextractor, 
                                      job, master):
      """! Create a catalog containing all data collected during shape measurement for this job """            


      # --- Build column list

      # Essential parameters
      param_names = model_dico["param_names"]
      col_list = ["iobj","GAL_id", "flag", "weight", "x0", "y0"]
      col_list.extend(param_names)
      col_list.extend(["gal_Xc", "gal_Yc", "psf_Xc", "psf_Yc",\
                       "GAL_Xc", "GAL_Yc",\
                     ])
      col_key_map = dict(zip(col_list, list(numpy.arange(0, len(col_list)))))     
      col_fmt_map = {"iobj":"%d","GAL_id":"%d","flag":"%d", "weight":"%.2f",\
                     "x0":"%9.3f", "y0":"%9.3f", "gal_Xc":"%9.3f", "gal_Yc":"%9.3f",\
                     "psf_Xc":"%9.3f", "psf_Yc":"%9.3f",\
                     "GAL_Xc":"%.18e", "GAL_Yc":"%.18e",\
                    }
      col_fmt_map.update( dict( zip(param_names, ["%15.9f"] * len(param_names)) ) )

      # Extra parameters (non SExtractor)
      for var in result_dico.keys():
         if var in col_list or "centroids" in var or "SE_" in var: 
            continue

         col_list.append(var)
         col_fmt_map[var] = "%15.9f"

      # SExtractor variables
      for var in result_dico.keys():
         if "centroids" in var: 
            continue

         if "SE_" in var:
            col_list.append(var)
            col_fmt_map[var] = "%15.9f"
   
      #for var in result_dico.keys():
      #   print var, result_dico[var], len(result_dico[var])

#      # Remove records with failed fits
#      result_dico_copy = deepcopy(result_dico)
#      result_dico_copy = self.helper._remove_failed_fits_from_result_dico(result_dico_copy, master)

      # --- Also record ellipticities of failed and filtered fits  
      result_dico_copy = deepcopy(result_dico)

      failed_marking_value = master.config.get_as_float("FAILED_ELLIPTICITY_VALUE", 
                                                        "SHAPE_MEASUREMENT")
      if failed_marking_value in result_dico_copy["e1"] or \
         failed_marking_value in result_dico_copy["e2"]:

         # --- Get filtered ellipticity values
         for (ID, e1) in zip(filter_dico["GAL_id"], filter_dico["e1"]):
            result_dico_copy["e1"][result_dico_copy["GAL_id"].index(ID)] = e1       
         for (ID, e2) in zip(filter_dico["GAL_id"], filter_dico["e2"]):
            result_dico_copy["e2"][result_dico_copy["GAL_id"].index(ID)] = e2 

         # --- Get failed ellipticity values
         for (ID, e1) in zip(error_dico["GAL_id"], error_dico["e1"]):
            result_dico_copy["e1"][result_dico_copy["GAL_id"].index(ID)] = e1       
         for (ID, e2) in zip(error_dico["GAL_id"], error_dico["e2"]):
            result_dico_copy["e2"][result_dico_copy["GAL_id"].index(ID)] = e2 

      # --- Save records to file
#       self.helper.save_from_list_dico(result_dico_copy, output_directory, output_filename, col_list,
#                                       key_index_map=col_key_map, key_fmt_map=col_fmt_map)
      self.helper.create_from_list_dico(result_dico, output_directory, output_filename, \
                                       job, master, col_list, 
                                       key_index_map=col_key_map, key_fmt_map=col_fmt_map, 
                                       default_fmt="%.9e",
                                       is_ascii=is_ascii, is_sextractor=is_sextractor, hdu_no=1)


   # -----------------------------------------------------------------------------------------------
   def _create_galaxy_SE_dataset(self, master):
      """! Create a MpfxDataset object to query Sextractor-related data """

      dataset_name = master.config.get_as_string("NAME", "SEXTRACTOR_DATASET_GALAXY")
      dataset_dir  = master.config.get_as_string("BASE_DIR", "SEXTRACTOR_DATASET_GALAXY")

      return MpfxDataset(master, dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _create_psf_SE_dataset(self, master):
      """! Create a MpfxDataset object to query Sextractor-related data """

      dataset_name = master.config.get_as_string("NAME", 
                                                 "SEXTRACTOR_DATASET_PSF")
      dataset_dir  = master.config.get_as_string("BASE_DIR",
                                                 "SEXTRACTOR_DATASET_PSF")

      return MpfxDataset(master, dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _create_psf_image_dataset(self, master):
      """! Create a MpfxDataset object to query PSF-related data at the positions of galaxies """

      dataset_name = master.config.get_as_string("NAME", "PSF_DATASET")
      dataset_dir  = master.config.get_as_string("BASE_DIR", "PSF_DATASET")

      return MpfxDataset(master, dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _create_galaxy_sky_dataset(self, master):
      """! Create a MpfxDataset object to get the sky background image of galaxies 
           (if specified) 
      """

      dataset_obj = None
      
      dataset_name = "gal_sky_dataset"
      dataset_dir = "." 

      if master.config.get_as_boolean("SUBSTRACT_GALAXY_SKY", "SHAPE_MEASUREMENT"):
         if master.config.has_section("SKY_MODEL_GALAXY"):
         
            sky_gal_model_name =  master.config.get_as_string("SKY_MODEL_NAME", "SKY_MODEL_GALAXY")
            if sky_gal_model_name == "EXT_IMAGE":
               # --- A dataset for the sky image has been specified   
               section_name = "SKY_MODEL_GALAXY.{0}".format(sky_gal_model_name)
               section_dico = master.config.get_section_data(section_name)
               if "NAME" in section_dico:
                  dataset_name = master.config.get_as_string("NAME", section_name)
               if "BASE_DIR" in section_dico:
                  dataset_dir = master.config.get_as_string("BASE_DIR", section_name)
               dataset_obj = MpfxDataset(master, dataset_name, dataset_dir)   

      return dataset_obj

   # -----------------------------------------------------------------------------------------------
   def _create_psf_sky_dataset(self, master):
      """! Create a MpfxDataset object to get the sky background image of PSFs (if specified) """

      dataset_obj = None
      
      dataset_name = "psf_sky_dataset"
      dataset_dir = "." 
      if master.config.get_as_boolean("SUBSTRACT_PSF_SKY", "SHAPE_MEASUREMENT"):

         if master.config.has_section("SKY_MODEL_PSF"):
            sky_psf_model_name =  master.config.get_as_string("SKY_MODEL_NAME", "SKY_MODEL_PSF")
            if sky_psf_model_name == "EXT_IMAGE":
               # --- A dataset for the sky image has been specified   
               section_name = "SKY_MODEL_PSF.{0}".format(sky_psf_model_name)
               section_dico = master.config.get_section_data(section_name)
               if "NAME" in section_dico:
                  dataset_name = master.config.get_as_string("NAME", section_name)
               if "BASE_DIR" in section_dico:
                  dataset_dir = master.config.get_as_string("BASE_DIR", section_name)
               dataset_obj = MpfxDataset(master, dataset_name, dataset_dir)   

      return dataset_obj

   # -----------------------------------------------------------------------------------------------
   def _create_galaxy_noise_dataset(self, master):
      """! Create a MpfxDataset object to get the rms map of the noise in galaxy images 
           (if specified) 
      """

      dataset_obj = None
      
      dataset_name = "gal_noise_dataset"
      dataset_dir = "." 

      if master.config.has_section("GALAXY_NOISE"):
         
            gal_noise_model_name =  master.config.get_as_string("NOISE_MODEL_NAME", "GALAXY_NOISE")
            if gal_noise_model_name == "EXT_NOISE_MAP":
               # --- A dataset for the rms noise map has been specified   
               section_name = "GALAXY_NOISE.{0}".format(gal_noise_model_name)
               section_dico = master.config.get_section_data(section_name)
               if "NAME" in section_dico:
                  dataset_name = master.config.get_as_string("NAME", section_name)
               if "BASE_DIR" in section_dico:
                  dataset_dir = master.config.get_as_string("BASE_DIR", section_name)
               dataset_obj = MpfxDataset(master, dataset_name, dataset_dir)   

      return dataset_obj

   # -----------------------------------------------------------------------------------------------
   def _create_psf_noise_dataset(self, master):
      """! Create a MpfxDataset object to get the rms map of the noise in PSF images 
           (if specified) 
      """

      dataset_obj = None
      
      dataset_name = "psf_noise_dataset"
      dataset_dir = "." 

      if master.config.has_section("PSF_NOISE"):
         
            psf_noise_model_name =  master.config.get_as_string("NOISE_MODEL_NAME", "PSF_NOISE")
            if psf_noise_model_name == "EXT_NOISE_MAP":
               # --- A dataset for the rms noise map has been specified   
               section_name = "PSF_NOISE.{0}".format(psf_noise_model_name)
               section_dico = master.config.get_section_data(section_name)
               if "NAME" in section_dico:
                  dataset_name = master.config.get_as_string("NAME", section_name)
               if "BASE_DIR" in section_dico:
                  dataset_dir = master.config.get_as_string("BASE_DIR", section_name)
               dataset_obj = MpfxDataset(master, dataset_name, dataset_dir)   

      return dataset_obj

   # -----------------------------------------------------------------------------------------------
   def _locate_galaxy_se_files(self, se_dataset, job, master):
      """! 
         @param se_dataset a secondary dataset object pointing to the SExtractor galaxy data
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching SExtractor files
      """

      pattern_list = master.config.get_as_list("FILE_PATTERNS",
                                               "SEXTRACTOR_DATASET_GALAXY.CATALOG_PROPERTIES")
      dir_list = []
      dir_recurse = False      

      if master.config.has_key("DIR_INPUT_LIST", "SEXTRACTOR_DATASET_GALAXY"):
         dir_list = master.config.get_as_list("DIR_INPUT_LIST", "SEXTRACTOR_DATASET_GALAXY")

      if master.config.has_key("DIR_RECURSE", "SEXTRACTOR_DATASET_GALAXY"):
         dir_recurse = master.config.get_as_list("DIR_RECURSE", "SEXTRACTOR_DATASET_GALAXY")
      return self._get_matching_files(se_dataset, pattern_list, dir_list, dir_recurse, job, master) 

   # -----------------------------------------------------------------------------------------------
   def _locate_psf_se_files(self, se_dataset, job, master):
      """! 
         @param se_dataset a secondary dataset object pointing to the SExtractor STAR data
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching SExtractor files
      """

      pattern_list = master.config.get_as_list("FILE_PATTERNS", 
                                               "SEXTRACTOR_DATASET_PSF.CATALOG_PROPERTIES")

      dir_list = []
      dir_recurse = False      

      if master.config.has_key("DIR_INPUT_LIST", "SEXTRACTOR_DATASET_PSF"):
         dir_list = master.config.get_as_list("DIR_INPUT_LIST", "SEXTRACTOR_DATASET_PSF")
      if master.config.has_key("DIR_RECURSE", "SEXTRACTOR_DATASET_PSF"):
         dir_recurse = master.config.get_as_list("DIR_RECURSE", "SEXTRACTOR_DATASET_PSF")
      return self._get_matching_files(se_dataset, pattern_list, dir_list, dir_recurse, job, master) 

   # -----------------------------------------------------------------------------------------------
   def _locate_psf_image_files(self, psf_dataset, job, master):
      """!
         @param psf_dataset a secondary dataset object pointing to the PSF image data
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching PSF files
      """

      pattern_list = master.config.get_as_list("FILE_PATTERNS", "PSF_DATASET.IMAGE_PROPERTIES")

      dir_list = []
      dir_recurse = False      

      if master.config.has_key("DIR_INPUT_LIST", "PSF_DATASET"):
         dir_list = master.config.get_as_list("DIR_INPUT_LIST", "PSF_DATASET")

      if master.config.has_key("DIR_RECURSE", "PSF_DATASET"):
         dir_recurse = master.config.get_as_list("DIR_RECURSE", "PSF_DATASET")
      return self._get_matching_files(psf_dataset, pattern_list, dir_list, dir_recurse, job, master) 

   # -----------------------------------------------------------------------------------------------
   def _locate_sky_image_files(self, galaxy_dataset, section_name, job, master):
      """!
         @param galaxy_dataset a secondary dataset object pointing to the galaxy image data
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      pattern_list = master.config.get_as_list("FILE_PATTERNS", section_name)

      dir_list = []
      dir_recurse = False      

      if master.config.has_key("DIR_INPUT_LIST", section_name):
         dir_list = master.config.get_as_list("DIR_INPUT_LIST", section_name)

      if master.config.has_key("DIR_RECURSE", "SKY_MODEL_GALAXY.EXT_IMAGE"):
         dir_recurse = master.config.get_as_list("DIR_RECURSE", "SKY_MODEL_GALAXY.EXT_IMAGE")
      return self._get_matching_files(galaxy_dataset, pattern_list, dir_list, dir_recurse, 
                                      job, master) 

   # -----------------------------------------------------------------------------------------------
   def _locate_noise_map_files(self, dataset, section_name, job, master):
      """!
         @param dataset a secondary dataset object pointing to the galaxy or PSF noise maps
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      pattern_list = master.config.get_as_list("FILE_PATTERNS", section_name)

      dir_list = []
      dir_recurse = False      
      
      if master.config.has_key("DIR_INPUT_LIST", section_name):
         dir_list = master.config.get_as_list("DIR_INPUT_LIST", section_name)

      if master.config.has_key("DIR_RECURSE", section_name):
         dir_recurse = master.config.get_as_list("DIR_RECURSE", section_name)
      return self._get_matching_files(dataset, pattern_list, dir_list, dir_recurse, 
                                      job, master) 


   # -----------------------------------------------------------------------------------------------
   def _get_matching_files(self, dataset, file_patterns, dir_list, dir_recurse, job, master):
      """!
         @param dataset a dataset object
         @param file_patterns a list of Unix-like file patterns, like ["image*.fits", "*.cat"]
         @param dir_list a list of specific directories to search under base directory
         @param recurse tell whether to walk down directories (default @c True)
         @param job the GfitJob object corresponding to the image to analyze
         @param master master process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      # --- Locate the files matching the underlying criteria of the job (img_no, etc.)
      matching_dico = dataset.query(master, file_patterns, 
                                            dir_list=dir_list,
                                            image_list=[job.img_no], 
                                            image_range=[], 
                                            sort=True, recurse=dir_recurse)    

      # --- Aggregate all paths to handle multiepoch files
      # TODO: see how to handle multiple exposure files (multiepoch). May need another dataset
      #       containing multi-exposure files processed using a another tool.
      filepath_list = [] 
      for file_dico in matching_dico.values():
         filepath_list.extend(file_dico.values())

      # --- Build a dictionary with keys: (image_type, img_no, epoch) -> Value: file_path    
      file_path_dico = {}    
      for file_path in filepath_list:
         [file_type, img_no, epoch] = job.get_file_components_from_path(file_path)
         if not  (img_no, epoch) in file_path_dico:
            file_path_dico[(img_no, epoch)] = []
         # Store filetype_main, not file_type   
         file_path_dico[(img_no, epoch)].append(file_path)  

      # --- Find the 

#      # Build a dictionary with keys: (image_type, img_no, epoch) -> Value: file_path
#      # Note: it is assumed the SE catalog names are identical to the .FITS objedct file name from
#      #       which the centroids were estimated, except that the extension is different 
#      #       (e.g. .txt or .cat) => so we only store in <file_path_dico> the object filename without 
#      #       extension as file_type key.     
#      file_path_dico = {}    
#      for file_path in filepath_list:
#         [file_type, img_no, epoch] = job.get_file_components_from_path(file_path)
#         [filetype_main, filetype_ext] = os.path.splitext(file_type)
#         if not  (filetype_main, img_no, epoch) in file_path_dico:
#            file_path_dico[(filetype_main, img_no, epoch)] = []
#         # Store filetype_main, not file_type   
#         file_path_dico[(filetype_main, img_no, epoch)].append(file_path)  

      #print "FILEPATH DICO:", file_path_dico

      return file_path_dico

   # -----------------------------------------------------------------------------------------------
   def _get_selected_branches(self, master):

      return master.config.get_as_list("DIR_INPUT_LIST", "PRIMARY_DATASET")      

# -------------------------------------------------------------------------------------------------
class GfitJob(MpfxJob):      

   """! 
       A Gfit job. Based on mpfg.Job and mpfx.MpfxJob.  
   """

   def __init__(self, master, dataset, img_path_dico, *args):
      """!   
         Construct a new job for identifiying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:
         @code self.img_path_dico[type]@endcode with one of the aforementioned types.

         @param master master process instance
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized GfitJob object
      """

      MpfxJob.__init__(self, master, dataset, img_path_dico, *args)

      # --- PSF and SExtractor File path dictionaries
      self._psf_filepath_dico = {}
      self._galaxy_se_filepath_dico = {}
      self._star_se_filepath_dico   = {}

      # --- External sky and noise image file path dictionaries
      self._psf_sky_filepath_dico    = {}
      self._galaxy_sky_filepath_dico = {}
      
      self._psf_noise_map_filepath_dico    = {}
      self._galaxy_noise_map_filepath_dico = {}

      # --- Keep track of all defined datasets
      self._dataset_dico = {}  

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   # --- Getters

   @property
   def galaxy_se_filepath_dico(self):
      """! 
         Get the absolute paths of SExtractor catalogs linked to the GREAT3 galaxy files referenced
         by this job 
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of SExtractor matching files
      """
      return self._galaxy_se_filepath_dico

   @property
   def star_se_filepath_dico(self):
      """! 
         Get the absolute paths of SExtractor catalogs linked to the GREAT3 PSF files referenced
         by this job 
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of SExtractor matching files
      """
      return self._star_se_filepath_dico

   @property
   def psf_filepath_dico(self):
      """! 
         Get a the absolute paths of PSF image absolute paths matching the galaxies referenced by 
                 this job
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of PSF matching files
      """
      return self._psf_filepath_dico

   @property
   def galaxy_sky_filepath_dico(self):
      """! 
         Get a the absolute paths of external sky image absolute paths matching the galaxies
         referenced by this job
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of the matching files
      """
      return self._galaxy_sky_filepath_dico

   @property
   def psf_sky_filepath_dico(self):
      """! 
         Get a the absolute paths of external sky image absolute paths matching the PSFs
         referenced by this job
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of the matching files
      """
      return self._psf_sky_filepath_dico

   @property
   def galaxy_noise_map_filepath_dico(self):
      """! 
         Get a the absolute paths of external noise map paths matching the galaxies
         referenced by this job
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of the matching files
      """
      return self._galaxy_noise_map_filepath_dico

   @property
   def psf_noise_map_filepath_dico(self):
      """! 
         Get a the absolute paths of external noise map absolute paths matching the PSFs
         referenced by this job
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of the matching files
      """
      return self._psf_noise_map_filepath_dico

   @property
   def dataset_dico(self):
      """! 
         Get the dataset dictionary containing all created datasets
         @return a dictionary containing all registered datasets
      """

#   @property
#   def galaxy_stamp_size(self):
#      """! 
#         Get the original PSF postage stamp size
#      """
#      return self._galaxy_stamp_size

#   @property
#   def psf_stamp_size(self):
#      """! 
#         Get the original PSF postage stamp size
#      """
#      return self._psf_stamp_size


   # --- Setters

   @galaxy_se_filepath_dico.setter
   def galaxy_se_filepath_dico(self, galaxy_se_filepath_dico):
      """! 
         Specify the absolute paths to SExtractor catalog linked to the GREAT3 galaxy files 
         referenced by this job
         @param galaxy_se_filepath_dico filepaths dictionary containing for each file_type, 
                image number and epoch the absolute path of SExtractor matching files
      """
      self._galaxy_se_filepath_dico = galaxy_se_filepath_dico

   @star_se_filepath_dico.setter
   def star_se_filepath_dico(self, star_se_filepath_dico):
      """! 
         Specify the absolute paths of SExtractor catalog linked to the PSF GREAT3 files referenced
         by this job
         @param star_se_filepath_dico dictionary containing for each file_type, image number
                and epoch the absolute path of SExtractor matching files
      """
      self._star_se_filepath_dico = star_se_filepath_dico

   @psf_filepath_dico.setter
   def psf_filepath_dico(self, psf_filepath_dico):
      """! 
         Specify the absolute paths of PSF image absolute paths matching the galaxies referenced by 
         this job
         @param psf_filepath_dico dictionary containing for each file_type, image number 
                and epoch the absolute path of PSF matching files
      """
      self._psf_filepath_dico = psf_filepath_dico

   @galaxy_sky_filepath_dico.setter
   def galaxy_sky_filepath_dico(self, galaxy_sky_filepath_dico):
      """! 
         Specify the absolute paths of the sky image absolute paths matching the galaxies 
         referenced by this job
         @param psf_sky_filepath_dico dictionary containing for each file_type, image number 
                and epoch the absolute path of matching files
      """
      self._galaxy_sky_filepath_dico = galaxy_sky_filepath_dico

   @psf_sky_filepath_dico.setter
   def psf_sky_filepath_dico(self, psf_sky_filepath_dico):
      """! 
         Specify the absolute paths of the sky image absolute paths matching the PSFs 
         referenced by this job
         @param psf_filepath_sky_dico dictionary containing for each file_type, image number 
                and epoch the absolute path of matching files
      """
      self._psf_sky_filepath_dico = psf_sky_filepath_dico

   @galaxy_noise_map_filepath_dico.setter
   def galaxy_noise_map_filepath_dico(self, galaxy_noise_map_filepath_dico):
      """! 
         Specify the absolute paths of the noise map absolute paths matching the galaxies 
         referenced by this job
         @param psf_noise_map_filepath_dico dictionary containing for each file_type, image number 
                and epoch the absolute path of matching files
      """
      self._galaxy_noise_map_filepath_dico = galaxy_noise_map_filepath_dico

   @psf_noise_map_filepath_dico.setter
   def psf_noise_map_filepath_dico(self, psf_noise_map_filepath_dico):
      """! 
         Specify the absolute paths of the noise map absolute paths matching the PSFs 
         referenced by this job
         @param psf_filepath_noise_map_dico dictionary containing for each file_type, image number 
                and epoch the absolute path of matching files
      """
      self._psf_noise_map_filepath_dico = psf_noise_map_filepath_dico

   @dataset_dico.setter
   def dataset_dico(self, dataset_dico):
      """! 
         Specify a dataset dictionary containing all created datasets
         @param dataset_dico dictionary containing all created datasets
      """
      self._dataset_dico = dataset_dico


# -------------------------------------------------------------------------------------------------
class GfitJobResult(MpfxJobResult):      

   """! 
       The result of a GREAT3 job. 
       @see mpfx.MpfxJobResult
       @see mpf.JobResult.
   """

   def __init__(self, job, result, worker):   
      """! Construct a GfitJobResult job result object """
      MpfxJobResult.__init__(self, worker, job, result)

# -- EOF gfit_job.py
