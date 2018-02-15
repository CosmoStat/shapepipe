"""! 
   mkp_job.py - Job Processing  
"""

# -- Python imports
import os, sys
import time
import numpy
from operator import itemgetter, attrgetter

# -- External imports
from mpfg3.mpfg3_job import *    # base job processing for GREAT3
from mpfg3.mpfg3_data import *   # data access to GREAT3

# --- Module-specific imports
from mkp_psf  import *           # PSF reconstruction        
from mkp_plot import *           # plotter 
from mkp_help import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class MkpJobProcessor(Mpfg3JobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpf.JobProcessor.   
   """

   def __init__(self, master):
      """! Job Processor constructor """

      Mpfg3JobProcessor.__init__(self, master)
      
      self._helper  = MkpHelper()         # helper utility functions
      self._plotter = MkpPlotter()        # plotter

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   # --- Getters

   @property
   def helper(self):
      """! @return the MkpHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the Plotter instance. """
      return self._plotter

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, dataset, *args):
      """! 
         Factory method for creating a MkpJob object. 
         @param master the master process
         @param dataset data source used to create the job
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                data necessary to process the job
         @return a new MkpJob job object
         @note overrides Mpfg3JobProcessor.create_job() 
         @see Job, Mpfg3Job
      """
      
      return MkpJob(master, dataset, *args)  # here we create a Mpfg3Job instance, based on parent class mpf.Job

   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """!
         Invoked by the Worker to perform some optional pre-processing on the job. 
         It may be for instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class MkpJob to process
         @param worker instance of the worker process

         @return an object of class MkpJobResult containing the data of the processed job

         @note overrides Mpfg3JobProcessor.preprocess_job() 
               to perform any processing required <b>before</b> process_job() is called
         @see preprocess_job, postprocess_job, Job, JobResult, Mpfg3JobResult
      """
      pass  # we no nothing here in mkp

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """! 
         Process a job of class MkpJob and return the corresponding results in the form of a 
         MkpJobResult object to the Master.       

         @param job object of class MkpJob with processed data
         @param worker Worker object instance

         @return an object of class MkpJobResult containing the data of the processed job

         @note overrides Mpfg3JobProcessor.process_job() 
         @see Job, Mpfg3Job, JobResult, Mpfg3JobResult 
      """

      # --- This is where we add the code to process the job

      # --- Depending on the object selection criteria, we may have to process either catalogs
      #     matching "*catalog-???.fits" or images matching "*image-???-??.fits"
      #     e.g. "star_catalog-002-0.fits" or "starfield_image-002-0.fits" for regular stars.
      #     The class attribute: job.img_path_dico contains all these paths as a dictionary
      #     (see description of class: MkpJob below).  

      # --- Here we only want to process images and catalogs. We don't any other files accesible
      #     from this job, like deep_images or other dither files. 

      object_per_type_dico = {}  # will contain the job results per image type (images)

      # --- Iterate over the path dictionary content to generate and process SEXtractor catalogs
      try:

         for file_type in job.get_file_types():    # check the type of file (looking for images) 

            image_filepath = job.get_file_path(file_type)
            if job.dataset.is_image(worker, image_filepath):   

               catalog_file_type = job.get_catalog_from_image(file_type)   # cat. with .txt ext.

               if worker.logging_enabled():
                  worker.logger.log_info_p(
                         "{0} - Processing job - /{1}/image-{2:03d}-{3:1d} - {4}...".format(
                            worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type))

               object_per_type_dico[file_type] = self._build_PSF_field(
                                                     file_type, catalog_file_type, job, worker)

            else:
               # Another file type (e.g. a catalog)              
               object_per_type_dico[file_type] = {}

      except:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                        "{0} - Some error occurred while processing job: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                sys.exc_info()[1]))

      # --- Create a MkpJobResult object with the results from all image types
      return MkpJobResult(worker, job, object_per_type_dico)


   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.
         @param job_result object of class MkpJobResult with processed data
         @param worker instance of the worker process

         @note overrides Mpfg3JobProcessor.preprocess_job() 
               to perform any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult, Mpfg3JobResult
      """
      pass  # we no nothing here in mkp

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job.       

         @param job_result object of class MkpJobResult with processed data
         @param master Master object instance
         @note overrides Mpfg3JobProcessor.process_job_result()
         @see Job, JobResult, Mpfg3JobResult
      """
      
      pass

   # -----------------------------------------------------------------------------------------------
   def all_jobs_processed(self, master):
      """! 
         This method is called by the Master once all the jobs have been processed. 
         @param master instance of the Master
         @note Job results can be selectively obtained by calling the get_job_results() methods,
               specifying the relevant query criterion: branch, observation type, data type, image 
         @note the entire list of Job results can aldo be obtained with by calling 
               JobProcessor.job_result_list()
         @see Job, JobResult, get_job_result()
      """

      # --- Create mosaic .FITS files from all results, for each image type, if requested
      if master.config.get_as_boolean("CREATE_PSF_MOSAIC", "DEBUGGING"):
         self._create_PSF_mosaic(master)
         

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   def _create_PSF_mosaic(self, master):
      """! 
         Create a mosaic FITS file containing one star selected from each of the 200 star fields. 
         The dimensions of the mosaic are 10 x 20 = 200.
         @param master master process
      """            

      # --- Regroup job results per branch tree and file types 
      all_jobs_results = self.get_job_results()
 
      type_per_branch_dico = {}
      for job_result in all_jobs_results:
         job = job_result.job
         branch_key = tuple(job.branch_dirs)
         if not branch_key in type_per_branch_dico:
            type_per_branch_dico[branch_key] = {}
         file_types = job_result.result.keys()
         for file_type in file_types:
            filepath = job.get_file_path(file_type)
            if filepath is None:
               continue
            if job.dataset.is_image(master, filepath):
               if not file_type in type_per_branch_dico[branch_key]:
                  type_per_branch_dico[branch_key][file_type] = self.get_job_results(
                                 branch=job.branch)

      # --- Create and populate the mosaic files for each branch/obs_type/data_type combinations
      #     in the appropriate output directories
      for branch_key in type_per_branch_dico.keys():
         for file_type in type_per_branch_dico[branch_key].keys():
            filepath = job.get_file_path(file_type)
            if filepath is None:
               continue
            if job.dataset.is_image(master, filepath):
 
               # --- Mosaic output filename and directory
               branch_tree = os.path.join(branch_key[0], os.path.join(branch_key[1], branch_key[2]))
               output_dir = os.path.join(master.result_output_dir, branch_tree)
               mosaic_filename = file_type.replace(
                                          "starfield_image", "psfmosaic_image_{0}_{1}_{2}".format(
                                                       branch_key[0], branch_key[1], branch_key[2]))     
               mosaic_filepath = os.path.join(output_dir, mosaic_filename)
 
               # --- Only keep job results with epoch index 0 and sort by img_no
               job_results = type_per_branch_dico[branch_key][file_type]
               job_results = sorted(job_results, key=attrgetter("job.img_no"))   # sorted / img_no
               job_results = [r for r in job_results \
                                if r.job.epoch is not None and r.job.epoch == 0 ]
 
               if len(job_results) > 0:
                  stamp_size = job_results[0].job.get_stamp_size()   # find stamp_size
 
                  # Populate mosaic
                  height = 20 * stamp_size
                  width  = 10 * stamp_size  
                  mosaic_data = numpy.zeros((height, width)).astype('float32')
     
                  row = 0
                  col = 0
                  for job_result in job_results:
                     if len(job_result.result.keys()) > 0:
                        # --- Populate the mosaic from the results
                        if file_type in job_result.result.keys():
                           mosaic_data[row:row+stamp_size, col:col+stamp_size] =\
                                                          job_result.result[file_type]["psf_image"]
                           col += stamp_size                 
                           if col >= width:
                              col = 0
                              row += stamp_size
                           if row >= height:
                              break 
 
                        #print (row, col)  
 
               # --- Create .FITS file
               self.helper.write_as_fits(mosaic_data, mosaic_filepath)
 
         # --- end for file
 
      # --- end for branch_key



   # -----------------------------------------------------------------------------------------------
   def _build_PSF_field(self, image_file_type, catalog_file_type, job, worker):

      # --- Load PSF recontruction module based on configuration info 
      #     *** Note: this cannot be done in __init__() because we don't have the Master's instance 
      #               to access the configutation object 

      build_method_name = worker.config.get_as_string("BUILD_METHOD", "METHODS")
      base_method_dir = worker.config.get_as_string("BASE_METHOD_DIR", "METHODS")
      module = self.helper.import_module(build_method_name, base_method_dir)      
      build_method = module.PSFBuildMethod(build_method_name, base_method_dir)

      # --- Invoke the PSFBuildMethod.make_PSF() method
      return build_method.make_PSF_field(image_file_type, catalog_file_type, job, worker)      
      
# -------------------------------------------------------------------------------------------------
class MkpJob(Mpfg3Job):      

   """! 
       Represents a Mkp job. Based on mpf.Job and mpfg3.Mpfg3Job.  
   """

       
   def __init__(self, master, dataset, img_path_dico, *args):

      """! 
         Construct a new job for identifying a set of files in the dataset to process together.

         To obtain the path associated with a given type, just use:

         @code self.img_path_dico[type]@endcode with one of the aforementioned types.
         @param master instance of the Master
         @param dataset data source used to create the job
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param args a list with the branch tree, image number and epoch number

         @return an initialized MkpJob object 
      """

      # --- Construct a new Job
      MpfxJob.__init__(self, master, dataset, img_path_dico, *args)

# -------------------------------------------------------------------------------------------------
class MkpJobResult(Mpfg3JobResult):      

   """! 
       Represents the result of a GREAT3 job. See mpfg3.Mpfg3JobResult and mpf.JobResult.
   """

   def __init__(self, worker, job, result):   
      """! Construct a Mpfg3JobResult job result object """

      MpfxJobResult.__init__(self, worker, job, result)

      # --- Cached Job result data and associated stats
#      self._data_dico  = self._collect_catalog_data()
#      self._stats_dico = self._compute_stats(self._data_dico)

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def data_dico(self):
      """! @return the data dictionary """
      return self._data_dico

   @property
   def stats_dico(self):
      """! @return the statistics dictionary """
      return self._stats_dico

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

#   # -----------------------------------------------------------------------------------------------
#   def _collect_catalog_data(self):
#      data_dico = {}

#      # --- Forbudden quantities
#      forbidden_col_names = ["NUMBER", "xc", "yc", "x", "y", "A_IMAGE", "B_IMAGE", "FLAGS", \
#                             "CLASS_STAR", "X_IMAGE", "Y_IMAGE", "MAGERR_BEST", "FLUXERR_BEST"]

#      # --- Collect data for each file type (image)
#      object_per_type_dico = self.result
#      for file_type in object_per_type_dico.keys():
#         se_run_dico = object_per_type_dico[file_type]["se_run_dico"]

#         # --- Open SE Catalog and collect column data
#         se_output_cat_filepath = se_run_dico["se_output_cat_filepath"]

#         se_catalog = SExCatalog(se_output_cat_filepath)
#         se_catalog.open()

#         # --- Extra, derived column names
#         data_dico[file_type] = {}
#         data_dico[file_type]["SNR"] =\
#                     se_catalog.get_col_data("FLUX_BEST") / se_catalog.get_col_data("FLUXERR_BEST")

#         # --- Allowed SExtractor column data
#         for col_name in se_catalog.get_col_names():
#            if not col_name in forbidden_col_names:
#               data_dico[file_type][col_name] = se_catalog.get_col_data(col_name)

#         se_catalog.close()

#      return data_dico

#   # -----------------------------------------------------------------------------------------------
#   def _compute_stats(self, data_dico):
#      """! Compute statistical operations on the SExtractor catalog data """

#      stats_dico = {}

#      # --- Statistics to calculate
#      stats_oper = [('Len','len'),('Sum','numpy.sum'),('Avg','numpy.mean'),('Med','numpy.median'),\
#                    ('Stdev','numpy.std'),('Min','numpy.min'),('Max','numpy.max')]
#      stats_oper_dict = dict(stats_oper)
#      stats_oper_keys = stats_oper_dict.keys()
#      stats_oper_keys.sort()

#      # --- Collect data for each file type (image)
#      object_per_type_dico = self.result
#      for file_type in object_per_type_dico.keys():

#         # --- Compute stats for variables recorded in dictionary
#         var_names = data_dico[file_type].keys()
#         var_names.sort()

#         stats_dico[file_type] = {}
#         for var_name in var_names:
#            stats_data = numpy.asarray(data_dico[file_type][var_name]) 
#            oper_result_dico = {}
#            for oper in stats_oper_keys:
#               o_ptr = eval(stats_oper_dict[oper])
#               oper_result_dico[oper] = o_ptr(stats_data)
#            stats_dico[file_type][var_name] = oper_result_dico

#      return stats_dico


# -- EOF mkp_job.py
