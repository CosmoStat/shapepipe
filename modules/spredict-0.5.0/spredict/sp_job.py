"""! 
   sp_job.py - Job Processing  
"""

# -- Python imports
import os, sys
import time
import numpy
from operator import itemgetter, attrgetter
from copy import deepcopy
import shutil

# -- External imports
from sconfig import SConfig      # configuration file management
from scatalog import *           # catalog management
from mpfg3.mpfg3_job import *    # base job processing for GREAT3
from mpfg3.mpfg3_data import *   # data access to GREAT3

# --- Module-specific imports
from sp_interp import *       # base interpolation classes
from sp_map import *          # mapper
from sp_plot import *         # plotter 
from sp_helper import *       # helper utility functions



# -------------------------------------------------------------------------------------------------
class SpredictJobProcessor(Mpfg3JobProcessor):
   
   """! 
      Job processor: submit jobs and process associated job results. Based on mpfg3.JobProcessor.   
   """

   def __init__(self, master):
      """! Job Processor constructor """

      Mpfg3JobProcessor.__init__(self, master)
      
      self._interpolator = SpatialInterpolator(master)         # interpolator (stateless)

      self._mapper  = SpredictMapper()                         # mapper (stateless)
      self._plotter = SpredictPlotter()                        # plotter (stateless)
      self._helper  = SpredictHelper()                         # helper (stateless)


   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def interpolator(self):
      """! @return the SpatialInterpolator instance. """
      return self._interpolator

   @property
   def mapper(self):
      """! @return the SpredictMapper instance. """
      return self._mapper

   @property
   def plotter(self):
      """! @return the SpredictPlotter instance. """
      return self._plotter

   @property
   def helper(self):
      """! @return the SpredictHelper instance. """
      return self._helper

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def create_jobs(self, master):
      """! 
         Locate all objects to process and create the corresponding jobs. 
         @param master Master object instance
         @return the list of created jobs 
      """

      # --- Name of the model whose parameters will be interpolated
      #     Note: the interpolation methods to use will depend on the parameters, as specified
      #           in file "interpolate.cfg" 
      interp_config = self._load_interp_config(master) 

      # --- GREAT3 Dataset
      dataset_name = master.config.get_as_string("NAME", "DATASET")
      dataset_base_dir = master.config.get_as_string("BASE_DIR", "DATASET")
      dataset_branches = master.config.get_as_list("BRANCHES", "DATASET")
      dataset_data_types = master.config.get_as_list("DATA_TYPES", "DATASET")
      dataset_obs_types  = master.config.get_as_list("OBS_TYPES", "DATASET")
      dataset_image_list = master.config.get_as_list("IMAGE_LIST", "DATASET")
      dataset_image_range = master.config.get_as_list("IMAGE_RANGE", "DATASET")

      dataset = Great3Dataset(dataset_name, dataset_base_dir)

      # --- Search for matching files to process

      if master.logging_enabled():

         # --- Dump interpolation configuration
         master.helper.dump_config(interp_config, master.logger, "INTERPOLATION ")

         master.logger.log_info_p("{0} - Looking for files to process in directory: {1}...".format(
                                                             master, dataset_base_dir))     
         master.logger.log_info_p("{0} - Branches  : {1}".format(master, dataset_branches))     
         master.logger.log_info_p("{0} - obs types : {1}".format(master, dataset_obs_types))     
         master.logger.log_info_p("{0} - data types: {1}".format(master, dataset_data_types))     

      dataset_pattern_list = master.config.get_as_list("FILE_PATTERNS", "DATASET")

      filepath_dico = dataset.query(dataset_pattern_list, 
                                    branches=dataset_branches,
                                    obs_types=dataset_obs_types,
                                    data_types=dataset_data_types,
                                    image_list=dataset_image_list,
                                    image_range=dataset_image_range,
                                    sort=True, recurse=True)


      # --- For the multiepoch branch, make sure the catalog paths are also included
      #     for images with non-zero epoch indice
      for branch_key in filepath_dico: 
         (branch, obs_type, data_type, img_no, epoch) = branch_key
         if len(filepath_dico[branch_key].keys()) < 2:
            # --- Extract catalog file_types and paths from image with 0 epoch index
            #     and insert them in filepath_dico[branch_key] where they are missing
            branch_items = filepath_dico[(branch, obs_type, data_type, img_no, 0)].items()
            for (file_type, file_path) in branch_items:
               if "catalog" in file_type:
                  filepath_dico[branch_key][file_type] = file_path

            #print "==>", filepath_dico[(branch, obs_type, data_type, img_no, epoch)]
            
      # --- Create job objects, each of them specifying:
      #     - paths to catalogs referenced by the job
      #     - interpolation configuration object, where, among other things are specified:
      #       - the name of the model whose parameters will be interpolated (e.g. "moffat")
      #       - an interpolation method (e.g. RBF) for interpolating one or more model parameters
      #         (e.g. Moffat ["e1", "e2"])
      self._job_list = self._create_jobs(master, interp_config, filepath_dico)

      # --- Create output directories for job outputs:
      #      (<output_dir>/<branch>/<obs_type>)/<data_type>/
      output_dir_dico = master.config.get_section_data("DIR.OUTPUT")
      for (dir_key, dir_value) in output_dir_dico.items():
         if dir_key.startswith("OUTPUT_"):
            if dir_key.find("LOG") == -1:  # no tree for log dirs
               base_dir = os.path.join(master.run_output_dir, dir_value)
               self._create_dir_tree(base_dir, dataset_branches, 
                                            dataset_obs_types, dataset_data_types)
      
      return sorted(self._job_list, key=attrgetter('img_no'))   # sort by img_no


   # -----------------------------------------------------------------------------------------------
   def create_job(self, master, model_name, method_name, param_name, *args):
      """! 
         Factory method for creating a SpredictJob object. 
         @param master the master process
         @param model_name name of the model whose parameters are to be interpolated
         @param method_name name of the interpolation methods to use
         @param param_name name of the parameter (or list of grouped parameters) to process
         @param args a list of arguments chosen by the caller and whose nature depends on the 
                     data necessary to process the job
         @return a new SpredictJob job object
         @note overrides Mpfg3JobProcessor.create_job() 
         @see Job, Mpfg3Job
      """

      # --- Parameter list
      if type(param_name) == list:
         param_names = param_name     # already a list
      else:
         param_names = [param_name]   # a one-parameter list

      # --- Create a Mpfg3Job instance, based on parent class mpf.Job
      job = SpredictJob(model_name, method_name, param_names, *args)  

      # --- Specify the names of the model and interpolation method required to process the job  
      job.model_name  = model_name
      job.method_name = method_name

      # --- Attach the relevant filepath dictionaries required to process the job
      job.input_catalog_filepath_dico  = self._locate_input_catalogs(job, master)
      job.target_catalog_filepath_dico = self._locate_target_catalogs(job, master)
      job.model_catalog_filepath_dico  = self._locate_model_catalogs(job, master)

      #print job.input_catalog_filepath_dico
      #print job.target_catalog_filepath_dico
      #print job.model_catalog_filepath_dico 

      return job  # here we create a SpredictJob instance, based on parent class mpfg3.Job
      
   # -----------------------------------------------------------------------------------------------
   def preprocess_job(self, job, worker):
      """!
         Invoked by the Worker to perform some optional pre-processing on the job. 
         It may be for instance some format conversion. or other preparation step before the actual
         processing takes place in process_job().

         @param job an object of class SpredictJob to process
         @param worker instance of the worker process

         @return an object of class SpredictJobResult containing the data of the processed job

         @note overrides Mpfg3JobProcessor.preprocess_job() 
               to perform any processing required <b>before</b> process_job() is called
         @see preprocess_job, postprocess_job, Job, JobResult, Mpfg3JobResult
      """
      pass  # we no nothing here in spredict

   # -----------------------------------------------------------------------------------------------
   def process_job(self, job, worker):
      """! 
         Process a job of class SpredictJob and return the corresponding results in the form of a 
         SpredictJobResult object to the Master.       

         @param job object of class SpredictJob with processed data
         @param worker Worker object instance

         @return an object of class SpredictJobResult containing the data of the processed job

         @note overrides Mpfg3JobProcessor.process_job() 
         @see Job, Mpfg3Job, JobResult, Mpfg3JobResult 
      """

      # --- We want to process catalogs that contain the input, target coordinates and
      #     input parameter values at input coordinates. The result of the processing are the 
      #     predicted values at target coordinates.

      object_per_type_dico = {}  # will contain all interpolated results per image type   

      try:

         # --- Load model
         model  = self._get_model(worker, job.model_name)

         # --- Dump model configuration
         if worker.logging_enabled():
            worker.helper.dump_config(model.config, worker.logger, "MODEL ")

         # --- For each file type referenced by the job...
         file_types = job.get_file_types()

         for file_type in file_types: 

            # --- We are interested in input (deep)star catalogs
            if "star_catalog" in file_type:  
               if worker.logging_enabled():
                  file_main, file_ext = os.path.splitext(file_type)
                  worker.logger.log_info_p(
                    "{0} - /{1}/{2}-{3:03d}.{4} - Interpolating parameter(s) {5}\t({6})...".format(
                       worker.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
                       job.param_names, job.method_name))

               try:

                  # --- Setup interpolation method
                  method = self._get_method(worker, job.method_name)

                  if worker.logging_enabled():

                        # --- Dump model configuration
                        worker.helper.dump_config(method.method_config, worker.logger, "METHOD ")

                        worker.logger.log_info(
                          "{0} - /{1}/{2}-{3:03d}.{4} - Reading catalogs...".format(
                             worker.name, job.get_branch_tree(), file_main, job.img_no, file_ext))

                  # --- Create a InterpolationRequest object to process by the Interpolator
                  interp_request = InterpolationRequest(file_type, model, method, self.helper, 
                                                        job, worker)

                  if worker.logging_enabled():
                     worker.logger.log_info(
                      "{0} - /{1}/{2}-{3:03d}.{4} - Applying {5} interpolation algorithm "\
                      "on parameter(s) {6}...".format(
                         worker.name, job.get_branch_tree(), file_main, 
                         job.img_no, file_ext, job.method_name, job.param_names))

                  # --- Process the InterpolationRequest object and produce interpolation results
                  object_per_type_dico[file_type] = self.interpolator.interpolate(
                                                                              interp_request,
                                                                              job, worker)

                  

                  if len(object_per_type_dico[file_type]) > 0:
                     if worker.logging_enabled():
                           worker.logger.log_info_p(
                             "{0} - /{1}/{2}-{3:03d}.{4} - Interpolation complete".format(
                              worker.name, job.get_branch_tree(), file_main, job.img_no, file_ext))

                     # --- Compute statistics
                     # TODO

                     # --- Create Input & Target plots
                     if worker.config.get_as_boolean("CREATE_PLOTS", "PLOTS.INPUT") or \
                        worker.config.get_as_boolean("CREATE_PLOTS", "PLOTS.TARGET"):        
                        self.plotter.create_plots(model, method, 
                                                  object_per_type_dico[file_type]["interpolated"], 
                                                  file_type, job, worker)                  

                     # --- Create Input & Target Maps
                     if worker.config.get_as_boolean("CREATE_MAPS", "MAPS.INPUT") or \
                        worker.config.get_as_boolean("CREATE_MAPS", "MAPS.TARGET"):        
                        self.mapper.create_maps(model, method, 
                                                object_per_type_dico[file_type]["interpolated"], 
                                                file_type, job, worker)

               except SpatialInterpolator.InterpolationError as detail:
                  if worker.logging_enabled():
                     worker.logger.log_error_p(detail)

                  continue    # nevertheless, proceed with the next job
      

      except Exception as detail:

         if worker.logging_enabled():
            worker.logger.log_error_p(
                        "{0} - Some error occurred while processing job: {1} ({2})".format(
                                                                                worker.name, 
                                                                                job,
                                                                                detail))

      # --- Create a SpredictJobResult object with the results from all image types
      return SpredictJobResult(job, object_per_type_dico)


   # -----------------------------------------------------------------------------------------------
   def postprocess_job(self, job_result, worker):
      """! 
         Invoked by the Worker to perform some optional post-processing on the job.
         @param job_result object of class SpredictJobResult with processed data
         @param worker instance of the worker process

         @note overrides Mpfg3JobProcessor.preprocess_job() 
               to perform any processing required <b>after</b> process_job() is called
         @see process_job, postprocess_job, JobResult, Mpfg3JobResult
      """
      pass  # we no nothing here in spredict

   # -----------------------------------------------------------------------------------------------
   def process_job_result(self, job_result, master):
      """! 
         Process the result associated with a processed job.       

         @param job_result object of class SpredictJobResult with processed data
         @param master Master object instance
         @note overrides Mpfg3JobProcessor.process_job_result()
         @see Job, JobResult, Mpfg3JobResult
      """
      
#      # --- Model, Method, Param(s)
#      job = job_result.job
#      model  = self._get_model(master, job.model_name)
#      method = self._get_method(master, job.method_name)

#      # --- Interpolation config section model-method section name
#      config_section_name = "{0}.{1}".format(model.name.upper(), method.name) 

#      # --- Iterate through the results of all file types in the job...  
#      object_per_type_dico = job_result.result

#      for file_type in object_per_type_dico.keys():

#         if master.logging_enabled():
#            file_main, file_ext = os.path.splitext(file_type)
#            master.logger.log_info(
#              "{0} - /{1}/{2}-{3:03d}.{4} - Processing job results for "\
#               "parameter(s) {5}\t({6})...".format(
#                 master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
#                 job.param_names, job.method_name))

#         result_dico = object_per_type_dico[file_type]

#         # --- Gives a status of what has been processed at job level

#         # --- Compute statistics
      
      pass


   # -----------------------------------------------------------------------------------------------
   def record_job_result(self, job_result, master):
      """! 
         Add one or more job results to the list of job results. This method is invoked by the
         Master. 
         @param master instance of the Master
         @param job_result object of class JobResult containing the data of the processed job  
         @see Job, JobResult       
      """

      branch = job_result.job.branch
      obs_type  = job_result.job.obs_type
      data_type = job_result.job.data_type
      img_no = job_result.job.img_no
      epoch = job_result.job.epoch

      if not branch in self.job_result_dico:
          self.job_result_dico[branch] = {}
      if not obs_type in self.job_result_dico[branch]:
          self.job_result_dico[branch][obs_type] = {}
      if not data_type in self.job_result_dico[branch][obs_type]:
          self.job_result_dico[branch][obs_type][data_type] = {}
      if not img_no in self.job_result_dico[branch][obs_type][data_type]:
          self.job_result_dico[branch][obs_type][data_type][img_no] = {}

      if not epoch in self.job_result_dico[branch][obs_type][data_type][img_no]:
         self.job_result_dico[branch][obs_type][data_type][img_no][epoch] = []

      # --- Record one or more job results with the relevant branch tree structure   
      #print "***SPREDICT: Recording job results:", (branch, obs_type, data_type, img_no), job_result.job.name
      self.job_result_dico[branch][obs_type][data_type][img_no][epoch].append(job_result)


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
         @see Job, JobResult, Mpfg3JobResult
      """

      # *** Note: all job rsults belonging to a given image-epoch are now grouped in a separate list    
      all_jobs_results = self.get_job_results()

      branches = master.config.get_as_list("BRANCHES", "DATASET")
      data_types = master.config.get_as_list("DATA_TYPES", "DATASET")
      obs_types  = master.config.get_as_list("OBS_TYPES", "DATASET")

      # --- Create result catalogs, regrouping jobs linked to the same image/catalog
      for branch in branches:  # control, etc.
         for obs_type in obs_types:  # ground vs. space
            for data_type in data_types:   # constant vs. variable

               # --- Get a list of all job results on a given image, i.e. all jobs results that
               #     relate to the model parameters interpolated for that image
               job_results_list = self.get_job_results(branch=branch,
                                                       obs_type=obs_type, data_type=data_type)

               # --- Do whatever process is required on the results after interpolation
               for job_results in job_results_list:  

                  if len(job_results) > 0:
                     
                     # --- Invoke model-specific method to possible alter interpolated results
                     ### self._adjust_results(job_results, master)

                     # --- Create a catalog with the G3 columns (x, y, ID, etc.) plus one column per
                     #     interpolated parameter data
                     self._create_result_catalog(job_results, master)

            # --- end for data_type
         # --- end for obs_type
      # --- end for branch

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _load_interp_config(self, master):

      # --- Config name is derived from name FittingElement          
      
      config_dir = master.arg_options.options['-d']
      config_filepath = os.path.join(config_dir, 'interpolate.cfg') 
      try:
         return SConfig(config_filepath)
      except:
         # No configuration file found => The user must manually provide the relevant information
         self.helper.print_error("error reading configuration file {0}.".format(config_filepath)) 
         return None

   # -----------------------------------------------------------------------------------------------
   def _create_jobs(self, master, interp_config, filepath_dico):

      # --- Parse the content of the interpolation configuration file and 
      #     - decide which model parameter(s) should be interpolated with which method
      #     - create the SpredictJob objects jobs and return them in a list

      job_list = []

      model_name = master.config.get_as_string("MODEL_NAME", "INTERPOLATION").lower()

      if interp_config.has_section(model_name.upper()):

         # --- For each method...
            
         for subsection_name in interp_config.get_subsections(model_name.upper()): 

            # --- Retrieve the list of model parameters specified for each interpolation  method
            [_, method_name] = subsection_name.split('.')
            param_names = interp_config.get_as_list("PARAMS", subsection_name)

            # --- create the jobs
            for param_name in param_names:            
               job_list.extend([self.create_job(master, model_name, method_name, param_name, 
                                                 filepath_dico[k], *k) \
                                 for k in sorted(filepath_dico.keys())] )
      
      else:
            self.helper.print_error("no configuration found in {0} for model {1}.".format(
                                                                                 "interpolate.cfg",
                                                                                 self.method.name))   

      return job_list


   # -----------------------------------------------------------------------------------------------
   def _import_module(self, module_name, module_dir):
      """! Find and load a module <module_name> from a directory <module_dir> """

      try:

         file_obj, filename, data = imp.find_module(module_name, [module_dir])
         return imp.load_module(module_name, file_obj, filename, data)

      except Exception as detail:
         self.helper.print_error("Could not load module {0} from {1} {2}".format(
                                                                  module_name, module_dir, detail))

   # -----------------------------------------------------------------------------------------------
   def _get_model(self, worker, model_name):
      """! Load and instantiate a Model object for method <model_name> """   

      model_obj = None 

      if model_name is not None:

         model_config_dir = os.path.join(worker.arg_options.options['-d'], "models")
         model_module_dir = os.path.join(worker.arg_options.options['-m'], "models")

         module = self._import_module(model_name, model_module_dir)
         if not module is None:
            model_obj = module.Model(config_dir=model_config_dir)
         else:
            self.helper.print_error("Could not setup model {0}".format(model_name))
      else:
         self.helper.print_error("Could not setup model {0}".format(model_name))

      return model_obj


   # -----------------------------------------------------------------------------------------------
   def _get_method(self, worker, method_name):
      """! Load and instantiate a Method object for method <method_name> """   

      interp_config_dir = worker.arg_options.options['-d']
      method_config_dir = os.path.join(worker.arg_options.options['-d'], "methods")
      method_module_dir = os.path.join(worker.arg_options.options['-m'], "methods")

      module = self._import_module(method_name, method_module_dir)

      return module.Method(method_config_dir=method_config_dir, 
                                 interp_config_dir=interp_config_dir)


   # -----------------------------------------------------------------------------------------------
   def _get_all_interpolated_param_values(self, job_results, master):

      # *** Note: by construction (see record_job_results()) all jobs results mentioned here 
      #           relate to the same branch and image. Nevertheless, we have multiple job 
      #           results, each containing the interpolated data of the parameters they manage.        

      # --- Find all file types and parameter names in job results
      file_types = []
      param_names = [] 
      method_names = []
      for job_result in job_results:
         file_types.extend(job_result.result.keys())      
         param_names.extend(job_result.job.param_names)
         method_names.append(job_result.job.method_name)

      file_types   = numpy.unique(file_types)
      param_names  = numpy.unique(param_names)

      # --- Dictionary that will contain the catalog data
      col_data_dico = {}
      for file_type in file_types:
         col_data_dico[file_type] = {}

      # --- Retrieve the interpolated data for each parameter(s) processed by the jobs  
      for file_type in file_types:
         for job_result in job_results:
            result_dico = job_result.result[file_type]

            for param_name in job_result.job.param_names:
               col_data_dico[file_type][param_name] = numpy.asarray(
                                         self._flatten_tiled_param_dico(result_dico["interpolated"], 
                                                                        param_name))

      return file_types, param_names, method_names, col_data_dico

#   # -----------------------------------------------------------------------------------------------
#   def _adjust_results(self, job_results, master):
#      """! 
#         Invoke model-specific method to possibly alter parameter results 
#         @param job_results list of SpredictJobResult objects for the same branch, inage, epoch
#         @param master process   
#      """            

#      # *** Note: by construction (see record_job_results()) all jobs results mentioned here 
#      #           relate to the same branch and image. Nevertheless, we have multiple job 
#      #           results, each containing the interpolated data of the parameters they manage.        

#      # --- Find all file types and parameter names in job results
#      file_types = []
#      param_names = [] 
#      method_names = []
#      for job_result in job_results:
#         file_types.extend(job_result.result.keys())      
#         param_names.extend(job_result.job.param_names)
#         method_names.append(job_result.job.method_name)

#      file_types   = numpy.unique(file_types)
#      param_names  = numpy.unique(param_names)

#      # --- Dictionary that will contain the catalog data
#      col_data_dico = {}
#      for file_type in file_types:
#         col_data_dico[file_type] = {}

#      model  = self._get_model(master, job_results[0].job.model_name)

#      # --- Retrieve the interpolated data for each parameter(s) processed by the jobs  
#      for file_type in file_types:
#         for job_result in job_results:
#            result_dico = job_result.result[file_type]

#            model.adjust_results(result_dico["interpolated"], master)


   # -----------------------------------------------------------------------------------------------
   def _create_result_catalog(self, job_results, master):
      """! 
         Create a catalog on disk with additional columns for interpolated data per parameter 
         @param job_results list of SpredictJobResult objects for the same branch, inage, epoch
         @param master process   
      """            

      # *** Note: by construction (see record_job_results()) all jobs results mentioned here 
      #           relate to the same branch and image. Nevertheless, we have multiple job 
      #           results, each containing the interpolated data of the parameters they manage.        

      # --- Get all interpolated parameter values (in <col_data_dico>) plus extra info
      file_types, param_names, method_names, col_data_dico =\
                                      self._get_all_interpolated_param_values(job_results, master)

      # --- Assign fixed positions and formats to some columns
      col_list = ["x", "y", "ID", "x_tile_index", "y_tile_index",\
                  "tile_x_pos_deg", "tile_y_pos_deg", "x_field_true_deg", "y_field_true_deg"]
      col_key_map = {"x":0, "y":1, "ID":2, "x_tile_index":3, "y_tile_index":4,\
                     "tile_x_pos_deg":5, "tile_y_pos_deg":6,\
                     "x_field_true_deg":7, "y_field_true_deg":8}
      col_fmt_map = {}

      # --- Retrieve the data related to target positions, IDs, etc. in the target catalog, and
      #     create the result catalog

      for file_type in file_types:

         # Just take the first job in the job_results list, to get the branch & image_no
         job = job_results[0].job   

         # --- Filepath of result catalog to create 
         file_main, file_ext = os.path.splitext(file_type)
         output_filename = "{0}-{1:03}.txt".format(file_main, job.img_no)
         output_directory = os.path.join(master.result_output_dir, job.get_branch_tree())
         new_catalog_filepath = os.path.join(output_directory, output_filename)

         if master.logging_enabled():
            file_main, file_ext = os.path.splitext(file_type)
            master.logger.log_info_p(
              "{0} - /{1}/{2}-{3:03d}.{4} - Generating result catalog: {5} "\
               "parameter(s) {6}\t({7})...".format(
                 master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
                 new_catalog_filepath, param_names, method_names))

         # --- Open catalog with target positions
         original_catalog_file_extension = master.config.get_as_string(
                                                                    "CATALOG_EXTENSION", 
                                                                    "CATALOG_PROPERTIES.TARGET")
         original_hdu_no = master.config.get_as_int("CATALOG_HDU_NO", 
                                                    "CATALOG_PROPERTIES.TARGET")
         original_catalog_is_SE = master.config.get_as_boolean("CATALOG_IS_SEXTRACTOR",
                                                               "CATALOG_PROPERTIES.TARGET") 
         original_catalog_filepath = self.helper.get_catalog_filepath(
                                                             job.target_catalog_filepath_dico, 
                                                             original_catalog_file_extension,
                                                             job)

         if original_catalog_file_extension == ".fits":
            catalog = FITSCatalog(original_catalog_filepath, hdu_no=original_hdu_no, mem_map=True)
         elif original_catalog_file_extension == ".txt":
            if original_catalog_is_SE:
               catalog = SExCatalog(original_catalog_filepath)
            else:
               catalog = TextCatalog(original_catalog_filepath)

         # --- Read target catalog column data
         catalog.open()

         for col_name in catalog.get_col_names():
            col_data_dico[file_type][col_name] = catalog.get_named_col_data(col_name)

         catalog.close()

         # --- Create the result catalog as a plain .txt file

         col_list.extend(param_names)
         col_key_map = dict(zip(col_list, list(numpy.arange(0, len(col_list)))))
         col_fmt_map = {}

         self.helper.save_from_list_dico(col_data_dico[file_type], 
                                         output_directory, output_filename, col_list,
                                         key_index_map=col_key_map, key_fmt_map=col_fmt_map)
      
#   # -----------------------------------------------------------------------------------------------
#   def _flatten_tiled_dico(self, tiled_dico):

##      data_array = numpy.asarray([])
#      data_list = []

#      tile_tuples = sorted(tiled_dico.keys()) # sorted by tile tuple (0,0), (0,1)...(1,0), (1.1)...

#      for (x_tile_no, y_tile_no) in tile_tuples:
#         data_list.extend(tiled_dico[(x_tile_no, y_tile_no)])
#         #data_array = numpy.concatenate( (data_array, tiled_dico[(x_tile_no, y_tile_no)]) )

#      return data_list

   # -----------------------------------------------------------------------------------------------
   def _flatten_tiled_param_dico(self, tiled_dico, param_name):

      data_list = []

      tile_tuples = sorted(tiled_dico.keys()) # sorted by tile tuple (0,0), (0,1)...(1,0), (1.1)...
      for (x_tile_no, y_tile_no) in tile_tuples:
         data_list.extend(
                       tiled_dico[(x_tile_no, y_tile_no)][param_name]["results"]["target_values"])

      return data_list


   # -----------------------------------------------------------------------------------------------
   def _create_star_input_dataset(self, config):
      """! Create a Great3Dataset object to query input star data """

      dataset_name = config.get_as_string("NAME", "INPUT_DATASET")
      dataset_dir  = config.get_as_string("BASE_DIR", "INPUT_DATASET")

      return Great3Dataset(dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _create_star_target_dataset(self, config):
      """! Create a Great3Dataset object to query target star data """

      dataset_name = config.get_as_string("NAME", "TARGET_DATASET")
      dataset_dir  = config.get_as_string("BASE_DIR", "TARGET_DATASET")

      return Great3Dataset(dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _create_model_param_dataset(self, config):
      """! Create a Great3Dataset object to query model-dependent parameter data """

      dataset_name = config.get_as_string("NAME", "MODEL_PARAM_DATASET")
      dataset_dir  = config.get_as_string("BASE_DIR", "MODEL_PARAM_DATASET")

      return Great3Dataset(dataset_name, dataset_dir)

   # -----------------------------------------------------------------------------------------------
   def _locate_input_catalogs(self, job, master):
      """! 
         @param dataset a GREAT3 dataset object pointing to an input catalog
         @param job the SpredictJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      star_dataset = self._create_star_input_dataset(master.config) 
      pattern_list = master.config.get_as_list("FILE_PATTERNS", "INPUT_DATASET")
      return self._get_matching_files(star_dataset, pattern_list, job) 

   # -----------------------------------------------------------------------------------------------
   def _locate_target_catalogs(self, job, master):
      """! 
         @param dataset a GREAT3 dataset object pointing to a target catalog
         @param job the SpredictJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      star_dataset = self._create_star_target_dataset(master.config) 
      pattern_list = master.config.get_as_list("FILE_PATTERNS", "TARGET_DATASET")
      return self._get_matching_files(star_dataset, pattern_list, job) 

   # -----------------------------------------------------------------------------------------------
   def _locate_model_catalogs(self, job, master):
      """!
         @param dataset a GREAT3 dataset object pointing to the catalog with model parameters
         @param job the SpredictJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching PSF files
      """

      model_dataset = self._create_model_param_dataset(master.config)
      pattern_list = master.config.get_as_list("FILE_PATTERNS", "MODEL_PARAM_DATASET")

      return self._get_matching_files(model_dataset, pattern_list, job) 


   # -----------------------------------------------------------------------------------------------
   def _get_matching_files(self, dataset, file_patterns, job):
      """!
         @param dataset a GREAT3 dataset object
         @param file_patterns a list of Unix-like file patterns, like ["image*.fits", "*.cat"]
         @param job the SpredictJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary containing for each file_type, image number and epoch the absolute 
         path of matching files
      """

      # --- Locate the files matching the underlying criteria of the job (branch, obs_type, etc.)
      matching_dico = dataset.query(file_patterns, branches=[job.branch], 
                                                   obs_types=[job.obs_type], 
                                                   data_types=[job.data_type], 
                                                   image_list=[job.img_no], 
                                                   image_range=[], 
                                                   sort=True, recurse=True)    

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

      return file_path_dico

   # -----------------------------------------------------------------------------------------------
   def _get_selected_branches(self, master):
      branches = master.config.get_as_list("BRANCHES", "DATASET")
      data_types = master.config.get_as_list("DATA_TYPES", "DATASET")
      obs_types  = master.config.get_as_list("OBS_TYPES", "DATASET")

      branch_key_list = []      

      for branch in branches:  # control, etc.
         for obs_type in obs_types:  # ground vs. space
            for data_type in data_types:   # constant vs. variable
               branch_key_list.append((branch, obs_type, data_type))

      return branch_key_list   

# -------------------------------------------------------------------------------------------------
class SpredictJob(Mpfg3Job):      

   """! 
       Represents a Spredict job. Based on mpf.Job and mpfg3.Mpfg3Job.  
   """

   def __init__(self, model_name, method_name, param_names,
                      img_path_dico, branch, obs_type, data_type, img_no, epoch=None):
      """! 
         Construct a job with a name of the form "branch_obs_type, data_type_img_xxx" and a 
         dictionary of paths for each file type in the dataset.
         As of Sept. 2013, the possible types are the following:

         For deep images  (training):
         - "deep_epoch_dither"
         - "deep_galaxy_catalog"
         - "deep_image"
         - "deep_star_catalog"
         - "deep_starfield_image"
         - "deep_subfield_offset"
         
         For regular images:
         - "epoch_dither"
         - "galaxy_catalog"
         - "image"
         - "star_catalog"
         - "starfield_image"
         - "subfield_offset"

         To obtain the path associated with a given type, just use:

            @code self.img_path_dico[type]@endcode with one of the aforementioned types.
      
         @param model_name name of the model to interpolate (like "moffat")
         @param method_name name of the interpolation method/algorithm to use
         @param param_names list of parameter names to interpolate as part of the job      
         @param img_path_dico file absolute path dictionary  for every file types referenced in 
                the job
         @param branch the GREAT3 dataset branch name (control, multiepoch...)
         @param obs_type the GREAT3 dataset observation type (ground, space)
         @param data_type the GREAT3 dataset data type (constant, variale)
         @param img_no the image number (000 - 199)
         @param epoch the epoch index (0 - 5)         

      """

      Mpfg3Job.__init__(self, img_path_dico, branch, obs_type, data_type, img_no, epoch)

      # --- names of model, method and parameters
      self._model_name  = model_name   # name of model
      self._method_name = method_name  # name of interpolation method
      self._param_names = param_names  # list of model parameters to interpolate

      # --- Update the name of the job
      param_str_list = []
      self.name = "{0}-{1}-{2}-{3}".format(model_name, method_name, "-".join(param_names), self.name)      

      # --- Filepath dictionaries
      self._input_catalog_filepath_dico  = {}  # filepaths of catalog with input coords
      self._target_catalog_filepath_dico = {}  # filepaths of catalog with target coords
      self._model_catalog_filepath_dico  = {}  # filepaths of catalog with inputs coords and params 
      

   def __str__(self):
      """! 
         String representation of a SpredictJob object 
         @return a Python String describing the SpredictJob object
      """
      return Mpfg3Job.__str__(self)


   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   # --- Getters

   @property
   def model_name(self):
      """! 
         Get the names of the model whose parameter(s) must be interpolated
         @return name of the model
      """
      return self._model_name

   @property
   def method_name(self):
      """! 
         Get the names of the method to use to interpolated the job model parameters
         @return name of the interpolation method
      """
      return self._method_name

   @property
   def param_names(self):
      """! 
         Get the name of the model parameter(s) to interpolate as part of the job
         @return name of the model parameters to interpolate
      """
      return self._param_names

   @property
   def input_catalog_filepath_dico(self):
      """! 
         Get the absolute paths of input star catalog file paths linked to the GREAT3 files 
         referenced by this job 
         @return a dictionary containing absolute paths of input star catalog file paths
      """
      return self._input_catalog_filepath_dico

   @property
   def target_catalog_filepath_dico(self):
      """! 
         Get the absolute paths of target star catalog file paths linked to the GREAT3 files 
         referenced by this job 
         @return a dictionary containing absolute paths of target star catalog file path
      """
      return self._target_catalog_filepath_dico

   @property
   def model_catalog_filepath_dico(self):
      """! 
         Get the absolute paths of model catalogs linked to the GREAT3 PSF files referenced
         by this job 
         @return a dictionary containing absolute paths of model catalog file path
      """
      return self._model_catalog_filepath_dico


   # --- Setters

   @model_name.setter
   def model_name(self, model_name):
      """! 
         Specify the name of the model whose parameters must be interpolated
         @param model_name name of the model whose parameters must be interpolated
      """
      self._model_name = model_name

   @method_name.setter
   def method_name(self, method_name):
      """! 
         Specify the name of the interpolation method to apply
         @param method_name name of the interpolation method to apply
      """
      self._method_name = method_name

   @param_names.setter
   def param_names(self, param_names):
      """! 
         Specify the list of model parameter names to interpolate
         @param param_names list of model parameter names to interpolate
      """
      self._param_names = param_names

   @input_catalog_filepath_dico.setter
   def input_catalog_filepath_dico(self, input_catalog_filepath_dico):
      """! 
         Specify the absolute paths of input star catalog file paths linked to the GREAT3 files 
         referenced by this job 
         @param input_catalog_filepath_dico dictionary containing absolute paths of input star catalog 
                                         file paths
      """
      self._input_catalog_filepath_dico = input_catalog_filepath_dico

   @target_catalog_filepath_dico.setter
   def target_catalog_filepath_dico(self, target_catalog_filepath_dico):
      """! 
         Specify the absolute paths of input star catalog file paths linked to the GREAT3 files 
         referenced by this job 
         @param target_catalog_filepath_dico dictionary containing absolute paths of input star catalog 
                                          file paths
      """
      self._target_catalog_filepath_dico = target_catalog_filepath_dico

   @model_catalog_filepath_dico.setter
   def model_catalog_filepath_dico(self, model_catalog_filepath_dico):
      """! 
         Specify the absolute paths of model catalogs linked to the GREAT3 PSF files referenced
         by this job 
         @param model_catalog_filepath_dico dictionary containing absolute paths of  model catalogs 
                                            linked to the GREAT3 PSF files referenced by this job
      """   
      self._model_catalog_filepath_dico = model_catalog_filepath_dico


# -------------------------------------------------------------------------------------------------
class SpredictJobResult(Mpfg3JobResult):      

   """! 
       Represents the result of a GREAT3 job. See mpfg3.Mpfg3JobResult and mpf.JobResult.
   """

   def __init__(self, job, result):   
      """! Construct a SpredictJobResult job result object """

      Mpfg3JobResult.__init__(self, job, result)


# -- EOF sp_job.py
