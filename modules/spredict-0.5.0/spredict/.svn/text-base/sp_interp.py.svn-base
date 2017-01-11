"""! 
   sp_interpolate.py - Base PSF interpolation module
"""

# -- Python imports
import os, sys
import string
import time
import numpy

# --- Module-specific imports
from sp_helper import *       # helper utility functions
from sp_plot import *         # plotter 
from sp_outlier import *      # outlier detection and removal
from sp_validate import *     # validation by resampling

# --- External imports
from sconfig import SConfig
from scatalog import *

# -------------------------------------------------------------------------------------------------
class ModelParam(object):
   """! Represents a parameter in the model """

   def __init__(self, name, bounds=[None, None]):
      """!
         Construct a ModelParam object
         @param name name of the model
         @param bounds parameter min and max bounds in the form of list [min, max]   
         @note the minimum, maximum values or both may be set to @c None 
      """

      self._name = name             # name of parameter
      self._bounds = bounds         # allowed range (<min value>, <max value>)

   def __str__(self):
      """! 
         String representation of a ModelParam object 
         @return a Python String describing the ModelParam object
      """
      return "[{0},{1}]".format(self._name, self._bounds)

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters

   @property
   def name(self):
      """! @return the name of the parameter """   
      return self._name

   @property
   def bounds(self):
      """! @return the bounds of the parameter """   
      return self._bounds

   # --- Setters

   @name.setter 
   def name(self, name):
      """!     
         Set the name of the model parameter 
         @param name parameter name
      """
      self._name = name     

   @bounds.setter 
   def bounds(self, bounds):
      """!     
         Set the bounds of the model parameter 
         @param bounds parameter min and max bounds in the form of list [min, max]   
         @note the minimum, maximum values or both may be set to @c None
      """
      self._bounds = bounds     



# -------------------------------------------------------------------------------------------------
class BaseModel(object):
   """! 
      Represents a model to which the data must be interpolated
      @note this is the base class for all interpolation models
   """

   def __init__(self, name, config_dir):
      """!# --------------------------------- Default Interpolation Options ----------------------------------
         Constuct BaseModel object
         @param name name of the model   
         @param config_dir configuration directory secific to the model
         @param interp_config_dir configuration directory combining both model and method-specific 
                settings
      """

      # --- Private attributes (should not be accessed directly)
      self.__params = []         # the model parameters
      self.__param_dico = {}     # to access parameter by names instead of indice

      # --- Model name
      self.name = name

      # --- Configuration objects   
      self._config = self._load_config(config_dir)  # configuration data (as a SConfig object))

      # --- Model Parameters list: set from configuration or otherwise, manually
      model_params = self._get_params_from_config(self._config)
      if len(model_params) > 0:
         self.set_params(model_params)        # set parameters based on configuration          

   def __str__(self):
      """! 
         String representation of a BaseModel object 
         @return a Python String describing the BaseModel object
      """
      descr = "Model:\n"
      descr += "- name: {0}\n".format(self.name)
      if len(self.params) > 0:
         descr += "- params: {0}\n".format([str(p) for p in self.params])

      return descr

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters  

   @property
   def name(self):
      """! @return the name of the model """   
      return self._name

   @property
   def params(self):
      """! @return the model parameters in the form of a list of ModelParam objects """
      return self.get_params()

   @property
   def config(self):
      """! @return the configuration object of the model """
      return self._config

   # --- Setters  

   @name.setter
   def name(self, name):
      """!  
         Set the model name
         @param name the model name
      """ 
      self._name = name

   @params.setter
   def params(self, params):
      """!  
         Set the model parameters in the form of a list of ModelParam objects 
         @param params a list of ModelParam object
      """ 
      self.set_params(params)


   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_params(self):
      """! @return the model parameters in the form of a list of ModelParam objects  
           (natural order """
      return self.__params

   # -----------------------------------------------------------------------------------------------
   def set_params(self, params):
      """! 
         Set the model parameters in the form of a list of ModelParam objects  
         @param params list of ModelParam objects
      """
      self.__params = params
      for p in params:
         self.__param_dico[p.name] = p
      self.__param_names = [p.name for p in params]

   # -----------------------------------------------------------------------------------------------
   def reset_params_from_config(self, config):
      """! Reset parameters from configuration data """

      model_params = self.get_params_from_config(config)
      if len(model_params) > 0:
         self.set_params(model_params)        # set parameters based on configuration     
 
   # -----------------------------------------------------------------------------------------------
   def get_param_names(self):
      """! @return the name of the model parameters (natural order) """
      return [p.name for p in self.__params]   

   # -----------------------------------------------------------------------------------------------
   def get_param_by_name(self, name):
      """! 
         Access a model parameter by name 
         @param name the name of the parameter# --------------------------------- Default Interpolation Options ----------------------------------
       """
      return self.__param_dico.get(name)
      
   # -----------------------------------------------------------------------------------------------
   def add_param(self, param):
      """!  
         Register a parameter in the form of a ModelParam object
         @param param a ModelParam object
      """
      self.__params.append(param)
      self.__param_dico[param.name] = param

   # -----------------------------------------------------------------------------------------------
   def del_param(self, param):
      """! 
         Register a parameter in the form of a Mo('re', [3.0, None], 0.0)delParam object 
         @param param a ModelParameter object
      """
      self.__params.remove(param)
      if param in self.__param_dico:
         del self.__param_dico[param.name]

   # -----------------------------------------------------------------------------------------------
   def get_bounds(self):
      """! @return the parameter bounds in the form of a list [min, max] """     
      return [p.bounds for p in self.__params]
            
   # -----------------------------------------------------------------------------------------------
   def _get_params_from_config(self, config):
      """! @return the list of model parameters from the model configuration, assuming a key with
          name @c "params" has been defined in a section with the same name as the model
      """
      params = []

      if not config is None: 

         if config.has_section(self.name.upper()):

            for subsection in config.get_subsections("{0}.PARAMS".format(self.name.upper())): 

               # --- Param name
               name = config.get_as_string("NAME", subsection)

               # --- Param value range
               bounds = config.get_as_list("RANGE", subsection)
               if not bounds is None:
                  bounds = bounds
               else:
                  bounds = None    

               params.append(ModelParam(name, bounds)) 

      return params


   # -----------------------------------------------------------------------------------------------
   def adjust_results(self, data_dico, master):
         
      """! 
         Optionally post-process interpolation results
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param result_dico dictionary containing interpolation results
         @param job a SpredictJob object
         @param worker Worker process object
         @note this method should be overridden by a derived model class in order to customize the
               adjustmment process
      """
      return data_dico


   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _load_config(self, config_dir):

      # --- Config name is derived from name FittingElement          
      config_filepath = os.path.join(config_dir, self.name) + '.cfg'   
      try:
         return SConfig(config_filepath)
      except:
         # No configuration file found => The user must manually provide the relevant information
         print("Spredict *** Warning ***: error reading configuration file {0}.".format(
                                                                               config_filepath)) 
         return None


# -------------------------------------------------------------------------------------------------
class BaseMethod(object):
   
   """! 
      spredict base interpolation method
   """

   def __init__(self, name, method_config_dir, interp_config_dir):
      """! Interpolator constructor """

      self._name = name                            # name of the interpolation method
      self._is_local = False                       # True if a local interpolator, False otherwise

      self._helper  = SpredictHelper()             # helper utility functions
      self._plotter = SpredictPlotter()            # plotter

      self._method_config = self._load_method_config(method_config_dir) # method-specific config
      self._interp_config = self._load_interp_config(interp_config_dir) 

   def __str__(self):
      """! 
         String representation of a BaseInterpolator object 
         @return a Python String describing the BaseInterpolator object
      """
      descr = "Interpolator:\n"
      descr += "- Name: {0}\n".format(self.name)
      descr += "- Local Interpolator: {0}\n".format(self.is_local)
      descr += "- {0}\n".format(self.method_config)
      
      return descr

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def name(self):
      """! @return the name of the model """   
      return self._name

   @property
   def is_local(self):
      """! @return True for a local interpolator, False otherwise """   
      return self._is_local

   @property
   def method_config(self):
      """! @return the method-specific configuration object """   
      return self._method_config

   @property
   def interp_config(self):
      """! @return the extra configuration object (combining model+method configuration) """   
      return self._interp_config

   @property
   def helper(self):
      """! @return the SpredictHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the SpredictPlotter instance. """
      return self._plotter

   # --- Setters  

   @name.setter
   def name(self, name):
      """!  
         Set the model name
         @param name the model name
      """
      self._name = name

   @is_local.setter
   def is_local(self, is_local):
      """!  
         Specify if the interpolator is local or global
         @param is_local True if the interpolator is local, False if global
      """
      self._is_local = is_local

   @method_config.setter
   def method_config(self, method_config):
      """!  
         Set the method-specific configuration object
         @param method_config method-specific configuration object
      """
      self._method_config = method_config

   @interp_config.setter
   def interp_config(self, interp_config):
      """!  
         Set the extra configuration object (combining model+method configuration)
         @param interp_config extra configuration object
      """
      self._interp_config = interp_config


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

#   # -----------------------------------------------------------------------------------------------
#   def learn(self, job, worker, model, spatial_data=None):
#      """!
#         Learn about the data prior to predition
#         @param job the SpredictJob object corresponding to the image to analyze
#         @param worker worker process object
#         @return ???
#      """
#      pass
      

   # -----------------------------------------------------------------------------------------------
   def adjust_method_config(self, config_dico, model, param_name, job, worker):
      """!
         Overrides method-specific configuration with identical definitions (same keys) found
         in file "interpolate.cfg".
         @param config_dico dictionary of keys defined in the [xxx.OPTIONS] section where xxx is the
                            method name
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance
         @note In order to be overridden, the key of a (key, value) pair must be the same as in 
               the method-specific configuration. 
               For instance if the DEGREE key is defined in "polyfit.cfg", the DEGREE key will be
               overridden if found in "interpolate.cfg".
      """
               
      section_name = "{0}.{1}".format(model.name.upper(), self.name)

      for key in config_dico.keys():
         try:
            key_dico = self.interp_config.get_as_dict(key, section_name)
            if param_name in key_dico:
               config_dico[key] = eval("key_dico[param_name]")   
               #print "overriding:", param_name, key, config_dico[key], type(config_dico[key])

         except SConfig.KeyNotFound:
            # --- Config key not found, skip it
            pass

      return config_dico 

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _load_method_config(self, config_dir):
      """! Load the method-specific configuration file and return the configuration object """

      # --- Config name is derived from name FittingElement          
      config_filepath = os.path.join(config_dir, self.name) + '.cfg'   
      try:
         return SConfig(config_filepath)
      except:
         # No configuration file found => The user must manually provide the relevant information
         self.helper.print_error("error reading configuration file {0}.".format(config_filepath)) 
         return None

   # -----------------------------------------------------------------------------------------------
   def _load_interp_config(self, config_dir):
      """! Load the combined model-method configuration file and return the configuration object """

      # --- Config name is derived from name FittingElement          
      config_filepath = os.path.join(config_dir, 'interpolate.cfg') 
      try:
         return SConfig(config_filepath)
      except:
         # No configuration file found => The user must manually provide the relevant information
         self.helper.print_error("error reading configuration file {0}.".format(config_filepath)) 
         return None


# -------------------------------------------------------------------------------------------------
class SpatialInterpolator(object):
   
   """! 
      Spatial interpolation of one or more model parameters over a two-dimensional grid
   """

   def __init__(self, master):
      """! SpatialInterpolator constructor """

      self._helper  = SpredictHelper()             # helper utility functions
      self._plotter = SpredictPlotter()            # plotter
      self._cleaner = DataCleaner(master)          # outlier detection and removal 
      self._validator = Validator(master)          # validation by resampling

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the SpredictHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the SpredictPlotter instance. """
      return self._plotter

   @property
   def cleaner(self):
      """! @return the DataCleaner instance. """
      return self._cleaner

   @property
   def validator(self):
      """! @return the Validator instance. """
      return self._validator

   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def interpolate(self, request, job, worker):
      """!
         Spatially interpolate
         @param request a InterpolationRequest object
         @param job the SpredictJob object corresponding to the image to analyze
         @param worker worker process object
         @return a dictionary with interpolation results and extra information
      """
   
      # --- Check validity of request
      self._check_interpolation_request(request, job, worker)

      # --- Process the interpolation request. The result dictionary contains 
      #     interpolated model parameters values and extra information about the interpolation
      result_dico = self._process_interpolation_request(request, job, worker)

      return result_dico


   # -----------------------------------------------------------------------------------------------
   class InterpolationError(Exception):
      """! 
         Exception thrown when an interpolation error has occurred
      """
      def __init__(self, msg):
         """!
            @param msg error message 
         """
         self._msg = msg
         self._exc_info = sys.exc_info()[1]

      def __str__(self):
         """! String representation of the InterpolationError class """
         return "Spredict *** ERROR ***: {0}".format(self._msg)


   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _check_interpolation_request(self, request, job, worker):
      """! Check validity of the interpolation request, paths to catalogs in particular """

     
      # --- Check path of input catalog file paths  
      if request.input_catalog_filepath is None:
         msg = "{0} - Could not find input catalog for: /{1}/image-{2:03d}-{3:1d}"\
             " => Job will not be processed".format(
             worker.name, job.get_branch_tree(), job.img_no, job.epoch)  
         raise SpatialInterpolator.InterpolationError(msg)         

      # --- Check path of target catalog file paths  
      if request.target_catalog_filepath is None:
         msg = "{0} - Could not find target catalog for: /{1}/image-{2:03d}-{3:1d}"\
             " => Job will not be processed".format(
             worker.name, job.get_branch_tree(), job.img_no, job.epoch)  
         raise SpatialInterpolator.InterpolationError(msg)  

      # --- Check path of model catalog file paths  
      if request.model_catalog_filepath is None:
         msg = "{0} - Could not find model catalog for: /{1}/image-{2:03d}-{3:1d}"\
             " => Job will not be processed".format(
             worker.name, job.get_branch_tree(), job.img_no, job.epoch)  
         raise SpatialInterpolator.InterpolationError(msg)  


   # -----------------------------------------------------------------------------------------------
   def _process_interpolation_request(self, request, job, worker):
      """! Process an interpolation request """
      
      # --- For each tile, supply *true* field input coordinates (not tile coordinates) and Z 
      #     values to the interpolation method and obtain interpolated values for the model 
      #     parameter(s) bound to the job

      # --- Dictionary where the results will be contained
      result_dico = {}
      result_dico["interpolated"] = {}       # results from interpolation
      result_dico["cross-validated"] = {}    # results from cross-validation

      # --- Max number of tiles to process for debugging
      max_nb_tiles = worker.config.get_as_int("MAX_NB_TILES", "DEBUGGING")

      tile_count = 0
      exit_loop = False   

      # --- for each parameter managed by the job...
      for param_name in job.param_names:

         # --- for each tile found in the input catalog...
         for x_tile_no in xrange(request.input_nb_tiles):
            for y_tile_no in xrange(request.input_nb_tiles):
            
               if not (x_tile_no, y_tile_no) in result_dico["interpolated"]:
                  result_dico["interpolated"][(x_tile_no, y_tile_no)] = {}
                  result_dico["cross-validated"][(x_tile_no, y_tile_no)] = {}

               # --- Invoke interpolation method, passing field input positions and related Z 
               #     values of the model parameter(s)      
               input_coords  = request.input_field_positions_per_tile_dico[((x_tile_no, y_tile_no))] 
               target_coords = request.target_field_positions_per_tile_dico[((x_tile_no, y_tile_no))] 
               input_values  = request.input_values_per_tile_dico[param_name][((x_tile_no, y_tile_no))] 

               # --- Check the nb. of input values is the same as the nb. of coordinates 
               if len(input_coords) == len(input_values):

                  # --- If requested, detect and remove outliers from input values
                  model_outlier_section_name = "{0}.{1}.{2}".format(
                                                        request.model.name.upper(), 
                                                        "OUTLIERS", param_name)
                  if request.model.config.get_as_boolean("REMOVE_OUTLIERS", 
                                                         model_outlier_section_name):

                     input_coords, input_values, nb_removed = self.cleaner.remove_outliers(
                                                                               request.model,
                                                                               param_name,
                                                                               input_coords, 
                                                                               input_values, 
                                                                               job, worker)   
                     if worker.logging_enabled():
                        worker.logger.log_info_p(
                          "{0} - {1} - Tile: ({2},{3}) - Removed {4} outliers for parameter: "\
                          "{5}  ({6})...".format(
                             worker.name, request.method.name, x_tile_no, y_tile_no, nb_removed,
                             param_name, request.model.name))

                  
                  if worker.logging_enabled():
                     worker.logger.log_info_p(
                       "{0} - {1} - Tile: ({2},{3}) - Interpolating parameter: "\
                       "{4}  ({5}) on {6} coordinates...".format(
                          worker.name, request.method.name, x_tile_no, y_tile_no,
                          param_name, request.model.name, len(input_coords)))

                  # --- Prediction of target values at target coordinates
                  result_dico["interpolated"][(x_tile_no, y_tile_no)][param_name] =\
                                                                            request.method.predict( 
                                                                            request.model,
                                                                            param_name,  
                                                                            input_coords, 
                                                                            target_coords, 
                                                                            input_values, 
                                                                            job, worker)


                  # --- Cross-Validate results if requested
                  CV_section_name = "{0}.CV.{1}".format(request.model.name.upper(), param_name)
                  if request.model.config.get_as_boolean("CROSS_VALIDATE", CV_section_name):      

                     if worker.logging_enabled():
                        worker.logger.log_info_p(
                          "{0} - {1} - Tile: ({2},{3}) - Cross-Validating parameter: "\
                          "{4}  ({5}) on {6} coordinates...".format(
                             worker.name, request.method.name, x_tile_no, y_tile_no,
                             param_name, request.model.name, len(input_coords)))

                     result_dico["cross-validated"][(x_tile_no, y_tile_no)][param_name] =\
                                                                    self.validator.cross_validate(
                                                                          request.method,
                                                                          request.model,
                                                                          param_name,  
                                                                          x_tile_no, 
                                                                          y_tile_no,
                                                                          input_coords, 
                                                                          input_values, 
                                                                          job, worker)
             
               else:
                  if worker.logging_enabled():
                     worker.logger.log_info_p(
                       "{0} - /{1}/image-{2:03d}-{3:1d} - The nb. of input coordinates {4} "\
                       "should be the same as the nb. of input values {5}".format(
                          worker.name, job.get_branch_tree(), job.img_no, job.epoch, 
                          len(input_coords), len(input_values)))

               # --- For debugging   
               tile_count += 1
               if max_nb_tiles != -1 and tile_count >= max_nb_tiles: 
                  exit_loop = True
                  break 

            # --- end for

            if exit_loop:
               break

         # --- end for    
     
      #print "result_dico:", result_dico.keys()   

      return result_dico


# -------------------------------------------------------------------------------------------------
class InterpolationRequest(object):
   """! 
      Represent a request for the spatial interpolation of the parameters of a model over a 
      two-dimensionbal grid
   """

   def __init__(self, catalog_file_type, model, method, sp_helper, job, worker):
      """! 
         Construct an InterpolationRequest object 
         @param catalog_file_type an image file type (like starfield_image.fits or 
                                deep-starfield_image.fits)
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param method object of a type inheriting from the BaseMethod class (RBF)
         @param job a SpredictJob object
         @param worker Worker process object
      """         

      # --- Interpolation model and method
      self._model = model
      self._method = method

      # --- Preferred catalog file extension
      input_catalog_file_extension = worker.config.get_as_string("CATALOG_EXTENSION", 
                                                                 "CATALOG_PROPERTIES.INPUT")      
      target_catalog_file_extension = worker.config.get_as_string("CATALOG_EXTENSION", 
                                                                  "CATALOG_PROPERTIES.TARGET")      
      model_catalog_file_extension = worker.config.get_as_string("CATALOG_EXTENSION", 
                                                                 "CATALOG_PROPERTIES.MODEL")      

      # --- Format for .txt catalogs
      input_catalog_is_SE = worker.config.get_as_boolean("CATALOG_IS_SEXTRACTOR",
                                                         "CATALOG_PROPERTIES.INPUT")
      target_catalog_is_SE = worker.config.get_as_boolean("CATALOG_IS_SEXTRACTOR",
                                                          "CATALOG_PROPERTIES.TARGET") 
      model_catalog_is_SE = worker.config.get_as_boolean("CATALOG_IS_SEXTRACTOR",
                                                         "CATALOG_PROPERTIES.MODEL") 

      # --- Absolute file paths for the required images and catalogs
      self._input_catalog_filepath = sp_helper.get_catalog_filepath(
                                                          job.input_catalog_filepath_dico,
                                                          input_catalog_file_extension, job)
      self._target_catalog_filepath = sp_helper.get_catalog_filepath(
                                                          job.target_catalog_filepath_dico,
                                                          target_catalog_file_extension, job)
      self._model_catalog_filepath = sp_helper.get_catalog_filepath(
                                                          job.model_catalog_filepath_dico,
                                                          model_catalog_file_extension, job)

      # --- Number of tiles in input and target catalogs
      input_hdu_no  = worker.config.get_as_int("CATALOG_HDU_NO", "CATALOG_PROPERTIES.INPUT")
      target_hdu_no = worker.config.get_as_int("CATALOG_HDU_NO", "CATALOG_PROPERTIES.TARGET")
      model_hdu_no  = worker.config.get_as_int("CATALOG_HDU_NO", "CATALOG_PROPERTIES.MODEL")
      self._input_nb_tiles = sp_helper.get_nb_tiles(
                                                self._input_catalog_filepath, 
                                                input_catalog_file_extension, 
                                                input_catalog_is_SE,
                                                input_hdu_no, 
                                                job, worker)
      self._target_nb_tiles = sp_helper.get_nb_tiles(
                                                 self._target_catalog_filepath, 
                                                 target_catalog_file_extension, 
                                                 target_catalog_is_SE,
                                                 input_hdu_no, 
                                                 job, worker)

      # --- Labels of columns that represents X and Y input and target positions

      input_XY_labels  = worker.config.get_as_list("INPUT_XY_COLUMN_LABELS", "INTERPOLATION")
      target_XY_labels = worker.config.get_as_list("TARGET_XY_COLUMN_LABELS", "INTERPOLATION")

      input_field_XY_labels  = worker.config.get_as_list("INPUT_FIELD_XY_COLUMN_LABELS", 
                                                         "INTERPOLATION")
      target_field_XY_labels = worker.config.get_as_list("TARGET_FIELD_XY_COLUMN_LABELS", 
                                                         "INTERPOLATION")

#      # --- Load input and target tile index per coordinate
#      self._input_tile_index_dico = self._get_tile_index_dico(
#                                                              self._input_catalog_filepath, 
#                                                              input_catalog_file_extension,
#                                                              input_catalog_is_SE,
#                                                              input_hdu_no, job, worker)

#      self._target_tile_index_dico = self._get_tile_index_dico(
#                                                              self._target_catalog_filepath, 
#                                                              target_catalog_file_extension,
#                                                              target_catalog_is_SE,
#                                                              target_hdu_no, job, worker)

#      # --- Load input and target tile position per coordinate
#      self._input_tile_pos_dico = self._get_tile_pos_dico(self._input_catalog_filepath, 
#                                                          input_catalog_file_extension,
#                                                          input_catalog_is_SE,
#                                                          input_hdu_no, job, worker)
#      self._target_tile_pos_dico = self._get_tile_pos_dico(self._target_catalog_filepath, 
#                                                          target_catalog_file_extension,
#                                                          target_catalog_is_SE,
#                                                          target_hdu_no, job, worker)

#      # --- Load field position per coordinate
#      self._input_field_pos_dico = self._get_field_pos_dico(self._input_catalog_filepath, 
#                                                            input_catalog_file_extension,
#                                                            input_catalog_is_SE,
#                                                            input_hdu_no, job, worker)
#      self._target_field_pos_dico = self._get_field_pos_dico(self._target_catalog_filepath, 
#                                                             target_catalog_file_extension,
#                                                             target_catalog_is_SE,
#                                                             target_hdu_no, job, worker)

      # --- Load input and target positions per tile
      self._input_positions_per_tile_dico = sp_helper.get_positions_per_tile_dico(
                                                            self._input_catalog_filepath, 
                                                            input_catalog_file_extension,
                                                            input_catalog_is_SE, input_hdu_no, 
                                                            input_XY_labels, self._input_nb_tiles,
                                                            job, worker)
      self._target_positions_per_tile_dico = sp_helper.get_positions_per_tile_dico(
                                                            self._target_catalog_filepath, 
                                                            target_catalog_file_extension,
                                                            target_catalog_is_SE, target_hdu_no, 
                                                            target_XY_labels, self._target_nb_tiles,
                                                            job, worker)

      # --- load input model parameter values per tile (for param(s) handled by this job) 
      self._input_values_per_tile_dico = sp_helper.get_values_per_tile_dico(
                                                            job.param_names,
                                                            self._model_catalog_filepath, 
                                                            model_catalog_file_extension,
                                                            model_catalog_is_SE,
                                                            model_hdu_no, self._input_nb_tiles,
                                                            job, worker)

      # --- Load true field positions per tile
      self._input_field_positions_per_tile_dico = sp_helper.get_field_positions_per_tile_dico(
                                                            self._input_catalog_filepath, 
                                                            input_catalog_file_extension,
                                                            input_catalog_is_SE, input_hdu_no, 
                                                            input_field_XY_labels, 
                                                            self._input_nb_tiles,
                                                            job, worker)

      self._target_field_positions_per_tile_dico = sp_helper.get_field_positions_per_tile_dico(
                                                            self._target_catalog_filepath, 
                                                            target_catalog_file_extension,
                                                            target_catalog_is_SE, target_hdu_no, 
                                                            target_field_XY_labels, 
                                                            self._target_nb_tiles,
                                                            job, worker)


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def model(self):
      """!
         Get the model whose parameters must be interpolated
         @return model whose parameters must be interpolated 
      """
      return self._model

   @property
   def method(self):
      """!
         Get the interpolation mthod object
         @return interpolation method of a type derived from BaseMethod 
      """
      return self._method

   @property
   def input_catalog_filepath(self):
      """!
         Get the input catalog absolute file paths
         @return input catalog absolute file paths 
      """
      return self._input_catalog_filepath

   @property
   def target_catalog_filepath(self):
      """!
         Get the target catalog absolute file paths
         @return target catalog absolute file paths 
      """
      return self._target_catalog_filepath

   @property
   def model_catalog_filepath(self):
      """!
         Get the model catalog absolute file paths
         @return model catalog absolute file paths 
      """
      return self._model_catalog_filepath

   @property
   def input_nb_tiles(self):
      """!
         Get the number of tiles in the input catalog
         @return number of tiles in the input catalog 
      """
      return self._input_nb_tiles

   @property
   def target_nb_tiles(self):
      """!
         Get the number of tiles in the target catalog
         @return number of tiles in the target catalog 
      """
      return self._target_nb_tiles

   @property
   def input_tile_index_dico(self):
      """! 
         Get the dictionary containing input tile indice per coordinates
         @return dictionary containing input tile indice per coordinates
      """
      return self._input_tile_index_dico

   @property
   def target_tile_index_dico(self):
      """! 
         Get the dictionary containing target tile indice per coordinates
         @return dictionary containing target tile indice per coordinates
      """
      return self._target_tile_index_dico

   @property
   def input_tile_pos_dico(self):
      """! 
         Get the dictionary containing input tile positions per coordinates
         @return dictionary containing input tile positions per coordinates
      """
      return self._input_tile_pos_dico

#   @property
#   def target_tile_pos_dico(self):
#      """! 
#         Get the dictionary containing target tile positions per coordinates
#         @return dictionary containing target tile positions per coordinates
#      """
#      return self._target_tile_pos_dico

#   @property
#   def input_field_pos_dico(self):
#      """! 
#         Get the dictionary containing input field positions per coordinates
#         @return dictionary containing input field positions per coordinates
#      """
#      return self._input_field_pos_dico

#   @property
#   def target_field_pos_dico(self):
#      """! 
#         Get the dictionary containing target field positions per coordinates
#         @return dictionary containing target field positions per coordinates
#      """
#      return self._target_field_pos_dico

   @property
   def input_tile_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing input tile positions per (x,y) tile pairs
         @return dictionary containing input tile positions per (x,y) tile pairs
      """
      return self._input_positions_per_tile_dico

   @property
   def target_tile_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing target tile positions per (x,y) tile pairs
         @return dictionary containing target tile positions per (x,y) tile pairs
      """
      return self._target_positions_per_tile_dico

   @property
   def input_field_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing input field positions per tile
         @return dictionary containing input field positions per tile
      """
      return self._input_field_positions_per_tile_dico

   @property
   def target_field_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing target field positions per tile
         @return dictionary containing target field positions per tile
      """
      return self._target_field_positions_per_tile_dico


   @property
   def input_stamp_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing input stamp positions per (x,y) tile pairs in an image
         @return dictionary containing input stamp positions per (x,y) tile pairs in an image
      """
      return self._input_stamp_positions_per_tile_dico

   @property
   def target_stamp_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing target stamp positions per (x,y) tile pairs in an image
         @return dictionary containing target stamp positions per (x,y) tile pairs in an image
      """
      return self._target_stamp_positions_per_tile_dico

   @property
   def input_field_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing input field positions per (x,y) tile pairs
         @return dictionary containing input field positions per (x,y) tile pairs
      """
      return self._input_field_positions_per_tile_dico

   @property
   def target_field_positions_per_tile_dico(self):
      """! 
         Get the dictionary containing target field positions per (x,y) tile pairs
         @return dictionary containing target field positions per (x,y) tile pairs
      """
      return self._target_field_positions_per_tile_dico

   @property
   def input_values_per_tile_dico(self):
      """! 
         Get the dictionary containing input values per (x,y) tile pairs
         @return dictionary containing input values per (x,y) tile pairs
      """
      return self._input_values_per_tile_dico

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~


# -- EOF sp_interp.py
