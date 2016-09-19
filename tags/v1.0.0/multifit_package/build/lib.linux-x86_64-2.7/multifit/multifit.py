"""! 
   @package multifit.multifit Model and Method-independent Fitting
   @author Marc Gentile
   @file multifit.py
   Model and Method-independent Fitting
""" 

# -- Python imports
import os, sys
import math
import string
import time
import imp
import numpy

# --- External imports
from sconfig import SConfig

from profilehooks import *
# -------------------------------------------------------------------------------------------------
class ModelParam(object):
   """! Represents a parameter in the model """

   def __init__(self, name, guess, bounds=[None, None]):
      """!
         Construct a ModelParam object
         @param name name of the model
         @param guess guess value
         @param bounds parameter min and max bounds in the form of list [min, max]   
         @note the minimum, maximum values or both may be set to @c None 
      """

      self._name = name             # name of parameter
      self._bounds = bounds         # allowed range (<min value>, <max value>)
      self._guess_value = guess     # guess value
      self._fitted_value = guess    # fitted value (initially equal to guess)

   def __str__(self):
      """! 
         String representation of a ModelParam object 
         @return a Python String describing the ModelParam object
      """
      return "<{0}, {1}, {2}, {3}>".format(self._name, 
                                        self._guess_value, self._bounds, self._fitted_value)

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

   @property
   def guess_value(self):
      """! @return the guess value of the parameter """   
      return self._guess_value

   @property
   def fitted_value(self):
      """! 
         @return the fitted value of the parameter once fitting is complete 
         @note the initial value of fitted_value() is the same as that returned by guess_value()
      """   
      return self._fitted_value

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

   @guess_value.setter 
   def guess_value(self, guess_value):
      """!     
         Set the guess value of the model parameter 
         @param guess_value guess value
      """
      self._guess_value = guess_value     

   @fitted_value.setter 
   def fitted_value(self, fitted_value):
      """! 
         Set the fitted value of the model parameter
         @param fitted_value fitted value
      """  
      self._fitted_value = fitted_value     


# --------------------------------------------------------------------------------------------------
class FittingElement(object):
   """! 
      A fitting model or algorithm/method. 
      @note this is the base class for a fitting method or model
   """

   def __init__(self, name, config_dir, config_filename):
      """!
         Construct a FittingElement object
         @name name of the element
         @config_dir directory of the configuration file for the element 
         @config_filename name of the configuration file for the element 
      """   
      self._name = name  # name of the model
      if name is not None:
         self._config = self._load_config(config_dir, config_filename)            # SConfig object
         self._default_config = self._load_config(config_dir, self.name + ".cfg") # SConfig object
      else:
         self._config = None
         self._default_config = None

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters   

   @property
   def name(self):
      """! @return the name of the element """   
      return self._name
      
   @property
   def config(self):
      """! @return the directory of the default configuration file for the element """   
      return self._config

   @property
   def default_config(self):
      """! @return the directory of the default configuration file for the element """   
      return self._default_config

   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _load_config(self, config_dir, config_filename):

      config_file = None

      # --- Config name is derived from name FittingElement          
      config_filepath = os.path.join(config_dir, config_filename)   
      default_config_filepath = os.path.join(config_dir, self.name) + '.cfg'   
      try:
         config_file = SConfig(config_filepath)
         if config_filepath != default_config_filepath:
            print("MultiFit *** Note *** using configuration file {0}".format(config_filepath))

      except:
         # No configuration file found => Try the default configuration name
         try:
            print("MultiFit *** Note *** using default configuration file {0}".format(
                                                                          default_config_filepath))
            config_file = SConfig(default_config_filepath)         
         except:
            # No configuration file found => The user must manually provide the relevant information
            print("MultiFit *** Warning ***: error reading configuration file {0}.".format(
                                                                         default_config_filepath))
      return config_file

# -------------------------------------------------------------------------------------------------
class BaseFittingModel(FittingElement):
   """! 
      Represents a parametric model to which the data must be fitted
      @note this is the base class for all fitting models
   """

   def __init__(self, name, config_dir, config_filename):
      """!
         Construct BaseFittingModel object
         @param name name of the model   
         @param config_dir configuration directory of the model
         @param config_filename configuration filename of the model
      """

      FittingElement.__init__(self, name, config_dir, config_filename) # initialize parent class

#      self._model_options = model_options    # options influencing how the model should be fitted
#      self._model_options = diag_options     # diagnose directives (logging, stats, etc.)

      # --- Private attributes (should not be accessed directly)
      self.__params = []         # the model parameters
      self.__param_dico = {}     # to access parameter by names instead of indice

      self._numerical_error_handler = None   # numerical error handler method
      self._out_of_bound_handler = None      # out-of-bound handler function

      # --- Model Parameters list: set from configuration or otherwise, manually
      model_params = self.get_params_from_config()

      if len(model_params) > 0:
         self.set_params(model_params)        # set parameters based on configuration          

      # --- Model function for evaluating the model given parameters
      self._func = self.get_model_func_from_config() # function for evaluating the model

      # --- Read catalog column mapping iformation if available
      self._col_mapping_dico = {}
      if self.config.has_section("MODEL.GUESS"):

         # Read catalog column mapping information
         if self.config.has_key("guess_col_mapping", "MODEL.GUESS"):
            self._col_mapping_dico = self.config.get_as_dict("guess_col_mapping", "MODEL.GUESS") 

      # --- Model check section if available
      if self.config.has_section("MODEL.CHECK"):

         # Read catalog column mapping information
         self._check_centroids = False
         self._max_centroid_shift = 2.0
         if self.config.has_key("check_centroids", "MODEL.CHECK"):
            self._check_centroids = self.config.get_as_boolean("check_centroids", "MODEL.CHECK")

            if self._check_centroids:
               if self.config.has_key("max_centroid_shift", "MODEL.CHECK"):
                  self._max_centroid_shift = self.config.get_as_int("max_centroid_shift", 
                                                                     "MODEL.CHECK")

#      #print "BaseFittingModel:", [(p.name, p.bounds, p.guess_value) for p in model_params]


   def __str__(self):
      """! 
         String representation of a BaseFittingModel object 
         @return a Python String describing the BaseFittingModel object
      """
      descr = "Fitting Model:\n"
      descr += "- name: {0}\n".format(self.name)
      if self.func is not None:
         descr += "- func: {0}\n".format(self.func.__name__)
      if len(self.params) > 0:
         descr += "- params: {0}\n".format([str(p) for p in self.params])

      return descr

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters  

   @property
   def params(self):
      """! @return the model parameters in the form of a list of ModelParam objects """
      return self.get_params()

   @property
   def func(self):
      """! @return the model function that builds a model given ModelParam objects """
      return self.get_func()

   @property
   def col_mapping_dico(self):
      """! 
         Return the dictionary for mapping catalog columns containing guess data (like SExtractor) 
         @return the dictionary for mapping catalog columns containing guess data
      """
      return self._col_mapping_dico

   @property
   def check_centroids(self):
      """! 
         Return whether fitted centroids should be checked or not
         @return True if fitted centroids must be checked 
      """
      return self._check_centroids

   @property
   def max_centroid_shift(self):
      """! 
         Return maximum allowed centroid shift if checking centroids
         @return the maximum allowed centroid shift if checking centroids
      """
      return self._max_centroid_shift

   @property
   def out_of_bound_handler(self):
      """! 
         Return the callback handler function to notify a parameter value came outside the 
                 allowed range of varation specified by Param.guess_value
         @return the callback handler function
      """
      #print "***MFIT: multifit, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """

      #print "***MFIT: multifit, getting  numerical_error_handler..."

      return self._numerical_error_handler


   # --- Setters  

   @params.setter
   def params(self, params):
      """!  
         Set the model parameters in the form of a list of ModelParam objects 
         @param params a list of ModelParam object
      """ 
      self.set_params(params)

   @func.setter
   def func(self, func):
      """!  
         Set the model function that builds a model given ModelParam objects 
         @param func the model function that builds a model given ModelParam objects
      """
      self._func = func


   @out_of_bound_handler.setter
   def out_of_bound_handler(self, out_of_bound_handler):
      """! 
         Set the callback handler function to notify a parameter value came outside the 
             allowed range of varation specified by Param.guess_value       
         @param out_of_bound_handler the out-of-bound callback function
      """
      self._out_of_bound_handler = out_of_bound_handler


   @numerical_error_handler.setter
   def numerical_error_handler(self, numerical_error_handler):
      """! 
         Set the callback function to notify a numerical error that occurred while computing 
             the parameter value       
         @param numerical_error_handler the numerical error callback function
      """
      self._numerical_error_handler = numerical_error_handler


   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def get_func(self):
      """! @return the model function that builds a model given ModelParam objects """
      return self._func

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
   def reset_params_from_config(self):
      """! Reset parameters from configuration data """

      model_params = self.get_params_from_config()
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
         @param name the name of the parameter
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
   def get_param_values(self):
      """! @return the parameter initial guess values """     
      return [p.guess_value for p in self.__params]

   # -----------------------------------------------------------------------------------------------
   def get_bounds(self):
      """! @return the parameter bounds in the form of a list [min, max] """     
      return [p.bounds for p in self.__params]

   # -----------------------------------------------------------------------------------------------
   def compute_model_func(self, func, params, kwargs):
      """! 
         Invoke a model function @c func with parameter list @c *params and arguments keywords
         @c **kwargs

         @param func the model function that builds a model given the ModelParam objects
         @param params a list of ModelParam objects
         @param kwargs a dictionary of parameter (key, value) pairs  
      """
      return func(*params, **kwargs) 

   # -----------------------------------------------------------------------------------------------
   def get_model_func_from_config(self):
      """! @return the model function from configuration """

      raise(ModelFitter.FittingError(
                "The method get_model_func_from_config() must be implemented in a derived class"))
            
   # -----------------------------------------------------------------------------------------------
   def get_params_from_config(self):
      """! @return the list of model parameters from the model configuration, assuming a key with
          name @c "params" has been defined in a section called @c "MODEL".
      """
      params = []

      # --- Read default config to make sure all parameters get defined
      if not self.default_config is None: 
         if self.default_config.has_section("MODEL.PARAMS"):

            for subsection in self.default_config.get_subsections("MODEL.PARAMS"): 

               # --- Check if subsection defined in the custom configuration, if not, take
               #     the valeus in the default one  
               if subsection in self.config.get_subsections("MODEL.PARAMS"): 

                  # --- Param name if not defined yet
                  name = self.config.get_as_string("name", subsection)

                  # --- Param value range
                  bounds = self.config.get_as_list("range", subsection)
                  if bounds is not None:

                     # --- Param guess value
                     if bounds[0] != bounds[1]: 
                        guess_str = self.config.get_as_string("guess", subsection)
                        if eval(guess_str) is not None:
                           guess = float(guess_str)
                        else:
                           guess = None    
                     else:
                        if not bounds[0] is None:
                           # --- Identical bounds: fix guess value to that of the boundary
                           guess = bounds[0]     
                           if bounds[0] is not None:
                              print("MultiFit *** Note *** parameter {0} " \
                                    "will remain constant during fitting with value {1:.2f}".format(
                                                                                  name, bounds[0]))
                        else:
                           guess_str = self.config.get_as_string("guess", subsection)
                           if eval(guess_str) is not None:
                              guess = float(guess_str)
                           else:
                              guess = None    

               else:

                  # --- Take the default definition from the default configuration

                  # --- Param name if not defined yet
                  name = self.default_config.get_as_string("name", subsection)

                  # --- Param value range
                  bounds = self.default_config.get_as_list("range", subsection)
                  if bounds is not None:

                     # --- Param guess value
                     if bounds[0] != bounds[1]: 
                        guess_str = self.default_config.get_as_string("guess", subsection)
                        if not eval(guess_str) is None:
                           guess = float(guess_str)
                        else:
                           guess = None    
                     else:
                        # --- Identical bounds: fix guess value to that of the boundary
                        guess = bounds[0]      
                        if guess is not None:
                           print("MultiFit *** Note *** parameter {0} " \
                                 "will remain constant during fitting with value {1:.2f}".format(
                                                                                  name, bounds[0]))
               params.append(ModelParam(name, guess, bounds)) 

      return params


   # -----------------------------------------------------------------------------------------------
   def check_fitted_params(self, (x, y), actual_galaxy_image, 
                                         galaxy_info_dico, se_galaxy_info_dico, 
                                         guess_values, param_values):


      """!
         Check that the parameter values  @c param_values are within the allowed bounds
         @param param_values values of the parameters to check against the bounds
      """
      
      # --- Check parameter ranges
      lower_valid = [ v >= p.bounds[0] for (p, v) in zip(self.params, param_values) \
                                       if p.bounds[0] is not None ]
      upper_valid = [ v <= p.bounds[1] for (p, v) in zip(self.params, param_values) \
                                       if p.bounds[1] is not None ]
      valid = numpy.logical_and(numpy.all(lower_valid), numpy.all(upper_valid))
      if valid:
         if self.check_centroids:
            # --- Check parameter centroid shifts
            x0 = (actual_galaxy_image.shape[1]-1.0)/2.0
            y0 = (actual_galaxy_image.shape[0]-1.0)/2.0
            xc = param_values[self.get_param_names().index("xc")]
            yc = param_values[self.get_param_names().index("yc")]
            valid = (math.fabs(xc - x0) <= self.max_centroid_shift) and \
                    (math.fabs(yc - y0) <= self.max_centroid_shift) 

      return numpy.all(numpy.asarray(valid))


   # -----------------------------------------------------------------------------------------------
   def get_guess_value(self, coords, galaxy_info_dico, se_galaxy_info_dico, 
                              param_name, guess_col_name, default_value):
      """! Pick a guess value for a parameter using the data of a galaxy catalog """

      col_name = self.col_mapping_dico.get(guess_col_name, None)
      if not col_name is None:

         # --- Determine source dictionary
         if col_name.startswith("SE_"):
            data_dico = se_galaxy_info_dico
         else:
            data_dico = galaxy_info_dico

         try:
            guess_value = data_dico[col_name][coords]
         except Exception as detail:
            # --- Could not get guess value => will abort the fitting process
            print("*** Multifit: {0} model *** Could not retrieve guess value for {1} "\
                  "from guess_column {2} at coordinates {3} - {4}".format(
                     self.name, param_name, col_name, coords, detail))
            print("*** Multifit: {0} model *** Using '{2}' as default calue...".format(
                                                      self.name, param_name, default_value))
            guess_value = default_value
      else:
         print("*** Multifit: {0} model *** No mappped guess column found for {1} - "\
               "Using '{2}' as default value...".format(
                                                   self.name, param_name, default_value))
         guess_value = default_value

      return guess_value 


   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _import_module(self, module_name, module_dir):
      """! Find and load a module <module_name> from a directory <module_dir> """

      try:

         file_obj, filename, data = imp.find_module(module_name, [module_dir])
         imp.load_module(module_name, file_obj, filename, data)
         return imp.load_module(module_name, file_obj, filename, data)

      except:
         raise(ModelFitter.FittingError("Could not load module {0}.".format(module_name)))



# -------------------------------------------------------------------------------------------------
class BaseFittingMethod(FittingElement):

   """! 
      Represent a fitting method/algorithm (like Levenberg-Marquardt) 
      @note this is the base class for all fitting methods
   """

   def __init__(self, name, 
                      method_config_dir=None, 
                      fitting_config_dir=None, 
                      method_config_filename=None):
      """!  
         Construct a BaseFittingMethod object
         @param name name of the method   
         @param method_config_dir configuration directory of the method
         @param fitting_config_dir directory of extra configuration files for the method
         @param method_config_filename name of method configuration file
      """

      FittingElement.__init__(self, name, method_config_dir, method_config_filename)

      self._objective_func  = None     # must be defined in the derived method class
      #self._model_fitting_func = None  # f(model) component of the objective function
      self._fitting_options = None     # options influencing fitting (maxfev, etc.)
      self._diag_options = None        # diagnose directives (logging, stats, etc.)

      self._fitting_status = BaseFittingMethod.FittingStatus(self) # method-indep. fitting status
      self._fitting_info  = None       # method-specific fitting information (meta-data like maxfev)       
      self._fitting_data  = None       # method-specific collected fitting data if any collected      
      self._fitting_stats = None       # method-specific collected statistics if any computed      
      self._fitted_params = None       # hold fitted value once fitting is complete

      self._fitting_config = None      # fitting config object if any

      self._fitting_config_dir = fitting_config_dir   # fitting directives

   def __str__(self):
      """! 
         String representation of a BaseFittingMethod object 
         @return a Python String describing the BaseFittingMethod object
      """
      descr = "Fitting Method {0}:\n".format(self.name)
      if self.fitting_options is not None:
         descr += str(self.fitting_options)
      if self.diag_options is not None:
         descr += str(self.diag_options)
      descr += str(self.fitting_status)
      if self.fitting_info is not None:
         descr += str(self.fitting_info)
      if self.fitted_params is not None:
        descr += str(self.fitted_params)

      return descr

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~
   
   # --- Getters 

   @property
   def fitting_config(self):
      """! @return configuration object for fitting """
      return self._fitting_config

   @property
   def fitting_config_dir(self):
      """! @return configuration directory for fitting """
      return self._fitting_config_dir

   @property
   def objective_func(self):
      """! @return the objective function used for fitting """
      return self._objective_func

#   @property
#   def model_fitting_func(self):
#      """! @return the model fiting function """
#      return self._model_fitting_func

   @property
   def fitted_params(self):
      """! @return the fitted values or all parameters in the form of a list of ModelParam objects,
           once fitting is complete """
      return self._fitted_params

   @property
   def fitting_options(self):
      """! 
         @return fiting options influencing the fitting process (max iterations, accuracy, etc.)
         @note the returned objects will be derived from BaseFittingOptions, the exact type
               depending on the fitting method chosen  
      """
      return self._fitting_options

   @property
   def diag_options(self):
      """! @return diagnostic options for monitoring the fitting process (logs, stats, etc.). """
      return self._diag_options

   @property
   def fitting_status(self):
      """! 
         @return the status of the fitting process (errno, messages, warnings) in the form of a 
                 BaseFittingMethod.FittingStatus object
      """
      return self._fitting_status

   @property
   def fitting_info(self):
      """! @return method-specific fitting information in the form of an object derived from 
           BaseFittingInfo object once fitting is completed
         @note the exact class of the object returned will depend on the chosen method 
      """
      return self._fitting_info

   @property
   def fitting_data(self):
      """! 
         @return method-specific collected fitting data if any, collected once fitting is 
                complete 
         @note return @c None if no fitting data available   
      """
      return self._fitting_data

   @property
   def fitting_stats(self):
      """! 
         @return method-specific collected statistics if any, computed once fitting is completed. 
         @note return @c None if no fitting statistics available   
      """
      return self._fitting_stats

   # --- Setters  

   @fitting_config.setter
   def fitting_config(self, fitting_config):
      """! 
         Set the configuration object for fitting
         @param fitting_config configuration object for fitting
      """
      self._fitting_config = fitting_config

   @objective_func.setter
   def objective_func(self, objective_func):
      """! 
         Set the objective function used for fitting
         @param objective_func objective function
      """
      self._objective_func = objective_func

#   @model_fitting_func.setter
#   def model_fitting_func(self, model_fitting_func):
#      """! 
#         Set the model_fitting function
#         @param model_fitting_func model_fitting function
#      """
#      self._model_fitting_func = model_fitting_func

   @fitting_status.setter
   def fitting_status(self, fitting_status):
      """! 
         Set the status of the fitting process (errno, messages, warnings) in the form of a 
         BaseFittingMethod.FittingStatus object 
         @param fitting_status a BaseFittingMethod.FittingStatus object
      """
      self._fitting_status = fitting_status

   @fitting_options.setter
   def fitting_options(self, fitting_options):
      """!  
         Set fitting options influencing the fitting process (max iterations, accuracy, etc.)
         @param fitting_options  
         @note the parameter must of a type derived from BaseFittingOption 
      """ 
      self._fitting_options = fitting_options

   @diag_options.setter
   def diag_options(self, diag_options):
      """! 
         Set diagnostic options for monitoring the fitting process (logs, stats, etc.) 
         @param diag_options  
         @note the parameter must of a type derived from BaseDiagOption 
      """
      self._diag_options = diag_options


   # ~~~~~~~~~~~~~
   # Inner Classes
   # ~~~~~~~~~~~~~
   
   # -----------------------------------------------------------------------------------------------
   class FittingStatus(object):
      """! 
         Hold the status of a fitting operation
         @note this class is independent of the fitting method selected
      """

      def __init__(self, method, residuals=None, errno=0, msg="", warn_msgs=[], custom_data=None):
         """!
            Construct a FittingStatus object for a given method
            @param method an object derived from BaseFittingMethod
            @param residuals fitting residuals
            @param errno an error code
            @param msg an error message 
            @param warn_msgs a list of warning messages
            @param custom_data custom data
            @note @c errno will be negative in case of error, zero otherwise
         """

         self._method = method   # owner method instance 
         self._errno = errno     # execution status: 0: success, < 0:failure
         self._msg = msg         # execution status message
         self._warn_msgs = warn_msgs      # warning messages, if any
         self.residuals = residuals       # flattened residuals
         self._custom_data = custom_data  # custom data, if any 

      def __str__(self):
         """! 
            String representation of a FittingStatus object 
            @return a Python String describing the FittingStatus object
         """
         descr = "Fitting Status for {0}:\n".format(self.method.name)
         descr += "- errno: {0}\n".format(self.errno)
         descr += "- message: {0}\n".format(self.msg)
         descr += "- warnings: {0}\n".format([w for w in self.warn_msgs])
         descr += "- residuals: {0}\n".format(self.residuals)
         if self.custom_data is not None:
            descr += "- custom data: {0}\n".format(self.custom_data)
         else:
            descr += "- custom data: None\n" 
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~

      # --- Getters ---

      @property
      def method(self):
         """! @return the associated method instance as an object derived from BaseFittingMethod """
         return self._method

      @property
      def errno(self):
         """! 
            Returns the fitting execution error code
            @retval " >= 0" success
            @retval "< 0" failure
         """
         return self._errno

      @property
      def msg(self):
         """! 
            @return the fitting execution status message 
         """
         return self._msg

      @property
      def warn_msgs(self):
         """! @return a list of warning messages (may be empty) """
         return self._warn_msgs

      @property
      def custom_data(self):
         """! @return custom data collected by model fiting if any, may be None """
         return self._custom_data


      # --- Setters ---

      @errno.setter
      def errno(self, errno):
         """! 
            Returns the fitting execution error code 
            param errno error code (0 for success, < 0 for failure)
         """   
         self._errno = errno

      @msg.setter
      def msg(self, msg):
         """! Set the fitting execution status message """
         self._msg = msg

      @warn_msgs.setter
      def warn_msgs(self, warn_msgs):
         """! Set the fitting execution status warning message """
         self._warn_msgs = warn_msgs

      @custom_data.setter
      def custom_data(self, custom_data):
         """! Set custom data collected by model fiting if any """
         self._custom_data = custom_data

      # ~~~~~~~~~~~~~~~
      # Public methods
      # ~~~~~~~~~~~~~~~

      # --------------------------------------------------------------------------------------------
      def add_warn_msg(self, warn_msg):
         """! 
            Add a warning message 
            @param warn_msg the warning message to add 
         """
         if not warn_msg in self._warn_msgs:
            self._warn_msgs.append(warn_msg)         

      # --------------------------------------------------------------------------------------------
      def set_fitting_status(self, errno, msg):
         """! 
            Set the fititng status of the fit
            @param errno an error code
            @param msg an error message 
         """
         self._errno = errno
         self._msg = msg

      # --------------------------------------------------------------------------------------------
      def clear(self):
         """! Clear the fitting status """
         self._errno = 0         # execution status:  >= 0: success, < 0:failure
         self._msg = ''          # execution status message
         del self._warn_msgs[:]  # list of warning messages (may be empty)



# --------------------------------------------------------------------------------------------------
class BaseFittingOptions(object):

   """! 
      Represents options that may influence the fitting process 
      @note this is the base class for all method-specific fitting options
   """

   def __init__(self, method):
      """! Construct a BaseFittingOptions object """
      self._method = method    # the fitting method

   def __str__(self):
      """!
         String representation of a BaseFittingOptions object 
         @return a Python String describing the BaseFittingOptions object
      """
      return "Fitting options for {0}:\n".format(self.method.name)

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~

   @property
   def method(self):
      """! @return the associated fitting method as an object derived from BaseFittingMethod """
      return self._method


# --------------------------------------------------------------------------------------------------
class BaseDiagOptions(object):

   """! 
      Represents options that can be used to collect diagonotic information about the fitting 
      process 
      @note this is the base class for all method-specific diagnostic options
   """
   def __init__(self, method):
      """! Construct a BaseDiagOptions object """
      self._method = method    # the fitting method

      # --- Private properties  
      self._diag_flags = 0
      self._base_output_dir = "."        
      self._log_output_dir   = os.path.join(self._base_output_dir, 'logs')
      self._plot_output_dir  = os.path.join(self._base_output_dir, 'plots')
      self._stats_output_dir = os.path.join(self._base_output_dir, 'stats')
      self._log_file_name = 'multifit_{0}.log'.format(time.strftime("%d.%m.%y_%H.%M.%S",
                                                      time.localtime()))

      # --- Initialize options from configuration
      self._init_from_config(method.config)

   def __str__(self):
      """!
         String representation of a BaseDiagOptionss object 
         @return a Python String describing the BaseDiagOptionsobject
      """

      descr =  "Diagnostic options for {0}:\n".format(self.method.name)  
      descr += "- diagnostic flags:{0}\n".format(self._diag_flags) 
      descr += "- log file name: {0}\n".format(self._log_file_name)
      descr += "- base output dir: {0}\n".format(self._base_output_dir)
      descr += "- log output dir: {0}\n".format(self._log_output_dir)
      descr += "- plot output dir: {0}\n".format(self._plot_output_dir)
      descr += "- stats output dir: {0}\n".format(self._stats_output_dir)

      return descr

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~   

   # --- Getters 

   @property
   def method(self):
      """! @returns the associated fitting method. """
      return self._method

   @property
   def diag_flags(self):
      """! 
         @return the diagnostic flags 
      """
      return self._diag_flags

   @property
   def base_output_dir(self):
      """! 
         @return the base output directory for logs, plots and stats 
      """
      return self._base_output_dir

   @property
   def log_output_dir(self):
      """! 
         @return the output directory for logging 
      """
      return self._log_output_dir

   @property
   def plot_output_dir(self):
      """! 
         @return the output directory for plotting 
      """
      return self._plot_output_dir

   @property
   def stats_output_dir(self):
      """! 
         @return the output directory for outputing statistics 
      """
      return self._stats_output_dir

   @property
   def log_file_name(self):
      """! 
         @return the name of the log file 
      """
      return self._log_file_name

   # --- Setters 

   @diag_flags.setter 
   def diag_flags(self, diag_flags):
      """! 
         Set the diagnostic flags
         @param diag_flags diagnostic flagsoutputing 
      """
      self._diag_flags = diag_flags

   @base_output_dir.setter 
   def base_output_dir(self, base_output_dir):
      """! 
         Set the base output directory for logs, plots and stats 
         @param base_output_dir base output directory
      """
      self._base_output_dir = base_output_dir

   @log_output_dir.setter 
   def log_output_dir(self, log_output_dir):
      """! 
         Set the output directory for logging 
         @param log_output_dir output directory for logging
      """
      self._log_output_dir = log_output_dir

   @plot_output_dir.setter 
   def plot_output_dir(self, plot_output_dir):
      """! 
         Set the output directory for plotting 
         @param plot_output_dir output directory for plotting
      """
      self._plot_output_dir = plot_output_dir

   @stats_output_dir.setter 
   def stats_output_dir(self, stats_output_dir):
      """! 
         Set the output directory for statistics 
         @param stats_output_dir output directory for outputing statistics
      """
      self._stats_output_dir = stats_output_dir

   @log_file_name.setter 
   def log_file_name(self, log_file_name):
      """! 
         Set the name of the log file  
         @param log_file_name name of the log file
      """
      self._log_file_name = log_file_name

   # ~~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~~

   def _init_from_config(self, config):

      # --- Default configuration
      if not config is None:
         method_section_name = self.method.name.upper()
         model_options_section_name = "{0}.OPTIONS".format(method_section_name)

         if config.has_section(method_section_name):
            if config.has_key("log_file_name", model_options_section_name ):
               self._log_file_name = config.get_as_string("log_file_name", 
                                                          model_options_section_name)            
            if config.has_key("base_output_dir", model_options_section_name):
               self._base_output_dir = config.get_as_string("base_output_dir", 
                                                            model_options_section_name)            
            if config.has_key("log_output_dir", model_options_section_name):
               self.log_output_dir = config.get_as_string("log_output_dir", 
                                                          model_options_section_name)            
            if config.has_key("plot_output_dir", model_options_section_name):
               self._plot_output_dir = config.get_as_string("plot_output_dir", 
                                                            model_options_section_name)            
            if config.has_key("stats_output_dir", model_options_section_name):
               self._stats_output_dir = config.get_as_string("stats_output_dir", 
                                                             model_options_section_name)         
      else:
         if not self.method.name is None:
            print("MultiFit *** Warning ***: configuration file notfound for {0}.".format(
                                                                                 self.method.name))   

# --------------------------------------------------------------------------------------------------
class BaseFittingInfo(object):

   """! 
       Fitting information available from a method once fitting is complete. 
   """

   def __init__(self, method):
      
      self._method = method    # the fitting method

   def __str__(self):

      descr = "Fitting Info for {0}:\n".format(self.method.name)
      descr += "- Fitted parameters: {0}\n".format(self.method.fitted_params)
      return descr

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property
   def method(self):
      """! @return the associated fitting method as an object derived from BaseFittingMethod """
      return self._method



# -------------------------------------------------------------------------------------------------
class ModelFitter(object):
   """!
      Fit a parametric model using a fitting method
   """

   def __init__(self, model_name=None, method_name=None, 
                      model_config_filename=None, 
                      method_config_filename=None,
                      model_config_dir="./config/models", 
                      method_config_dir="./config/methods", 
                      fitting_config_dir = "./config/fitting", 
                      model_module_dir="./modules/models", 
                      method_module_dir="./modules/methods"):
      """!
         Construct a ModelFitter object given a model derived from BaseFittingModel and a method
         derived from BaseFittingmethod
         @param model_name name of a fitting model
         @param method_name name of a method 
         @param model_config_filename name of model configuration name 
         @param method_config_filename name of method configuration name 
         @param model_config_dir base configuration directory for models
         @param method_config_dir base configuration directory for methods
         @param fitting_config_dir configuration directory where extra fitting configuration reside
         @param model_module_dir base directory whhere the modules for specific fitting models are 
                located
         @param method_module_dir base directory whhere the modules for specific fitting models are 
                located
      """      

      # --- Properties
      self._model_config_dir = model_config_dir             # model configuration directory name   
      self._method_config_dir =  method_config_dir          # model configuration directory name  
      self._fitting_config_dir = fitting_config_dir         # fitting configuration directory name  
      self._model_module_dir = model_module_dir             # module directory for models   
      self._method_module_dir = model_module_dir            # module directory for methods  

      if model_name is None:
         raise(ModelFitter.InvalidArgument(
                "Invalid model name", model_name))
      if method_name is None:
         raise(ModelFitter.InvalidArgument(
                "Invalid method name", method_name))

      if model_config_filename is not None and \
         self._file_exists(os.path.join(model_config_dir, model_config_filename)):
            self._model_config_filename  = model_config_filename  # modele configuration file name
      else:
         self._model_config_filename  = model_name + ".cfg"

      if method_config_filename is not None and \
         self._file_exists(os.path.join(method_config_dir, method_config_filename)):
            self._method_config_filename  = method_config_filename  # method configuration file name
      else:
         self._method_config_filename  = method_name + ".cfg"
      

      self._model  = self._setup_model(model_name, self._model_config_filename, 
                                       model_config_dir, model_module_dir) 

      self._method = self._setup_method(method_name, 
                                        self._method_config_filename, 
                                        method_config_dir, fitting_config_dir, method_module_dir) 

      self._fitted_params = None    # fitted parameters: populated when fitting is complete
      self._fitting_status = self._method.fitting_status   # method-independent status

      self._version = "0.5.8"

      # --- Private atttributes

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~   

   # --- Getters  

   @property
   def model(self):
      """! @return the model to fit as an object derived from BaseFittingModel """
      return self._model

   @property
   def method(self):
      """! @return the fitting method to apply as an object derived from BaseFittingMethod """
      return self._method

   @property
   def model_config_filename(self):
      """! @return the base configuration file name for models. """
      return self._model_config_filename

   @property
   def method_config_filename(self):
      """! @return the base configuration file name for methods. """
      return self._method_config_filename

   @property
   def model_config_dir(self):
      """! @return the base configuration directory for models. """
      return self._model_config_dir

   @property
   def method_config_dir(self):
      """! @return the base configuration directory for methods. """
      return self._method_config_dir

   @property
   def fitting_config_dir(self):
      """! @return the base fitting configuration directory. """
      return self._fitting_config_dir

   @property
   def model_module_dir(self):
      """! @return the base module directory for models. """
      return self._model_module_dir

   @property
   def method_module_dir(self):
      """! @return the base module directory for methods. """
      return self._method_module_dir

   @property
   def fitted_params(self):
      """! @return the fitted values after fit() has returned. """
      return self.method._fitted_params

   @property
   def fitting_status(self):
      """! 
         @return method-independent fitting status as a BaseFittingMethod.FittingStatus object once
                 fitting is complete
      """
      return self.method.fitting_status

   @property
   def diag_options(self):
      """! 
         @return the method diagnostic options
      """
      return self.method.diag_options


   @property
   def fitting_info(self):
      """! 
         @return method-specific fitting information as a BaseFittingInfo object once fitting is 
         complete 
      """
      return self.method.fitting_info

   @property
   def version(self):
      """! @return the version of the ModelFitter class """
      return self._version


   # --- Setters  

   @model.setter 
   def model(self, model):
      """! 
         Secify the model to fit 
         @param model a model object derived from BaseFittingModel   
      """
      self._model = model 

   @method.setter 
   def method(self, method):
      """! 
         Specify the fitting method to apply 
         @param method a method object derived from BaseFittingMethod   
      """
      self._method = method 

   # ~~~~~~~~~~~~~~
   # Public Methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def fit(self, args):
      """! 
         Fit data using input arguments and the specified parametric model & fitting method. 
         Return a BaseFittingMethod.FittingStatus object that contains the outcome of the fit
         @param args a list of arguments that can be exploited for performing the fit
         @return FittingStatus object 
      """

      if len(self.model.params) == 0: 
         # Nothing to fit
         raise(ModelFitter.FittingError(
             "The model does not have any parameter to fit model {0}.".format(self.model.name)))

      # Invoke method <method> to fit model <model>    
      self._fitting_status = self.method.fit(self.model, args)

      return self._fitting_status

   # -----------------------------------------------------------------------------------------------
   def get_version(self):
      return self.version

   # ~~~~~~~~~~~~~~~
   # Private Methods 
   # ~~~~~~~~~~~~~~~
   
   # -----------------------------------------------------------------------------------------------
   def _import_module(self, module_name, module_dir):
      """! Find and load a module @c module_name from a directory @c module_dir """

      try:

         file_obj, filename, data = imp.find_module(module_name, [module_dir])
         return imp.load_module(module_name, file_obj, filename, data)

      except:
         raise(ModelFitter.FittingError("Could not load module {0} from {1}".format(
                                                                       module_name, module_dir)))

   # -----------------------------------------------------------------------------------------------
   def _setup_method(self, method_name, method_config_filename,
                           method_config_dir, fitting_config_dir, method_module_dir):
      """! Load and instantiate a FittingMethod object for method @c method_name """   

      module = self._import_module(method_name, method_module_dir)

      return module.FittingMethod(method_config_dir=method_config_dir, 
                                  fitting_config_dir=fitting_config_dir,
                                  method_config_filename=method_config_filename)

   # -----------------------------------------------------------------------------------------------
   def _setup_model(self, model_name, model_config_filename, 
                          model_config_dir, model_module_dir):
      """! Load and instantiate a FittingModel object for module @c model_name """   

      model = None 
      if model_name is not None:
         module = self._import_module(model_name, model_module_dir)
         model = module.FittingModel(config_dir=model_config_dir, 
                                     config_filename=model_config_filename)

      return model

   # -----------------------------------------------------------------------------------------------
   def _default_objective_func(self):
      print("Multifit - In default dummy objective function. Please specify a valid objective function.")


   # -----------------------------------------------------------------------------------------------
   def _file_exists(self, filepath):
      """! 
         Safe way of checking if a file (not a directory) does exist. 
         @param filepath file path to chech for existence
         @retval True if file exists
         @retval False otherwise
      """
      try:
        with open(filepath) as file:
            return True
      except IOError as e:
        return False

   # ~~~~~~~~~~
   # Exceptions 
   # ~~~~~~~~~~   

#   # ------------------------------------------------------------------------------
#   class ModelFitterException(Exception):
#      def __str__(self):
#         return "MultiFit *** ERROR ***: {0}".format(self.message)

   # ------------------------------------------------------------------------------
   class InvalidArgument(Exception):
      """! 
         Exception thrown in case of invalid argument(s) supplied
      """
      def __init__(self, message, argument):
         """!
            @param message message to print
            @param argument argument passed
         """
         self._message = message
         self._argument = argument

      def __str__(self):
         """! String representation of the InvalidArgument class """
         return "MultiFit *** ERROR ***: {0} - Invalid argument: {1} .".format(
                                                                    self._message, self._argument)

   # ------------------------------------------------------------------------------
   class ConfigError(Exception):
      """! 
         Exception thrown in case of the configuration file specified could not be found on disk
      """
      def __init__(self, filepath):
         """!
            @param filepath file path of the configuration file 
         """
         self._filepath = filepath

      def __str__(self):
         """! String representation of the ConfigError class """
         return "MultiFit *** ERROR ***: configuration file {0} not found.".format(self._filepath)

   # ------------------------------------------------------------------------------
   class NotImplementedError(Exception):
      """! 
         Exception thrown when a method that is not yet implented has been invoked
      """
      def __init__(self, method_name):
         Exception.__init__(self)

         """!
            @param method_name method name that has been invoked
         """
         self._method_name = method_name

      def __str__(self):
         """! String representation of the NotImplementedError class """
         return "MultiFit *** ERROR ***: method {0} must be implemented in derived classes.".format(
                                                                                       method_name) 

   # -----------------------------------------------------------------------------------------------
   class FittingError(Exception):
      """! 
         Exception thrown when a fitting error has occurred
      """
      def __init__(self, msg):
         """!
            @param msg error message 
         """
         self._msg = msg
         self._exc_info = sys.exc_info()[1]

      def __str__(self):
         """! String representation of the FittingError class """
         if self._exc_info is not None:
            return "MultiFit *** ERROR ***: {0} ({1})".format(self._msg, self._exc_info)
         else:
            return "MultiFit *** ERROR ***: {0}".format(self._msg)
   

   
