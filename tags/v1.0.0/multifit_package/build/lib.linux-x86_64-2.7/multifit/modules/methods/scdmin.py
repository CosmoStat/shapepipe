"""! 
   @package multifit.scdmin Fitting using SCDM
   @author Marc Gentile
   @file scdmin.py
   Fitting using SCDM
""" 

# --- Python imports
import math
import numpy
import numpy.linalg

# -- External imports
import scdm

# --- MultiFit
from multifit import *
from profilehooks import *

# --------------------------------------------------------------------------------------------------
class FittingMethod(BaseFittingMethod):

   """! 
      Fitting using the @c scdm coordinate fitting
      @note extends parent BaseFittingMethod class
   """

   def __init__(self, method_config_dir='./config/methods', 
                      fitting_config_dir='./config/fitting',
                      method_config_filename="scdmin.cfg"):
      """! 
         Construct a FittingMethod object 
         @param method_config_dir configuration directory of the method
         @param fitting_config_dir directory of extra configuration files for the method
         @param method_config_filename configuration filename of the method
      """     
 
      BaseFittingMethod.__init__(self, name="scdmin", 
                                       method_config_dir=method_config_dir,
                                       fitting_config_dir=fitting_config_dir,
                                       method_config_filename=method_config_filename)
      
      # --- Public properties  

#      # --- Objective function
#      if self.objective_func is None:
#         self.objective_func = self.chi2_objective_func  # template chi2 function by default

      # --- Diagnostic Options
      self._fitting_options = FittingMethod.FittingOptions(self)

      scdm_fitting_options =  scdm.FittingOptions(max_iter=self.fitting_options.max_iter,
                                                  max_fev=self.fitting_options.max_fev,
                                                  max_func_tol=self.fitting_options.max_func_tol,
                                                  mini_type=self.fitting_options.minimization_type)

      # --- Diagnostic Options
      self._diag_options = FittingMethod.DiagOptions(self)
      scdm_diag_options = scdm.DiagOptions(self.diag_options.diag_flags, 
                                           self.diag_options.base_output_dir)
      scdm_diag_options.log_output_dir = self.diag_options.log_output_dir
      scdm_diag_options.plot_output_dir = self.diag_options.plot_output_dir
      scdm_diag_options.stats_output_dir = self.diag_options.stats_output_dir
      scdm_diag_options.log_file_name = self.diag_options.log_file_name

      # --- Extra parameter information for fitting
      self.__model_params_dico = self._prepare_sdcm_parameters(self.name, fitting_config_dir)

      # --- Minimizer (private)
      self.__minimizer = scdm.Minimizer(fitting_options=scdm_fitting_options, 
                                        diag_options=scdm_diag_options)

   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~
   
   # --- Getters

   @property
   def fitting_options(self):
      """! @return the fitting options. """
      return self._fitting_options

   @property
   def diag_options(self):
      """! @return the diagnostic options. """
      return self._diag_options

   # --- Setters

   @diag_options.setter 
   def diag_options(self, diag_options):
      """! 
         Set the the diagnostic options.
         @note refer to the corresponding @c SCDM documentation
      """

      scdm_diag_options = scdm.DiagOptions(diag_options.diag_flags, 
                                           diag_options.base_output_dir)
      scdm_diag_options.log_output_dir = diag_options.log_output_dir
      scdm_diag_options.plot_output_dir = diag_options.plot_output_dir
      scdm_diag_options.stats_output_dir = diag_options.stats_output_dir
      scdm_diag_options.log_file_name = diag_options.log_file_name
      self.__minimizer.diag_options = scdm_diag_options


   @fitting_options.setter 
   def fitting_options(self, fitting_options):
      """! 
         Set the the fitting options.
         @note refer to the corresponding @c SCDM documentation
      """
      scdm_fitting_options =  scdm.FittingOptions(fitting_options.max_iter,
                                                  fitting_options.max_fev,
                                                  fitting_options.max_func_tol,
                                                  fitting_options.minimization_type)
      scdm.Minimizer.fitting_options = scdm_FittingOptions

   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   #@profile(immediate=True)
   def fit(self, model, args):
      """! 
         Fit a model @c model using input arguments @c args
         @param model an object derived from BaseFittingModel specifying the parametric model to fit   
         @param args a list of arguments that can be exploited for performing the fit
         @return a multifit.multifit.BaseFittingMethod.FittingStatus object
      """

      if __debug__:
         if self.objective_func is None:      
            raise(ModelFitter.FittingError(
                "An objective function must be specified to fit model {0}.".format(model.name)))

      # --- Setup minimizer object: dynamically set guess values. 
      #     Make sure the original parameter order is preserved

      if model.name in self.__model_params_dico:

         # --- Set extra parameter fitting data
         for param in model.params:

            # --- Guess values
            self.__model_params_dico[model.name]["params"][param.name].guess_value =\
                                                                              param.guess_value

            # --- Check whether the parameter should vary (i.e. be fitted) or not 
            if not self.__model_params_dico[model.name]["params"][param.name].is_constant:

               # --- Tell the minimizer if the parameter is constant or variying
               self.__model_params_dico[model.name]["params"][param.name].is_constant =\
                   param.bounds[0] is not None and param.bounds[0]==param.bounds[1]

            # --- Target params
            if param.name in self.__model_params_dico[model.name]["target_params"]:

               # --- Target guess values
               self.__model_params_dico[model.name]["target_params"][param.name].guess_value =\
                                                                                 param.guess_value

               # --- Tell the minimizer if the parameter is constant or variying
               if not self.__model_params_dico[model.name]["target_params"][param.name].is_constant:
                  self.__model_params_dico[model.name]["target_params"][param.name].is_constant =\
                     param.bounds[0]==param.bounds[1] and param.bounds[0] is not None

         # --- Initialize SCDM parameter list 
         self.__minimizer.params = scdm.Params(
            [self.__model_params_dico[model.name]["params"][p.name] for p in model.params])

         self.__minimizer.target_params = scdm.Params(
            [self.__model_params_dico[model.name]["target_params"][p.name] \
                  for p in model.params \
                   if p.name in self.__model_params_dico[model.name]["target_params"]])            

         #print "**MFIT: fit() target params", [p for p in self.__minimizer.target_params.param_list]

      else:

         # --- No extra parameter fitting information
         self.__minimizer.params = scdm.Params( 
                       [scdm.Param(p.name, p.guess_value, p.bounds) for p in model.params] )         

      # --- Numerical & out-of-bound handlers at parameter level
      for scdm_param in self.__minimizer.params.param_list:

         # --- Numerical & out-of-bound handlers
         if model.out_of_bound_handler is not None:
            scdm_param.out_of_bound_handler = model.out_of_bound_handler
         if model.numerical_error_handler is not None:
            scdm_param.numerical_error_handler = model.numerical_error_handler

         #print "*** MFIT: scdm_param.out_of_bound_handler:", scdm_param.out_of_bound_handler

      self.__minimizer.args = args
      self.__minimizer.func = self.objective_func

      #print "***SCDM: minimizer.target_params:", [p.name for p in self.__minimizer.target_params.param_list] 

      # --- Invoke minimizing routine
      [self._fitted_params, fitting_results] = self.__minimizer.minimize()

      # --- FittingInfo 
      self._fitting_info = FittingMethod.FittingInfo(self, 
                                                     self.__minimizer.fitting_status.errno,
                                                     fitting_results.nb_iter,
                                                     fitting_results.get_nb_func_evals(),
                                                     self.__minimizer.fitting_status.msg,
                                                     self.__minimizer.fitting_status.warn_msgs, 
                                                     fitting_results.fitting_data,
                                                     fitting_results.get_custom_data(),
                                                     self._fitted_params) 

      # --- Fitting status
      self.fitting_status = FittingMethod.FittingStatus(self) 
      self.fitting_status.msg = self.__minimizer.fitting_status.msg
      self.fitting_status.warn_msgs = self.__minimizer.fitting_status.warn_msgs

      self.fitting_status.custom_data = fitting_results.get_custom_data()
      self.fitting_status.residuals = self._fitting_info.residuals 

      errno = self.__minimizer.fitting_status.errno

      if errno >= 0:
         self.fitting_status.errno = errno        
      elif errno == -2 or errno == -1:
         # --- Maximum of iterations or function evaluations reached => still consider as a success
         self.fitting_status.errno = 5
      else:
         # Failure
         self.fitting_status.errno = -1

      return self.fitting_status


   # -----------------------------------------------------------------------------------------------
   def chi2_objective_func(self, params, fitting_func, model, obs_data, sigma_noise_sq, args):
      """!
         Chi2 template objective function 
      """

      model_data = fitting_func(params, model, *args)
      residuals = (obs_data - model_data)**2 / sigma_noise_sq

#      print "**** CHI2: residuals:", numpy.sum(residuals), "Sum Obs:", numpy.sum(obs_data), "Sum Data:", numpy.sum(model_data)

#      print "**** CHI2:params:", params, "residuals:", numpy.sum(residuals),\
#            "mean residuals:", numpy.mean(residuals), "Sum Obs:", numpy.sum(obs_data),\
#            "Sum Data:", numpy.sum(model_data)

      return [numpy.sum(residuals), [model_data, residuals]] 

   # -----------------------------------------------------------------------------------------------
   def chi2_r_objective_func(self, params, 
                                   fitting_func, model, obs_data, sigma_noise_sq, rlambda, args):
      """!
         Chi2 template objective function 
      """

      model_data = fitting_func(params, model, *args)
      if rlambda > 0:
         print rlambda, numpy.linalg.norm(numpy.asarray(params)), rlambda * numpy.linalg.norm(numpy.asarray(params))
         residuals = (obs_data - model_data)**2 / sigma_noise_sq + \
                     rlambda * numpy.linalg.norm(numpy.asarray(params))
      else:
         residuals = (obs_data - model_data)**2 / sigma_noise_sq

      return [numpy.sum(residuals), [model_data, residuals]] 

   # -----------------------------------------------------------------------------------------------
   #@profile(immediate=True)
   def chi2_w_objective_func(self, params, fitting_func, model, obs_data, sigma_noise_sq, *args):
      """!
         Chi2 template objective function: wavelet transform JLS
      """

      alphas = args[-4]
      model_data = fitting_func(params, model, *args)
#      residuals = (alphas * obs_data - model_data)**2 / sigma_noise_sq
      residuals = alphas * (obs_data - model_data)**2 / sigma_noise_sq

#      print "**** CHI2:params:", params, "residuals:", numpy.sum(residuals),\
#            "mean residuals:", numpy.mean(residuals), "Sum Obs:", numpy.sum(obs_data),\
#            "Sum Data:", numpy.sum(model_data)

      return [numpy.sum(residuals), [model_data, residuals]] 

   # -----------------------------------------------------------------------------------------------
   def get_chi2_objective_func(self):
      """! @return the chi squared template objective function. """
      return self.chi2_objective_func

   # -----------------------------------------------------------------------------------------------
   def get_chi2_r_objective_func(self):
      """! @return the chi squared template objective function: regularized """
      return self.chi2_r_objective_func

   # -----------------------------------------------------------------------------------------------
   def get_chi2_w_objective_func(self):
      """! @return the chi squared template objective function: wavelet transform """
      return self.chi2_w_objective_func

   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _prepare_sdcm_parameters(self, method_name, fitting_config_dir):

      model_param_dico = {}

      config_full_path = os.path.join(fitting_config_dir, "{0}_fitting.cfg".format(method_name))
      self.fitting_config = SConfig(config_full_path)

      model_section_names = self.fitting_config.get_subsections(recurse=False)
      for model_section_name in model_section_names:

         param_list = []
         target_param_list = [] 
         model_param_dico[model_section_name.lower()] = {}
         model_param_dico[model_section_name.lower()]["params"] = {}
         model_param_dico[model_section_name.lower()]["target_params"] = {}

         subsection_name = "{0}.PARAM_FITTING_INFO".format(model_section_name)
         param_subsections = self.fitting_config.get_subsections(subsection_name)
         for param_subsection in param_subsections:

            param_name = self.fitting_config.get_as_string("NAME", param_subsection)
            param_range = self.fitting_config.get_as_list("RANGE", param_subsection)

            param_step_size = self.fitting_config.get_as_list("STEP_SIZE", param_subsection)
            param_step_range  = self.fitting_config.get_as_list("STEP_RANGE", param_subsection)
            param_step_factor = self.fitting_config.get_as_float("STEP_FACTOR", param_subsection)
            #param_skip_func_err = eval(self.fitting_config.get_as_string("SKIP_FUNC_ERR", 
            #                                                        param_subsection))
            #param_skip_count = self.fitting_config.get_as_int("SKIP_COUNT", param_subsection)
            param_max_func_err = self.fitting_config.get_as_float("MAX_FUNC_ERR", param_subsection)

            if self.fitting_config.has_key("IS_TARGET", param_subsection):
               param_is_target   = self.fitting_config.get_as_boolean("IS_TARGET", param_subsection)
            else:
               param_is_target   = False
            if self.fitting_config.has_key("IS_CONSTANT", param_subsection):
               param_is_constant = self.fitting_config.get_as_boolean(
                                                                 "IS_CONSTANT", param_subsection)
            else:
               param_is_constant = (param_range[0]==param_range[1] and \
                                    param_range[0] is not None)

            model_param_dico[model_section_name.lower()]["params"][param_name] =\
              scdm.Param(name=param_name, param_range=param_range, 
                         step_size=param_step_size, step_range=param_step_range,
                         step_factor=param_step_factor, max_tol=param_max_func_err,
                         is_constant=param_is_constant)

            if param_is_target:
               model_param_dico[model_section_name.lower()]["target_params"][param_name] =\
                 scdm.Param(name=param_name, param_range=param_range, 
                            step_size=param_step_size, step_range=param_step_range,
                            step_factor=param_step_factor, max_tol=param_max_func_err,
                            is_constant=param_is_constant)
 
         # --- end for param_subsection            

      # --- end for model_section_name            

      return model_param_dico


   # ~~~~~~~~~~~~~
   # Inner Classes
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class FittingOptions(BaseFittingOptions):
      """! 
         Fitting options applicable to @c SCDM
         @note refer to the corresponding @c SCDM documentation 
         @see multifit.multifit.BaseFittingOptions
      """

      def __init__(self, method):
         """!
            Construct a FittingOptions class
            @param method a FittingMethod object to which the FittingOptions options apply
         """
         BaseFittingOptions.__init__(self, method)

         # --- Public properties
         self._max_iter = 10000
         self._max_fev  = 10000
         self._max_func_tol = 1.0e-8
         self._minimization_type = scdm.MinimizationType.minimize_param_delta | \
                                   scdm.MinimizationType.minimize_param_sigma | \
                                   scdm.MinimizationType.minimize_func_value

         # --- Set default options from configuration
         self._init_from_config(method.config)


      def __str__(self):
         """! 
            String representation of the FittingOptions object 
            @return a Python String describing the FittingOptions object
         """
         descr = BaseFittingOptions.__str__(self)
         descr += "- max iter: {0}\n".format(self.max_iter)
         descr += "- max fev: {0}\n".format(self.max_fev)
         descr += "- max func tol: {0}\n".format(self.max_func_tol)
         if self.minimization_type & scdm.MinimizationType.minimize_param_delta == scdm.MinimizationType.minimize_param_delta:
            descr += "- minimization: {0}\n".format("minimize parameter delta")
         elif self.minimization_type & scdm.MinimizationType.minimize_param_sigma == scdm.MinimizationType.minimize_param_sigma:
            descr += "- minimization: {0}\n".format("minimize parameter sigma")
         elif self.minimization_type & scdm.MinimizationType.minimize_func_value == scdm.MinimizationType.minimize_func_value:
            descr += "- minimization: {0}\n".format("minimize function value")
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters 

      @property
      def max_iter(self):
         """! 
            @return the maximum number of iterations allowed 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._max_iter

      @property
      def max_fev(self):
         """! 
            @return the maximum number of function evaluations allowed 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._max_fev

      @property
      def max_func_tol(self):
         """! 
            @return the maximum allowed error on residuals between model and data
            @note refer to the corresponding @c SCDM documentation
         """
         return self._max_func_tol

      @property
      def minimization_type(self):
         """! 
            @return the type of minimization requested (stopping criteria)
            @note: refer toi the SCDM documentation for the details
         """
         return self._minimization_type


      # --- Setters 

      @max_iter.setter 
      def max_iter(self, max_iter):
         """! 
            Set the maximum number of iterations allowed 
            @note refer to the corresponding @c SCDM documentation
         """
         self._max_iter = max_iter

      @max_fev.setter 
      def max_fev(self, max_fev):
         """! 
            Set the maximum number of function evaluations allowed 
            @note refer to the corresponding @c SCDM documentation
         """
         self._max_fev = max_fev

      @max_func_tol.setter 
      def max_func_tol(self, max_func_tol):
         """! 
            Set the maximum allowed error on residuals between model and data
            @note refer to the corresponding @c SCDM documentation
         """
         self._max_func_tol = max_func_tol

      @minimization_type.setter 
      def xtol(self, minimization_type):
         """! 
            Set the type of minimization requested (stopping criteria)
            @note refer to the corresponding @c SCDM documentatio
         """
         self._minimization_type = minimization_type


      # ~~~~~~~~~~~~~~~~
      # Private methods
      # ~~~~~~~~~~~~~~~~

      # -----------------------------------------------------------------------------------------------
      def _init_from_config(self, config):

         if not config is None:
            method_section_name = self.method.name.upper()
            method_options_section_name = "{0}.OPTIONS".format(method_section_name)

            if config.has_section(method_section_name):

               # --- Default fitting options
               if config.has_key("max_iter", method_options_section_name):
                  self._max_iter = config.get_as_int("max_iter", method_options_section_name)
               if config.has_key("max_fev", method_options_section_name):
                  self._max_fev  = config.get_as_int("max_fev",  method_options_section_name)
               if config.has_key("max_func_tol", method_options_section_name):
                  self._max_func_tol  = config.get_as_float("max_func_tol", 
                                                            method_options_section_name)

               if config.has_key("mini_type", method_options_section_name):

                  self._minimization_type = eval(config.get_as_string("mini_type",
                                                                      method_options_section_name))

                  #print "*** MFIT: minimization_type:", self._minimization_type

               self._all_minima = False    # not implemented

         else:
            if not self.method.name is None:
               # --- No configuration file found => taking suitable defaults             
               print("MultiFit *** Warning ***: configuration file not found for {0}.".format(
                                                                                 self.method.name))   

   # -----------------------------------------------------------------------------------------------
   class DiagOptions(BaseDiagOptions):
      """! 
         Diagnostic options for SCDM
         @note refer to the corresponding @c SCDM documentation 
         @see multifit.multifit.BaseDiagOptions
      """

      def __init__(self, method):
         """!
            Construct a DiagOptions class
            @param method a FittingMethod object to which the DiagOptions options apply
         """
         BaseDiagOptions.__init__(self, method)



         # --- Initialize options from configuration
         self._init_from_config(method.config)

      def __str__(self):
         """! 
            String representation of the DiagOptions object 
            @return a Python String describing the DiagOptions object
         """
         descr = BaseDiagOptions.__str__(self)
         if self._diag_flags & scdm.DiagFlag.logs == scdm.DiagFlag.logs:
            descr += "- log output enabled\n"
         if self._diag_flags & scdm.DiagFlag.plots == scdm.DiagFlag.plots:
            descr += "- plotting enabled\n"
         if self._diag_flags & scdm.DiagFlag.stats == scdm.DiagFlag.stats:
            descr += "- fitting statistics enabled\n"
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters 

      @property
      def diag_flags(self):
         """! 
            @return the diagnostic flags 
            @see scdm.DiagFlag
            @note refer to the corresponding @c SCDM documentation
         """
         return self._diag_flags


      # --- Setters 

      @diag_flags.setter 
      def diag_flags(self, diag_flags):
         """! 
            Set the diagnostic flags
            @param diag_flags diagnostic flagsoutputing 
            @see scdm.DiagFlag
            @note refer to the corresponding @c SCDM documentation
          """
         self._diag_flags = diag_flags


      # ~~~~~~~~~~~~~~~~
      # Private methods
      # ~~~~~~~~~~~~~~~~

      def _init_from_config(self, config):

         method_section_name = self.method.name.upper()
         method_options_section_name = "{0}.OPTIONS".format(method_section_name)

         if config.has_section(method_section_name):
            if config.has_key("diag_flags", method_options_section_name):
               self._diag_flags = eval(config.get_as_string("diag_flags", 
                                                            method_options_section_name))   

         else:
            if not self.method.name is None:
               print("MultiFit *** Warning ***: configuration file notfound for {0}.".format(
                                                                                 self.method.name))   



   # -----------------------------------------------------------------------------------------------
   class FittingInfo(BaseFittingInfo):

      """! 
          Fitting information available from SCDM, once fitting is complete. 
          @note refer to the corresponding @c SCDM documentation
      """

      def __init__(self, method, 
                         errno, nb_iter, nb_fev, msg, warn_msgs, 
                         fitting_data, custom_data, fitted_params):
         """!
            Construct a FittingInfo class

            @param method a FittingMethod object to which the FittingInfo options apply
            @param errno error code returned by SCDM
            @param nb_iter actual number of iterations
            @param nb_fev actual number of function evaluations
            @param msg message returned by SCDM, if any
            @param warn_msgs warning messages returned by SCDM, if any
            @param fitting_data data collected by SCDM, if any
            @param custom_data SCDM custom data, if any
            @param fitted_params best-fit parameters
         """         
         BaseFittingInfo.__init__(self, method)

         self._errno = errno                    # result code
         self._nb_iter = nb_iter                # number of iterations
         self._nb_fev = nb_fev                  # number of function evaluations
         self._msg = msg                        # message, if any
         self._warn_msgs = warn_msgs            # list of warning messages, if any
         self._fitting_data = fitting_data      # collected fitting data (if enabled)
         self._fitted_params = fitted_params    # best-fit parameters
         self._custom_data = custom_data        # specific additional fitting data, if any
         if not self._custom_data is None and len(self._custom_data) > 0:
            self._model_data = numpy.ravel(custom_data[0][0])
            self._residuals = numpy.ravel(custom_data[0][1])
         else:
            self._model_data = None
            self._residuals = None

      def __str__(self):
         """! 
            String representation of the FittingInfo object 
            @return a Python String describing the FittingInfo object
         """
         descr = BaseFittingInfo.__str__(self)
         descr += "- errno: {0}\n".format(self.errno)
         descr += "- nb iter: {0}\n".format(self.nb_iter)
         descr += "- nb fev: {0}\n".format(self.nb_fev)

         descr += "- msg: {0}\n".format(self.msg)
         if len(self.warn_msgs) > 0:
            descr += "- warn msgs:\n"
            for warn_msg in self.warn_msgs:
               descr += "- {0}\n".format(warn_msg)
         descr += "- best-fit parameters: {0}\n".format(self.fitted_params)
         descr += "- residuals: {0}\n".format(self.residuals)
         descr += "- model data: {0}\n".format(self.model_data)
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters          

      @property
      def errno(self):
         """! 
            @return the error code produced by SCDM 
            @retval ier = 1, 2, 3 or 4 if a solution was found (0: success, < 0: failure)
            @note refer to the corresponding @c SCDM documentation
         """
         return self._errno

      @property
      def nb_iter(self):
         """! 
            @return the number of iterations that were required for the fit 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._nb_iter

      @property
      def nb_fev(self):
         """! 
            @return the number of function evaluations that were required for the fit 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._nb_fev

      @property
      def msg(self):
         """! 
            @return the message returned by SCDM, if any 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._msg

      @property
      def warn_msgs(self):
         """! 
            @return the list of warning messages returned by SCDM, if any 
            @note refer to the corresponding @c SCDM documentation
         """
         return self._warn_msgs

      @property
      def fitting_data(self):
         """! 
            @return data collected by SCDM, if any
            @note refer to the corresponding @c SCDM documentation
         """
         return self._fitting_data

      @property
      def custom_data(self):
         """! 
            @return extra collected by the model fitting function, if any
            @note refer to the corresponding @c SCDM documentation
         """
         return self._custom_data

      @property
      def fitted_params(self):
         """! 
            @return the best-fit parameters
            @note refer to the corresponding @c SCDM documentation
         """
         return self._fitted_params

      @property
      def residuals(self):
         """! 
            @return the final residuals
            @note refer to the corresponding @c SCDM documentation
         """
         return self._residuals

      @property
      def model_data(self):
         """! 
            @return the final model data
            @note refer to the corresponding @c SCDM documentation
         """
         return self._model_data

# --- EOF scdm.py 
