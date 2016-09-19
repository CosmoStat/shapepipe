"""!
   @package multifit.kmpfit Fitting using kmpfit (Levenberg-Marquardt)
   [https://www.astro.rug.nl/software/kapteyn/kmpfit.html)]
   @author Marc Gentile
   @file kmpfit.py
   Fitting using kmpfit (Levenberg-Marquardt)
""" 

# --- Python imports
import math
import numpy
import numpy.linalg

# -- External imports
from kapteyn import kmpfit

# --- MultiFit
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingMethod(BaseFittingMethod):

   """! 
      Fitting using the @c kmpfit least-squares minimizer
      @note extends parent BaseFittingMethod class
   """

   def __init__(self, method_config_dir='./config/methods', 
                      fitting_config_dir='./config/fitting',
                      method_config_filename="kmpfit.cfg"):
      """! 
         Construct a FittingMethod object 
         @param method_config_dir configuration directory of the method
         @param fitting_config_dir directory of extra configuration files for the method
         @param method_config_filename configuration filename of the method
      """     
 
      BaseFittingMethod.__init__(self, name="kmpfit", 
                                       method_config_dir=method_config_dir,
                                       fitting_config_dir=fitting_config_dir,
                                       method_config_filename=method_config_filename)
      
      # --- Public properties  

      # --- Objective function
      if self.objective_func is None:
         self.objective_func = self.chi2_objective_func  # template chi2 function by default

      # --- Diagnostic Options
      self._fitting_options = FittingMethod.FittingOptions(self)

      # --- Diagnostic Options
      self._diag_options = FittingMethod.DiagOptions(self)

      # --- Extra parameter information for fitting
      self.__model_params_dico = self._prepare_kmpfit_parameters(self.name, fitting_config_dir)

      # --- Minimizer (private)
      self.__minimizer = kmpfit.Fitter(residuals=self.objective_func,
                                       ftol=self._fitting_options.ftol,
                                       xtol=self._fitting_options.xtol,
                                       gtol=self._fitting_options.gtol,
                                       epsfcn=self._fitting_options.epsfcn,
                                       stepfactor=self._fitting_options.stepfactor,
                                       covtol=self._fitting_options.covtol,
                                       maxiter=self._fitting_options.maxiter,
                                       maxfev=self._fitting_options.maxfev)

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
         @note refer to the corresponding @c kmpfit documentation
      """

      kmpfit_diag_options = kmpfit.DiagOptions(diag_options.diag_flags, 
                                           diag_options.base_output_dir)
      kmpfit_diag_options.log_output_dir = diag_options.log_output_dir
      kmpfit_diag_options.plot_output_dir = diag_options.plot_output_dir
      kmpfit_diag_options.stats_output_dir = diag_options.stats_output_dir
      kmpfit_diag_options.log_file_name = diag_options.log_file_name
      self.__minimizer.diag_options = kmpfit_diag_options


   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def fit(self, model, args):
      """! 
         Fit a model @c model using input arguments @c args
         @param model an object derived from BaseFittingModel specifying the parametric model to fit   
         @param args a list of arguments that can be exploited for performing the fit
         @return a multifit.multifit.BaseFittingMethod.FittingStatus object
         @note in case of sucess the @c errno code in 
               multifit.multifit.BaseFittingMethod.FittingStatus will be positive. 
               The actual error code returned by the underlying minimizer is accessible through
               the multifit.multifit.BaseFittingMethod.FittingInfo object
      """

      if __debug__:
         if self.objective_func is None:
            raise(ModelFitter.FittingError(
                "An objective function must be specified to fit model {0}.".format(model.name)))

      # --- Setup minimizer object: dynamically set guess values. 
      #     Make sure the original parameter order is preserved

      if model.name in self.__model_params_dico:

         #print "param dico:", self.__model_params_dico[model.name]["parinfo"]

         # --- Set extra parameter fitting data
         for param in model.params:

            # --- Check whether the parameter should vary (i.e. be fitted) or not 
            param_is_fixed = param.bounds[0]==param.bounds[1] and param.bounds[0] is not None
            self.__model_params_dico[model.name]["parinfo"][param.name]["fixed"] = param_is_fixed

            # --- Override the parameter range
#            if param.bounds[0] is None and param.bounds[1] is None:
#               param.bounds = None
            self.__model_params_dico[model.name]["parinfo"][param.name]["limits"] = param.bounds

         # --- Initialize kmpfit parameter info (parinfo) list 
         kmpfit_parinfo = [self.__model_params_dico[model.name]["parinfo"][p.name] \
                                                                            for p in model.params]


         # --- Parameter guess values
         kmpfit_guess_values = model.get_param_values()  

      else:

         # --- No extra parameter fitting information
         kmpfit_parinfo = [{} for p in model.params]

         # --- Parameter guess values
         kmpfit_guess_values = model.get_param_values()
      
      #print "*** kmpfit: parinfo:", kmpfit_parinfo
      #print "*** guess values:", kmpfit_guess_values

      # --- Invoke minimizing routine
      self.__minimizer.residuals = self.objective_func
      self.__minimizer.parinfo = kmpfit_parinfo
      self.__minimizer.data = args   

      self.fitting_status.errno = -1   

      try:      

         self.__minimizer.fit(params0 = kmpfit_guess_values)

      except Exception as detail:

         print "*** kmpfit fitting error: ", detail

      # --- Record best-fit parameters
      self._fitted_params = self.__minimizer.params

      # --- FittingInfo 
      self._fitting_info = FittingMethod.FittingInfo(self, 
                                                     self.__minimizer.version,
                                                     self.__minimizer.status,
                                                     self.__minimizer.message,
                                                     self.__minimizer.niter,
                                                     self.__minimizer.nfev,
                                                     self.__minimizer.nfree,
                                                     self.__minimizer.npegged,
                                                     self.__minimizer.dof,
                                                     self.__minimizer.xerror,
                                                     self.__minimizer.stderr,
                                                     self.__minimizer.covar,
                                                     self.__minimizer.chi2_min,
                                                     self.__minimizer.rchi2_min,
                                                     self.__minimizer.resid)

      # --- Fitting status
      self.fitting_status = FittingMethod.FittingStatus(self) 
      self.fitting_status.msg = self.__minimizer.message
      self.fitting_status.warn_msgs = []
      self.fitting_status.residuals = self.__minimizer.resid
      self.fitting_status.custom_data = None

      ier = self.__minimizer.status
      if ier > 0:
         self.fitting_status.errno = ier         
      else:
         # Failure
         self.fitting_status.errno = -1

      #print "info:", self._fitting_info
      #print "status:", self.fitting_status

      return self.fitting_status


   # -----------------------------------------------------------------------------------------------
   def chi2_objective_func(self, params, data):
      """!
         Chi2 template objective function 
      """
      fitting_func, model, obs_data, sigma_noise_sq, args = data

      #print "params:", params, "Sigma:", sigma_noise_sq

      #residuals = (obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq
      #print "**** CHI2: residuals:", numpy.sum(residuals)

      return numpy.ravel((obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq)

   # -----------------------------------------------------------------------------------------------
   def chi2_r_objective_func(self, params, data):
      """!
         Chi2 template objective function 
      """

      data = fitting_func, model, obs_data, sigma_noise_sq, rlambda, args

      #residuals = (obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq
      #print "**** CHI2: residuals:", numpy.sum(residuals)

      if rlambda > 0:
         rterm = rlambda * numpy.linalg.norm(numpy.asarray(params))
         return numpy.ravel((obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq) + \
                rlambda * numpy.linalg.norm(numpy.asarray(params))  
      else:
         return numpy.ravel((obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq) 


   # -----------------------------------------------------------------------------------------------
   def get_chi2_objective_func(self):
      """! @return the chi squared template objective function. """
      return self.chi2_objective_func



   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _prepare_kmpfit_parameters(self, method_name, fitting_config_dir):

      model_param_dico = {}

      config_full_path = os.path.join(fitting_config_dir, "{0}_fitting.cfg".format(method_name))
      self.fitting_config = SConfig(config_full_path)

      model_section_names = self.fitting_config.get_subsections(recurse=False)
      for model_section_name in model_section_names:

         param_list = []
         target_param_list = [] 
         model_sname = model_section_name.lower()
         
         model_param_dico[model_sname] = {}
         model_param_dico[model_sname]["parinfo"] = {}

         subsection_name = "{0}.PARAM_FITTING_INFO".format(model_section_name)
         param_subsections = self.fitting_config.get_subsections(subsection_name)
         keyword_list = ['fixed', 'limits', 'step', 'side', 'deriv_debug'] # order is important

         for param_subsection in param_subsections:

            param_name = self.fitting_config.get_as_string("NAME", param_subsection)
            param_range = eval(self.fitting_config.get_as_string("RANGE", param_subsection))

            #print "Param name:", param_name, "Range:", param_range

            if param_range is not None:
               param_range = self.fitting_config.get_as_list("RANGE", param_subsection)
#               if param_range[0] is None and param_range[1] is None:
#                  param_range = None

            param_step_size = eval(self.fitting_config.get_as_string("STEP_SIZE", param_subsection))
            if param_step_size is not None:
               param_step_size = self.fitting_config.get_as_float("STEP_SIZE", param_subsection)
               
            param_is_constant = self.fitting_config.get_as_boolean(
                                                        "IS_FIXED", param_subsection)

            if param_range is not None:   
               param_is_constant = param_is_constant or \
                                   (param_range[0]==param_range[1] and param_range[0] is not None)
            
            param_deriv_type  = self.fitting_config.get_as_float("DERIV_TYPE", param_subsection)
            param_deriv_debug = self.fitting_config.get_as_boolean("DERIV_DEBUG", param_subsection)

            model_param_dico[model_sname]["parinfo"][param_name] = \
                dict(zip(keyword_list,
                         [param_is_constant, param_range, param_step_size, \
                          param_deriv_type, param_deriv_debug]))

         # --- end for param_subsection            

      # --- end for model_section_name            

      return model_param_dico


   # ~~~~~~~~~~~~~
   # Inner Classes
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class FittingOptions(BaseFittingOptions):
      """! 
         Fitting options applicable to the @c kmpfit least-squares minimizer
         @note refer to the corresponding @c kmpfit documentation
         @see multifit.multifit.BaseFittingOptions
      """

      def __init__(self, method):
         """!
            Construct a FittingOptions class
            @param method a FittingMethod object to which the FittingOptions options apply
         """
         BaseFittingOptions.__init__(self, method)

         self._ftol = 1.0e-10          # Relative chi2 convergence criterium. Default: 1e-10
         self._xtol = 1.0e-10          # Relative parameter convergence criterium. Default: 1e-10
         self._gtol = 0.0              # Orthogonality convergence criterium. Default: 1e-10
         self._epsfcn  = 2.2204460e-16 # Finite derivative step size. Default: 2.2204460e-16
         self._covtol  = 1.0e-14       # Range tolerance for covariance calculation. Default: 1e-14
         self._stepfactor = 100.0      # Initial step bound. Default: 100.0
         self._maxfev  = 0             # Maximum number of function evaluations, 0: no limit
         self._maxiter = 200           # Max. number of iterations. Default: 200

         # --- Public properties  
   
         # --- Set default options from configuration
         self._init_from_config(method.config)

      def __str__(self):
         """! 
            String representation of the FittingOptions object 
            @return a Python String describing the FittingOptions object
         """
         descr = BaseFittingOptions.__str__(self)
         descr += "- ftol: {0}\n".format(self.ftol)
         descr += "- xtol: {0}\n".format(self.xtol)
         descr += "- gtol: {0}\n".format(self.gtol)
         descr += "- epsfcn: {0}\n".format(self.epsfcn)
         descr += "- covtol: {0}\n".format(self.covtol)
         descr += "- stepfactor: {0}\n".format(self.stepfactor)
         descr += "- maxfev: {0}\n".format(self.maxfev)
         descr += "- maxiter: {0}\n".format(self.maxiter)
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters 

      @property
      def ftol(self):
         """! 
            @return the @c ftol option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._ftol

      @property
      def gtol(self):
         """! 
            @return the @c gtol option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._gtol

      @property
      def xtol(self):
         """! 
            @return the @c xtol option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._xtol

      @property
      def maxfev(self):
         """! 
            @return the @c maxfev option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._maxfev

      @property
      def maxiter(self):
         """! 
            @return the @c maxiter option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._maxiter

      @property
      def epsfcn(self):
         """! 
            @return the @c epsfcn option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._epsfcn

      @property
      def covtol(self):
         """! 
            @return the @c covtol option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._covtol

      @property
      def stepfactor(self):
         """! 
            @return the @c stepfactor option value of @c kmpfit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._stepfactor


      # --- Setters 

      @ftol.setter 
      def ftol(self, ftol):
         """! 
            Set the @c ftol option value of @c kmpfit
            @param ftol the @c gtol option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._ftol = ftol

      @gtol.setter 
      def gtol(self, gtol):
         """! 
            Set the @c gtol option value of @c kmpfit
            @param gtol the @c gtol option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._gtol = gtol

      @xtol.setter 
      def xtol(self, xtol):
         """! 
            Set the @c xtol option value of @c kmpfit
            @param xtol the @c xtol option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._xtol = xtol

      @maxfev.setter 
      def maxfev(self, maxfev):
         """! 
            Set the @c maxfev option value of @c kmpfit
            @param maxfev the @c maxfev option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._maxfev = maxfev

      @maxiter.setter 
      def maxiter(self, maxiter):
         """! 
            Set the @c maxfev option value of @c kmpfit
            @param maxiter the @c maxiter option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._maxiter = maxiter

      @epsfcn.setter 
      def epsfcn(self, epsfcn):
         """! 
            Set the @c epsfcn option value of @c kmpfit
            @param epsfcn the @c epsfcn option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._epsfcn = epsfcn

      @covtol.setter 
      def covtol(self, covtol):
         """! 
            Set the @c factor option value of @c kmpfit
            @param covtol the @c covtol option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._covtol = covtol

      @stepfactor.setter 
      def stepfactor(self, stepfactor):
         """! 
            Set the @c stepfactor option value of @c kmpfit
            @param stepfactor the @c stepfactor option value
            @note refer to the corresponding @c kmpfit documentation
         """
         self._stepfactor = stepfactor


      # ~~~~~~~~~~~~~~~~
      # Private methods
      # ~~~~~~~~~~~~~~~~

      def _init_from_config(self, config):

         if not config is None:
            method_section_name = self.method.name.upper()
            model_options_section_name = "{0}.OPTIONS".format(method_section_name)
            if config.has_section(method_section_name):

               # --- Default configuration
               if config.has_key("ftol", model_options_section_name):
                  self.ftol = config.get_as_float("ftol", model_options_section_name)
               if config.has_key("gtol", model_options_section_name):
                  self.gtol = config.get_as_float("gtol", model_options_section_name)
               if config.has_key("xtol", model_options_section_name):
                 self.xtol = config.get_as_float("xtol", model_options_section_name)
               if config.has_key("maxfev", model_options_section_name):
                  self.maxfev = config.get_as_int("maxfev", model_options_section_name)
               if config.has_key("maxiter", model_options_section_name):
                  self.maxiter = config.get_as_int("maxiter", model_options_section_name)
               if config.has_key("epsfcn", model_options_section_name): 
                 self.epsfcn = config.get_as_float("epsfcn", model_options_section_name)  
               if config.has_key("covtol", model_options_section_name): 
                  self.covtol = config.get_as_float("covtol", model_options_section_name)      
               if config.has_key("stepfactor", model_options_section_name):
                  self.stepfactor = config.get_as_float("stepfactor", model_options_section_name)
         else:
            if not self.method.name is None:
               # --- No configuration file found => taking suitable defaults             
               print("MultiFit *** Warning ***: configuration file not found for {0}.".format(
                                                                                 self.method.name))

 

   # -----------------------------------------------------------------------------------------------
   class DiagOptions(BaseDiagOptions):
      """! 
         Diagnostic options for kmpfit
         @note refer to the corresponding @c kmpfit documentation 
         @see multifit.multifit.BaseDiagOptions
      """

      def __init__(self, method):
         """!
            Construct a DiagOptions class
            @param method a FittingMethod object to which the DiagOptions options apply
         """
         BaseDiagOptions.__init__(self, method)


      def __str__(self):
         """! 
            String representation of the DiagOptions object 
            @return a Python String describing the DiagOptions object
         """
         return BaseDiagOptions.__str__(self)



   # -----------------------------------------------------------------------------------------------
   class FittingInfo(BaseFittingInfo):

      """! 
          Fitting information available from kmpfit, once fitting is complete. 
          @note refer to the corresponding @c kmpfit documentation
      """

      def __init__(self, method, version, status, message, 
                         niter, nfev, nfree, npegged, dof, 
                         xerror, stderr, covar, chi2_min, rchi2_min, resid):

         """!
            Construct a FittingInfo class

            @param method a FittingMethod object to which the FittingInfo options apply
            @param version kmpfit version number
            @param status status code returned by kmpfit
            @param message message returned by kmpfit, if any
            @param niter actual number of iterations
            @param nfev actual number of function evaluations
            @param nfree number of free parameters
            @param npegged number of pegged parameters
            @param dof number of degrees of freedom
            @param xerror asymptotic error, i.e. 1-sigma parameter uncertainties per parameter 
                   (from covariance matrix)
            @param stderr standard errors assuming reduced Chi^2 = 1
            @param covar final covariance matrix
            @param chi2_min final chi^2
            @param rchi2_min final reduced chi^2
            @param resid final residuals
         """         
         BaseFittingInfo.__init__(self, method)

         self._errno = status                   # result code
         self._version = version                # mkfit version
         self._msg = message                    # message, if any
         self._nb_iter = niter                  # number of iterations
         self._nb_fev  = nfev                   # number of function evaluations
         self._nb_free = nfree                  # number of free parameters
         self._nb_pegged = npegged              # number of pegged parameters
         self._nb_dof = dof                     # number of degrees of freedom
         self._xerror = xerror                  # asymptotic error
         self._stderror = stderr                # Error assuming red.chi^2=1
         self._cov_matrix = covar               # Final parameter covariance matrix
         self._chi2_min = chi2_min              # Minimum Chi2
         self._rchi2_min = rchi2_min            # Reduced Chi^2
         self._residuals = resid                # Final residuals

      def __str__(self):
         """! 
            String representation of the FittingInfo object 
            @return a Python String describing the FittingInfo object
         """
         descr = BaseFittingInfo.__str__(self)
         descr += "- kmpfit version: {0}\n".format(self.version)
         descr += "- errno: {0}\n".format(self.errno)
         descr += "- msg: {0}\n".format(self.msg)
         descr += "- nb iter: {0}\n".format(self.nb_iter)
         descr += "- nb fev: {0}\n".format(self.nb_fev)
         descr += "- nb free: {0}\n".format(self.nb_free)
         descr += "- nb pegged: {0}\n".format(self.nb_pegged)
         descr += "- nb dof: {0}\n".format(self.nb_dof)
         descr += "- xerror: {0}\n".format(self.xerror)
         descr += "- stderror: {0}\n".format(self.stderror)
         descr += "- cov_matrix: {0}\n".format(self.cov_matrix)
         descr += "- chi^2: {0}\n".format(self.chi2_min)
         descr += "- reduced chi^2: {0}\n".format(self.rchi2_min)
         descr += "- residuals: {0}\n".format(self.residuals)
         return descr


      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters          

      @property
      def errno(self):
         """! 
            @return the error code produced by kmpfit 
            @retval ier = 1, 2, 3 or 4 if a solution was found (0: success, < 0: failure)
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._errno

      @property
      def version(self):
         """! 
            @return the kmpfit version number
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._version

      @property
      def msg(self):
         """! 
            @return the message returned by kmpfit, if any 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._msg

      @property
      def nb_iter(self):
         """! 
            @return the number of iterations that were required for the fit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._nb_iter

      @property
      def nb_fev(self):
         """! 
            @return the number of function evaluations that were required for the fit 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._nb_fev

      @property
      def nb_free(self):
         """! 
            @return the number of free parameters 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._nb_free

      @property
      def nb_pegged(self):
         """! 
            @return the number of peggued parameters 
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._nb_pegged

      @property
      def nb_dof(self):
         """! 
            @return the number of degrees of freedom
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._nb_dof

      @property
      def xerror(self):
         """! 
            @return the final 1-sigma parameter uncertainties per parameter (from covariance matrix)
            @note refer to the corresponding @c kmpfit documentation, especially regarding the 
                  validity of the interpretation of these errors
         """
         return self._xerror

      @property
      def stderror(self):
         """! 
            @return the standard errors assuming reduced Chi^2 = 1
            @note refer to the corresponding @c kmpfit documentation, especially regarding the 
                  validity of the interpretation of these errors  
         """
         return self._stderror

      @property
      def cov_matrix(self):
         """! 
            @return the parameter covariance matrix
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._cov_matrix

      @property
      def chi2_min(self):
         """! 
            @return the final chi^2
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._chi2_min

      @property
      def rchi2_min(self):
         """! 
            @return the final reduced chi^2
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._rchi2_min

      @property
      def residuals(self):
         """! 
            @return the final residuals
            @note refer to the corresponding @c kmpfit documentation
         """
         return self._residuals


# --- EOF kmpfit.py 
