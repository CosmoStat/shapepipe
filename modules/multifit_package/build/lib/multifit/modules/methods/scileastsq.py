"""! 
   @package multifit.scileastsq Fitting using scipy.optimize.leastsq (Levenberg-Marquardt)
   @author Marc Gentile
   @file scileastsq.py
   Fitting using scipy.optimize.leastsq (Levenberg-Marquardt)
""" 

# --- Python imports
import scipy.optimize
import numpy

# -- External imports
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingMethod(BaseFittingMethod):

   """! 
      Fitting using the @c scipy.optimize.leastsq implementation of Levenberg-Marquardt 
      @note extends parent BaseFittingMethod class
   """

   def __init__(self, method_config_dir='./config/methods', 
                      fitting_config_dir='./config/fitting',
                      method_config_filename="scileastsq.cfg"):
      """! 
         Construct a FittingMethod object 
         @param method_config_dir configuration directory of the method
         @param fitting_config_dir directory of extra configuration files for the method
         @param method_config_filename configuration filename of the method
      """      
      BaseFittingMethod.__init__(self, name="scileastsq", 
                                       method_config_dir=method_config_dir,
                                       fitting_config_dir=fitting_config_dir,
                                       method_config_filename=method_config_filename)
      
      # --- Public properties  

      # Override parent class options
      self._fitting_options = FittingMethod.FittingOptions(self)  # default fitting options
      self._diag_options = FittingMethod.DiagOptions(self)        # default diagnostic options

#      # --- Private properties 
#      try:
#         self.__module = __import__('scipy.optimize', fromlist=[''])  # load module scipy.optimize
#      except:
#         raise(ModelFitter.FittingError("Unable to load the scipy.optimize.leastsq module."))


   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

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

      # --- Invoke fitting routine
      results = scipy.optimize.leastsq(self.objective_func,
                                       numpy.asarray(model.get_param_values()),
                                       args=tuple(args), 
                                       ftol=self.fitting_options.ftol,  
                                       xtol=self.fitting_options.xtol, 
                                       gtol=self.fitting_options.gtol,
                                       col_deriv=self.fitting_options.col_deriv, 
                                       maxfev=self.fitting_options.maxfev, 
                                       diag=self.fitting_options.diag,
                                       full_output=self.diag_options.full_output)

      
      # --- Fitting info and diagnostic 
      if self.diag_options.full_output:
         # --- Store fitting results 
         self._fitted_params, cov_matrix, infodict, mesg, ier = results
         self._fitting_info = FittingMethod.FittingInfo(self, ier, infodict, mesg, cov_matrix)
      else:
         # --- Store fitting results 
         self._fitted_params, mesg, ier = results                          
         self._fitting_info = BaseFittingInfo(self)

      # --- Fitting status
      self.fitting_status = FittingMethod.FittingStatus(self) 
      self.fitting_status.msg = mesg
      self.fitting_status.warn_msgs = []
      self.fitting_custom_data = None
      self.fitting_status.residuals = infodict["fvec"]
      if ier in (1,2,3):
         self.fitting_status.errno = 0         
      elif ier == 4:
         # Failure
         self.fitting_status.errno = -2
      elif ier < 0:
         # Failure
         self.fitting_status.errno = -1

      return self.fitting_status


   # -----------------------------------------------------------------------------------------------
   def chi2_objective_func(self, params, fitting_func, model, obs_data, sigma_noise_sq, args):
      """!
         Chi2 template objective function 
      """

      #residuals = (obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq
      #print "**** CHI2: residuals:", numpy.sum(residuals)

      return numpy.ravel((obs_data - fitting_func(params, model, *args))**2 / sigma_noise_sq)

   # -----------------------------------------------------------------------------------------------
   def chi2_r_objective_func(self, params, 
                                   fitting_func, model, obs_data, sigma_noise_sq, rlambda, args):
      """!
         Chi2 template objective function 
      """

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

   # ~~~~~~~~~~~~~
   # Inner Classes
   # ~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   class FittingOptions(BaseFittingOptions):
      """! 
         Fitting options applicable to @c scipy.optimize.leastsq()
         @note refer to the corresponding @c scipy documentation 
         @see multifit.multifit.BaseFittingOptions
      """

      def __init__(self, method):
         """!
            Construct a FittingOptions class
            @param method a FittingMethod object to which the FittingOptions options apply
         """
         BaseFittingOptions.__init__(self, method)

         # --- Public properties  
         self._col_deriv = 0        # Jacobian func computes derivatives down the columns 
         self._ftol = 1.49012e-08   # Relative error desired in the sum of squares
         self._xtol = 1.49012e-08   # Relative error desired in the approximate solution
         self._gtol = 0.0           # Desired orthogonality
         self._maxfev = 0           # Maximum number of calls to the function
         self._epsfcn = None        # step length for the forward-difference approx of Jacobian
         self._factor = 100         # parameter determining the initial step bound
         self._diag = None          # N positive entries: serve as scale factors for variables
         self._full_output = 1      # if > 0, include full output information
   
         # --- Set default options from configuration
         self._init_from_config(method.config)

      def __str__(self):
         """! 
            String representation of the FittingOptions object 
            @return a Python String describing the FittingOptions object
         """
         descr = BaseFittingOptions.__str__(self)
         descr += "- col_deriv: {0}\n".format(self.col_deriv)
         descr += "- ftol: {0}\n".format(self.ftol)
         descr += "- xtol: {0}\n".format(self.xtol)
         descr += "- gtol: {0}\n".format(self.gtol)
         descr += "- maxfev: {0}\n".format(self.maxfev)
         if self.epsfcn is not None:
            descr += "- epsfcn: {0}\n".format(self.epsfcn)
         descr += "- factor: {0}\n".format(self.factor)
         if self.diag is not None:
            descr += "- diag: {0}\n".format(self.diag)
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters 

#      @property
#      def method(self):
#         """! 
#            @return the FittingMethod object to which the FittingOptions apply 
#            @note refer to the corresponding @c scipy documentation
#         """
#         return self._method

      @property
      def col_deriv(self):
         """! 
            @return the @c col_deriv option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._col_deriv

      @property
      def ftol(self):
         """! 
            @return the @c ftol option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._ftol

      @property
      def gtol(self):
         """! 
            @return the @c gtol option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._gtol

      @property
      def xtol(self):
         """! 
            @return the @c xtol option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._xtol

      @property
      def maxfev(self):
         """! 
            @return the @c maxfev option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._maxfev

      @property
      def epsfcn(self):
         """! 
            @return the @c epsfcn option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._epsfcn

      @property
      def factor(self):
         """! 
            @return the @c factor option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._factor

      @property
      def diag(self):
         """! 
            @return the @c diag option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._diag


      # --- Setters 

      @col_deriv.setter 
      def col_deriv(self, col_deriv):
         """! 
            Set the @c col_deriv option value of @c scipy.optimize.leastsq()
            @param col_deriv the @c col_deriv option value
            @note refer to the corresponding @c scipy documentation
         """
         self._col_deriv = col_deriv

      @ftol.setter 
      def ftol(self, ftol):
         """! 
            Set the @c ftol option value of @c scipy.optimize.leastsq()
            @param ftol the @c gtol option value
            @note refer to the corresponding @c scipy documentation
         """
         self._ftol = ftol

      @gtol.setter 
      def gtol(self, gtol):
         """! 
            Set the @c gtol option value of @c scipy.optimize.leastsq()
            @param gtol the @c gtol option value
            @note refer to the corresponding @c scipy documentation
         """
         self._gtol = gtol

      @xtol.setter 
      def xtol(self, xtol):
         """! 
            Set the @c xtol option value of @c scipy.optimize.leastsq()
            @param xtol the @c xtol option value
            @note refer to the corresponding @c scipy documentation
         """
         self._xtol = xtol

      @maxfev.setter 
      def maxfev(self, maxfev):
         """! 
            Set the @c maxfev option value of @c scipy.optimize.leastsq()
            @param maxfev the @c maxfev option value
            @note refer to the corresponding @c scipy documentation
         """
         self._maxfev = maxfev

      @epsfcn.setter 
      def epsfcn(self, epsfcn):
         """! 
            Set the @c epsfcn option value of @c scipy.optimize.leastsq()
            @param epsfcn the @c epsfcn option value
            @note refer to the corresponding @c scipy documentation
         """
         self._epsfcn = epsfcn

      @factor.setter 
      def factor(self, factor):
         """! 
            Set the @c factor option value of @c scipy.optimize.leastsq()
            @param factor the @c factor option value
            @note refer to the corresponding @c scipy documentation
         """
         self._factor = factor

      @diag.setter 
      def diag(self, diag):
         """! 
            Set the @c diag option value of @c scipy.optimize.leastsq()
            @param diag the @c diag option value
            @note refer to the corresponding @c scipy documentation
         """
         self._diag = diag

      # ~~~~~~~~~~~~~~~~
      # Private methods
      # ~~~~~~~~~~~~~~~~

      def _init_from_config(self, config):

         if not config is None:
            method_section_name = self.method.name.upper()
            model_options_section_name = "{0}.OPTIONS".format(method_section_name)
            if config.has_section(method_section_name):

               # --- Default configuration
               if config.has_key("col_deriv", model_options_section_name):
                  self.col_deriv = config.get_as_int("col_deriv", model_options_section_name)
               if config.has_key("ftol", model_options_section_name):
                  self.ftol = config.get_as_float("ftol", model_options_section_name)
               if config.has_key("gtol", model_options_section_name):
                  self.gtol = config.get_as_float("gtol", model_options_section_name)
               if config.has_key("xtol", model_options_section_name):
                 self.xtol = config.get_as_float("xtol", model_options_section_name)
               if config.has_key("max_fev", model_options_section_name):
                  self.maxfev = config.get_as_int("max_fev", model_options_section_name)
               if config.has_key("factor", model_options_section_name):
                  self.factor = config.get_as_float("factor", model_options_section_name)
               if config.has_key("epsfcn", model_options_section_name): # may be None
                 self.epsfcn = eval(config.get_as_string("epsfcn", model_options_section_name))  
               if config.has_key("diag", model_options_section_name): # may be None
                  self.diag = eval(config.get_as_string("diag", model_options_section_name))      
         else:
            if not self.method.name is None:
               # --- No configuration file found => taking suitable defaults             
               print("MultiFit *** Warning ***: configuration file not found for {0}.".format(
                                                                                 self.method.name))   
 
            

   # -----------------------------------------------------------------------------------------------
   class DiagOptions(BaseDiagOptions):
      """! 
         Diagnostic options applicable to scipy.optimize.leastsq()
         @note refer to the corresponding @c scipy documentation 
         @see multifit.multifit.BaseDiagOptions
      """

      def __init__(self, method):
         """!
            Construct a DiagOptions class
            @param method a FittingMethod object to which the DiagOptions options apply
         """
         BaseDiagOptions.__init__(self, method)

         # --- Public properties  
         self._full_output = 0    # non-zero to return all optional outputs
         
         # --- Initialize options from configuration
         self._init_from_config(method.config)

      def __str__(self):
         """! 
            String representation of the DiagOptions object 
            @return a Python String describing the DiagOptions object
         """
         descr = BaseDiagOptions.__str__(self)
         descr += "- full output: {0}\n".format(self.full_output)
         return descr

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~   

      # --- Getters 

      @property
      def method(self):
         """! 
            @return the FittingMethod object to which the DiagOptions apply 
            @note refer to the corresponding @c scipy documentation
         """
         return self._method

      @property
      def full_output(self):
         """! 
            @return the @c full_output option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._full_output

      # --- Setters 

      @full_output.setter 
      def full_output(self, full_output):
         """! 
            Set the @c full_output option value of @c scipy.optimize.leastsq()
            @param full_output the @c full_output option value
            @note refer to the corresponding @c scipy documentation
         """
         self._full_output = full_output

      # ~~~~~~~~~~~~~~~~
      # Private methods
      # ~~~~~~~~~~~~~~~~

      def _init_from_config(self, config):

         # --- Default configuration
         if not config is None:
            method_section_name = self.method.name.upper()
            model_options_section_name = "{0}.OPTIONS".format(method_section_name)
            self.full_output = config.get_as_int("full_output", model_options_section_name)
         else:
            if not self.method.name is None:
               print("MultiFit *** Warning ***: configuration file notfound for {0}.".format(
                                                                                 self.method.name))   



   # -----------------------------------------------------------------------------------------------
   class FittingInfo(BaseFittingInfo):

      """! 
          Fitting information available from scipy.optimize.leastsq(), once fitting is complete. 
          @note refer to the corresponding @c scipy documentation
      """

      def __init__(self, method, ier, infodict, mesg, cov_x):
         """!
            Construct a FittingInfo class

            @param method a FittingMethod object to which the FittingInfo options apply
            @param ier error code returned by scipy.optimize.leastsq()
            @param infodict extra optional information returned by scipy.optimize.leastsq()
            @param mesg message giving information about the cause of failure
            @param cov_x estimate of the jacobian around the solution
         """         
         BaseFittingInfo.__init__(self, method)

         self._ier = ier            # if 1, 2, 3 or 4, a solution was found
         self._infodict = infodict  # dictionary of optional outputs
         self._mesg = mesg          # string message giving information about the cause of failure
         self._cov_x = cov_x        # estimate of the jacobian around the solution            

      def __str__(self):
         """! 
            String representation of the FittingInfo object 
            @return a Python String describing the FittingInfo object
         """
         descr = BaseFittingInfo.__str__(self)
         descr += "- ier: {0}\n".format(self.ier)
         #descr += "- infodict: {0}\n".format(self.infodict)
         descr += "- infodict:\n"
         for key in self.infodict:
            descr += "  - key {0}: {1}\n".format(key, self.infodict[key])          
  
         descr += "- message (mesg): {0}\n".format(self.mesg)
         descr += "- Jacobian estimate (cov_x): {0}\n".format(self.cov_x)
         descr += "- residuals: {0}\n".format(self.residuals)

         return descr

      # ~~~~~~~~~~~
       # Properties
      # ~~~~~~~~~~~   

      # --- Getters          

      @property
      def ier(self):
         """! 
            @return the @c ier option value of @c scipy.optimize.leastsq() 
            @retval ier = 1, 2, 3 or 4 if a solution was found
            @note refer to the corresponding @c scipy documentation
         """
         return self._ier

      @property
      def mesg(self):
         """! 
            @return the @c mesg option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._mesg

      @property
      def infodict(self):
         """! 
            @return the @c infodict option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._infodict

      @property
      def cov_x(self):
         """! 
            @return the @c cov_x option value of @c scipy.optimize.leastsq() 
            @note refer to the corresponding @c scipy documentation
         """
         return self._cov_x

      @property
      def residuals(self):
         """! 
            @return the final residuals
            @note refer to the corresponding @c scipy documentation
         """
         return self.infodict["fvec"]

# --- EOF scileastsq.py 
