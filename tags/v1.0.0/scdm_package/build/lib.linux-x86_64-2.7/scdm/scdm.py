"""! 
   @package scdm.scdm Simple Coordinate Descent Minimizer
   @author Marc Gentile
   @file scdm.py
   A simple coordinate descent minimizer
   A Simple Coordinate Descent Minimizer
"""

import os, sys
import math
import time
import numpy as np
if __debug__:
   import slogger
   import pylab
   import matplotlib
   import matplotlib.figure
   import mpl_toolkits.mplot3d.axes3d as pylab3
   from matplotlib.font_manager import FontProperties
   from matplotlib.pylab import setp

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter of an objective function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Param(object):

   """! Represents a parameter of some model which has to be estimated """
   def __init__(self, name, guess_value=0, param_range=[None, None], 
                step_size=0.01, step_range=[None, None], step_factor=2,
                scaling_func=None, max_tol=1e-8, 
                out_of_bound_handler=None, numerical_error_handler=None, is_constant=False):
      """!
         Construct an object of class Param
         @param name name of the parameter
         @param guess_value initial guess value 
         @param param_range allowed range of variation of the parameter value as a 
                list [@c min, @c max]
         @param step_size initial step size for the variation of the parameter
         @param step_range allowed range of variation of Param.step_size as a 
                list [@c min, @c max]
         @param step_factor initial step factor for increasing or decreasing the parameter value
         @param scaling_func user-provided scaling function (e.g. power, log...)
         @param max_tol maximum tolerance allowed in the variation of the parameter value
         @param out_of_bound_handler callback function to notify a parameter value came outside the 
                allowed range of varation specified by Param.guess_value
         @param numerical_error_handler callback function to notify a numerical error that occurred
                while computing the parameter value
         @param is_constant if true, the parameter value will not vary (parameters not fitted)
      """

      self._name = name
      self._extrema_step_range  = [1.0e-16, 1.0e+16]
      self._extrema_param_range = [-1.0e16, 1.0e+16]
      self._param_range = param_range
      self._guess_value = guess_value
      self._step_range  = step_range
      self._step_size   = step_size
      self._step_factor = step_factor
      self._scaling_func = scaling_func
      self._is_constant = is_constant
      self._max_tol = max_tol               
      if out_of_bound_handler is None:
         self._out_of_bound_handler = self._default_out_of_bound_handler
      else:
         self._out_of_bound_handler = out_of_bound_handler

      if numerical_error_handler is None:
         self._numerical_error_handler = self._default_numerical_error_handler
      else:
         self._numerical_error_handler = numerical_error_handler

      self._check_input()   # validate input parameters


   def __str__(self):
      """! 
         String representation of the Param object 
         @return a Python String describing the Param object
      """
      return "<{0}, [{1:.2f}, {2:.2f}], {3:.2f}, {4}>".format(
                                     self._name, self._param_range[0], self._param_range[1], 
                                     self._guess_value, self._is_constant)

   # ~~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~~
   
   # --- Getters ---

   @property
   def name(self):
      """! 
         Return the parameter name
         @return the parameter name 
      """
      return self._name

   @property
   def param_range(self):
      """! 
         Return the allowed variation range of parameter values
         @return the allowed variation range 
      """
      return self._param_range

   @property
   def guess_value(self):
      """! 
         Return the initial parameter guess value
         @return the initial parameter guess value 
      """
      return self._guess_value

   @property
   def step_size(self):
      """!
         Return the step size for varying the parameter
         @return the step size for varying the parameter 
      """ 
      return self._step_size

   @property
   def step_range(self):
      """! Return the allowed variation range of the parameter step size 
         @return the allowed variation range of the parameter step size 
      """
      return self._step_range

   @property
   def step_factor(self):
      """! 
         Return the multiplicative factor applied when varying the parameter value 
         @return the multiplicative factor applied when varying the parameter value 
      """      
      return self._step_factor

   @property
   def scaling_func(self):
      """! 
         Return the user-provided scaling function (e.g. power, log...) for varying the parameter
                 values
         @return the user-provided scaling function
      """   
      return self._scaling_func   

   @property
   def max_tol(self):
      """! 
         Return the maximum tolerance allowed in the variation of the parameter value
         @return the maximum tolerance allowed 
      """
      return self._max_tol

   @property
   def out_of_bound_handler(self):
      """! 
         Return the callback handler function to notify a parameter value came outside the 
                 allowed range of varation specified by Param.guess_value
         @return the callback handler function
      """
      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      return self._numerical_error_handler

   @property
   def is_constant(self):
      """! 
         Tell wether the parameter value should remain constant during fitting, that is, if the 
         parameter must be fitted at all.
         @return True if the parameter value remain constant during fitting, 
                 False otherwise (default)
      """
      return self._is_constant

   # --- Setters ---

   @name.setter
   def name(self, name):
      """! 
         Set the name of the parameter
         @param name the parameter name
      """
      self._name = name 

   @param_range.setter
   def param_range(self, param_range):
      """! 
         Set the allowed variation range of parameter values
         @param param_range the allowed variation range
      """
      print "*** SCDM: setting param range to:", self._name, param_range
      self._param_range = param_range 

   @guess_value.setter
   def guess_value(self, guess_value):
      """! 
         Set the initial parameter guess value
         @param guess_value the initial parameter guess value
      """
      self._guess_value = guess_value

   @step_size.setter
   def step_size(self, step_size):
      """! 
         Set the step size for varying the parameter        
         @param step_size the step size
      """
      self._step_size = step_size

   @step_range.setter
   def step_range(self, step_range):
      """! 
         Set the allowed variation range of the parameter step size       
         @param step_range allowed step range
      """
      self._step_range = step_range

   @step_factor.setter
   def step_factor(self, step_factor):
      """! 
         Set the multiplicative factor applied when varying the parameter value      
         @param step_factor the step factor
         @note the step factor should be greater than unity
      """
      self._step_factor = step_factor

   @scaling_func.setter
   def scaling_func(self, scaling_func):
      """! 
         Set the user-provided scaling function (e.g. power, log...) for varying the parameter
              values       
         @param scaling_func the scaling function
      """
      self._scaling_func = scaling_func

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

   @is_constant.setter
   def is_constant(self, is_constant):
      """! 
         Set the parameter as constant (not varied)
         @param is_constant True if parameter set as constant (not varied), False otherwise 
      """
      self._is_constant = is_constant


   # ~~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~~

   def _check_input(self):
      """! Checks input parameter values """
      if self._step_factor <= 1:
         print('Warning: step factor {0:3.1f} for {1} must be greater than 1 "\
               "=> reset to 2.0'.format(self._step_factor, self._name))
         self._step_factor = 2.0

      # --- Step Range
      if self._step_range[0] is None:
         self._step_range[0] = self._extrema_step_range[0]
      if self._step_range[1] is None:
         self._step_range[1] = self._extrema_step_range[1]

      # --- Param range
      if self._param_range[0] is None:
         self._param_range[0] = self._extrema_param_range[0] 
      if self._param_range[1] is None:
         self._param_range[1] = self._extrema_param_range[1]  

   def _default_out_of_bound_handler(self, iiter, param, param_value, bound_index):
      """! Default parameter out of bound handler """
      # Default: reset parameter to original bound and continue fitting
      return [0, param_value, param.param_range[bound_index]]                             

   def _default_numerical_error_handler(self, iiter, param, value, value_name):
      """! Default parameter numerical error handler """
      return -1    # default: abort fitting in both Nan and Inf cases

# ~~~~~~~~~~~~~~~~~~~~~~~~~
# Hold all input parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~
class Params(object):

   """! 
      Represents the list of parameters of clas Param to estimate 
   """
   def __init__(self, param_list=[]):
      """!
         Construct an object of class Params
         @param param_list a list of Param objects 
         @see Param
      """
      self._param_list = param_list
      self._param_names = [p.name for p in param_list]

   def __str__(self):
      """! 
         String representation of the Params object 
         @return a Python String describing the Params object
      """
      return "[{0}]".format(", ".join(["{0}".format(p) for p in self._param_list]))


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   @property      
   def param_list(self):   
      """!
         Return the list of parameters objects
         @return the list of parameters as Param objects
         @see Param
      """
      return self._param_list

   @property      
   def param_names(self):   
      """!
         Return the list of parameters object names
         @return the parameter names of Params 
         @see Param
      """
      return self._param_names

   @param_list.setter
   def param_list(self, param_list):
      """!
         Set the list of parameters objects
         @param param_list a Params object 
         @see Param
      """
      self._param_list = param_list


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define what should be mimimized during fitting 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class MinimizationType(object):
   """!
      Define which quantity should be mimimized during fitting and hence the stopping criterion
      - @c minimize_func_value: minimize the mean variations of the objective function value over 
                                a window of the last @c N function value estimates; the estimate is
                                made against FittingOptions.max_func_tol 
      - @c minimize_param_delta: minimize the mean variation of each parameter over a window of 
                                the last @c N parameter value estimates; the estimate is
                                made against Param.max_tol 
      - @c minimize_param_sigma: minimize the standard deviation of the variations of each 
                                 parameter over a window of the last @c N parameter value estimates;
                                 the estimate is made against Param.max_tol
      @note one can specify a combination of the three aforementioned constants, like: 
            @c minimize_func_value | @c minimize_param_delta | @c minimize_param_sigma 
   """      
   (minimize_func_value, minimize_param_delta, minimize_param_sigma) = (1,2,4)


# ~~~~~~~~~~~~~~~
# Diagnose level 
# ~~~~~~~~~~~~~~~
class DiagFlag(object):
   """! 
      Define which diagnostic information should be collected
      - @c logs: produce diagnostic log files 
      - @c plots: produce diagnostic plots
      - @c plots: produce diagnostic statistics 
      @note one can specify a combination of the three aforementioned constants, like: 
            @c logs | @c plots | @c stats
   """   
   (logs, plots, stats) = (1,2,4)

# ~~~~~~~~~~~~~~~~
# Fitting options
# ~~~~~~~~~~~~~~~~
class FittingOptions(object):

   """!
      Class for specifying options applicable to the fitting process   
   """
   def __init__(self, max_iter=1000, max_fev=1000, max_func_tol=1e-8, all_minima=False,
                      mini_type=MinimizationType.minimize_func_value |
                      MinimizationType.minimize_param_delta |
                      MinimizationType.minimize_param_sigma | 
                      MinimizationType.minimize_func_value):

      # --- Public Options
      self._max_iter = max_iter              # maximum number of iterations allowed: beyond, fitting will stop and paraameter estimates may not be accurate
      self._max_fev  = max_fev               # maximum number of function evaluations (fev). There might me as many fev per iteration as there are parameters, or less
      self._max_func_tol = max_func_tol      # maximum allowed error on the residuals between model and data 
      self._all_minima  = all_minima         # if True, attempts to locate all minima (TODO) 
      self._mini_type = mini_type            # what to minimize: minimise global residuals (all parameters), error on individual parameters or both 

      # --- Internal options
      self._dcount_step_inc_threshold = 1    # number of steps in the same direction beyond which the step size should start increasing by a factor <param.step_factor> 
      self._stop_delta_func_win_size  = 10   # Size of stack of global delta_func values to check before stopping the fit         
      self._stop_delta_param_win_size = 10   # Size of stack of individual param delta_func values to check before stopping the fit
      self._stop_sigma_param_win_size = 10   # Size of stack of individual param parameter values on with to calculate the standard deviation and possibly stop the fit
      self._stop_step_param_win_size  = 10   # Size of stack of individual param parameter values on with to evaluate the last step values and possibly stop the fit
      self._data_window_size = max(self._stop_delta_func_win_size, self._stop_delta_param_win_size, self._stop_sigma_param_win_size)  

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters ---

   @property
   def max_iter(self):
      return self._max_iter

   @property
   def max_fev(self):
      return self._max_fev

   @property
   def max_func_tol(self):
      """! 
         Return the maximum tolerance allowed in the variation of the objective function value
         @return maximum tolerance in the variation of the objective function value
      """
      return self._max_func_tol

   @property
   def all_minima(self):
      return self._all_minima

   @property
   def mini_type(self):
      return self._mini_type

   @property
   def data_window_size(self):
      return self._data_window_size

   # --- Setters ---

   @max_iter.setter
   def max_iter(self, max_iter):
      self._max_iter = max_iter

   @max_fev.setter
   def max_fev(self, max_fev):
      self._max_fev = max_fev

   @max_func_tol.setter
   def max_func_tol(self, max_func_tol):
      self._max_func_tol = max_func_tol

   @all_minima.setter
   def all_minima(self, all_minima):
      self._all_minima = all_minima

   @mini_type.setter
   def mini_type(self, mini_type):
      self._mini_type = mini_type

   @data_window_size.setter
   def data_window_size(self, data_window_size):
      self._data_window_size = data_window_size


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diagnostic options
# [to output logs, plots and statistics to tune fitting or diagnose fitting failure]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class DiagOptions(object):
   """! 
     Diagnostic helper class to output logs, plots and statistics to 
     tune fitting or diagnose fitting failure 
   """

   def __init__(self, diag_flags=DiagFlag.logs | DiagFlag.plots | DiagFlag.stats, base_output_directory='.'):

      self._diagnose  = False            
      self._diag_flags = diag_flags
      if  diag_flags is None:
         diag_flags = 0  
      self._base_output_dir  = base_output_directory
      self._log_output_dir   = os.path.join(self.base_output_dir, 'logs')
      self._plot_output_dir  = os.path.join(self.base_output_dir, 'plots')
      self._stats_output_dir = os.path.join(self.base_output_dir, 'stats')
      self._log_file_name = 'scdm.log'
      

   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters ---

   @property
   def diag_flags(self):
      return self._diag_flags

   @property
   def base_output_dir(self):
      return self._base_output_dir

   @property
   def log_output_dir(self):
      return self._log_output_dir

   @property
   def plot_output_dir(self):
      return self._plot_output_dir

   @property
   def stats_output_dir(self):
      return self._stats_output_dir

   @property
   def log_file_name(self):
      return self._log_file_name


   # --- Setters ---

   @diag_flags.setter
   def diag_flags(self, diag_flags):
      self._diag_flags = diag_flags

   @base_output_dir.setter
   def base_output_dir(self, base_output_dir):
      self._base_output_dir = base_output_dir

   @log_output_dir.setter
   def log_output_dir(self, log_output_dir):
      self._log_output_dir = log_output_dir

   @plot_output_dir.setter
   def plot_output_dir(self, plot_output_dir):
      self._plot_output_dir = plot_output_dir

   @stats_output_dir.setter
   def stats_output_dir(self, stats_output_dir):
      self._stats_output_dir = stats_output_dir

   @log_file_name.setter
   def log_file_name(self, log_file_name):
      self._log_file_name = log_file_name


   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~

   def is_diagnose(self):
      return self.logging_enabled() and self.plotting_enabled() and self.stats_enabled()

   def logging_enabled(self):
      return self._diag_flags & DiagFlag.logs == DiagFlag.logs

   def plotting_enabled(self):
      return self._diag_flags & DiagFlag.plots == DiagFlag.plots

   def stats_enabled(self):
      return self._diag_flags & DiagFlag.stats == DiagFlag.stats


# ~~~~~~~~~~~~~~~
# Minmizer class
# ~~~~~~~~~~~~~~~
class Minimizer(object):
   """! 
       Minimization of an objective function using a variant of the Coordinate Descent algorithm 
   """

   def __init__(self, params=Params(), target_params=Params(), func=None, args=[], fitting_options=FittingOptions(), diag_options=DiagOptions()):

      # Options
      self._fitting_options = fitting_options      # fitting options
      self._diag_options = diag_options            # diagnostic/tuning options

      # Fitting
      self._params = params                # all parameters to fit

      self._target_params = target_params  # stop fitting as soon as parameters in list <targetParams> are fitted 
      self._args = args                    # arguments to pass to evaluation function
      if func is None:
         self._func = self._dummy_func(self._params)
      else:
         self._func = func

      self._setup_options()

      self._fitting_cache = self._FittingCache()    # hold cached data with parameter tuple (p0, p1, ... pn) as key
      self._fitting_results = None                  # hold fitting results for all parameters 
      self._fitting_status = self._FittingStatus()  # result status of fitting

      # Helpers
      self._logger = None
      self._plot_helper  = None
      self._stats_helper = None


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~

   # --- Getters ---
   
   @property
   def params(self):
      return self._params

   @property
   def target_params(self):
      return self._target_params
   
   @property
   def args(self):
      return self._args

   @property
   def func(self):
      return self._func

   @property
   def fitting_options(self):
      return self._fitting_options

   @property
   def diag_options(self):
      return self._diag_options

   @property
   def fitting_data(self):
      return self._fitting_data

   @property
   def fitting_results(self):
      """! Returns available fitting results """
      if not self._fitting_results is None:
         del self._fitting_results
      self._fitting_results = Minimizer._FittingResult(self._fitting_data)

      return self._fitting_results

   @property
   def fitting_cache(self):
      return self._fitting_cache

   @property
   def fitting_status(self):
      return self._fitting_status

   @property
   def base_output_dir(self):
      return self._base_output_dir

   @property
   def logger(self):
      return self._logger

   @property
   def plot_helper(self):
      return self._plot_helper

   @property
   def stats_helper(self):
      return self._stats_helper


   # --- Setters ---

   @params.setter
   def params(self, params):
      self._params = params
      # Keep an history of size <data_window_size> 

      self._fitting_data  = self._FittingData(self._params, 
                                              self.fitting_options._data_window_size) 

   @target_params.setter
   def target_params(self, target_params):
      self._target_params = target_params

   @args.setter
   def args(self, args):
      self._args = args

   @func.setter
   def func(self, func):
      self._func = func

   @fitting_options.setter
   def fitting_options(self, fitting_options):
      self._fitting_options = fitting_options
      self._setup_options()

   @diag_options.setter
   def diag_options(self, diag_options):
      self._diag_options = diag_options
      self._setup_options()

   @logger.setter
   def logger(self, logger):
      self._logger = logger

#   @fitting_data.setter
#   def fitting_data(self, fitting_data):
#      self._fitting_data = fitting_data

#   @fitting_status.setter
#   def fitting_status(self, fitting_status):
#      self._fitting_status = fitting_status

#   @fitting_results.setter
#   def fitting_results(self, fitting_results):
#      self._fitting_results = fitting_results

#   @fitting_cache.setter
#   def fitting_cache(self, fitting_cache):
#      self._fitting_cache = fitting_cache

#   @base_output_dir.setter
#   def base_output_dir(self, base_output_dir):
#      self._base_output_dir = base_output_dir

   @logger.setter
   def logger(self, logger):
      self._logger = logger

#   @plot_helper.setter
#   def plot_helper(self, plot_helper):
#      self._plot_helper = plot_helper

#   @stats_helper.setter
#   def stats_helper(self, stats_helper):
#      self._stats_helper = stats_helper


   # ~~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~~

   def is_diagnose(self):
      return self._diag_options.is_diagnose()

   def is_logging_enabled(self):
      return self._diag_options.logging_enabled()

   def is_plotting_enabled(self):
      return self._diag_options.plotting_enabled()

   def is_stats_enabled(self):
      return self._diag_options.stats_enabled()

   def get_log_directory(self):
      return self._diag_options.log_output_dir

   def get_log_filename(self):
      return self._diag_options.log_file_name

   def get_plot_directory(self):
      return self._diag_options.plot_output_dir
   
   def get_stats_directory(self):
      return self._diag_options.stats_output_dir

   def get_nb_iter(self):
      return self._fitting_results.nb_iter

   def get_nb_func_eval(self):
      return self._fitting_results.get_nb_func_evals()

   def get_fitted_param_values(self):
      return self._fitting_results.fitted_param_values

   # ~~~~~~~~~~~~~~~~~~~~~~~~
   # Fit a set of parameters
   # ~~~~~~~~~~~~~~~~~~~~~~~~
   def fit_parameters(self):

      # Stopping condition
      stop_fitting = False

      # Total number of fitting operations
      tot_nb_oper = 0

      # Total number of function evaluations
      tot_nb_fev = 1

      # Stuff to do at iteration 0
      prev_param_values = [p.guess_value for p in self.params.param_list]  # first parameter set initialised with input guess values
      nb_params = len(prev_param_values)

      try:

         # Initial call to custom fitting user function
         result = self.func(prev_param_values, *self.args)
         if (type(result) == type([])):
            if len(result) == 2: # user passed custom data as 2nd list element  
               [last_func_value, custom_object] = result
            elif len(result) > 2:
               print('Error: evaluation function may return at most 2 parameters')
               self.fitting_status.set_fitting_status(-10, 'Error - Incorrect number of parameter returned in evaluation function') 
               return [prev_param_values, self.fitting_results]            
            elif len(result) == 1:
               [last_func_value, custom_object] = [result[0], None]
               print('Error: custom function not found as second returned parameter of evaluation function')
               self.fitting_status.set_fitting_status(-10, 'Error - Incorrect number of parameter returned in evaluation function') 
               return [prev_param_values, self.fitting_results]            
         else:
            [last_func_value, custom_object] = [result, None]  # assumed to be a number

         has_custom_data = (not custom_object is None)

         # Target param list: actually we may only want to fit parameters in that list and not all specified parameters
         #                    => we stop fitting as soon values for target parameters are within TOL
         # Check that the list is not empty. If so, replace it by the full list of parameters
         if len(self.target_params.param_names) == 0:
            self.target_params = self.params

         # Global fitting data (all parameters)
         fitting_data = self._fitting_data
         fitting_data.clear_data()
         self._fitting_results = None

         # Cached fitting data (key is parameter tuple (p0, p1, ..., pn))
         cached_data = self.fitting_cache
         cached_data.clear()

         # Clear status information
         self.fitting_status.clear()

         # Create an initial record of fitting execution
         cached_data.add_record(prev_param_values, last_func_value)

         # Start fitting...
         for iiter in xrange(1, self.fitting_options.max_iter+1):
            
            if __debug__:
               if self.is_logging_enabled(): 
                  #self.logger.log_blank_line()
                  self.logger.log_info(
                      '=' * 76 + ' Fitting iteration: {0} '.format(str(iiter) + ' ' + '=' * 76) )
      
            # Vary each parameter...      
            for param in self.params.param_list:

               rank = 1 + self.params.param_list.index(param)

               # Fitting data for target parameter
               param_fitting_data = fitting_data.get_param_fitting_data_for(param.name)

               # Total number of operations (for all parameters)
               tot_nb_oper += 1

               # Check if we still must vary this parameter
               if param_fitting_data.stop_varying:
                  if __debug__:
                     if self.is_logging_enabled(): 
                        self.logger.log_info('*** Parameter: {0} not varied anymore...'.format(param.name)) 
                        self.logger.log_info('-' * 177)

                  # But we still keep a record of the data...

                  # Record global fitting data for this iteration and parameter
                  if has_custom_data:
                     fitting_data.record_fitting_data(index=tot_nb_oper, feval=fitting_data.get_last_fitting_data('feval'),
                                                                    delta_func=fitting_data.get_last_fitting_data('delta_func'),
                                                                    custom_obj=fitting_data.get_last_fitting_data('custom_objs'))
                  else:
                     fitting_data.record_fitting_data(index=tot_nb_oper, feval=fitting_data.get_last_fitting_data('feval'), 
                                                                    delta_func=fitting_data.get_last_fitting_data('delta_func'),
                                                                    custom_obj=custom_object)

                  # Record parameter fitting data for this iteration and parameter
                  param_fitting_data.record_fitting_data(pval=param_fitting_data.get_last_fitting_data('params'), feval=param_fitting_data.get_last_fitting_data('feval'), 
                                                     deval=param_fitting_data.get_last_fitting_data('deval'), step=param_fitting_data.get_last_fitting_data('steps'),
                                                     delta_func=param_fitting_data.get_last_fitting_data('delta_func'), dcount=param_fitting_data.get_last_fitting_data('dcount'), iiter=iiter)

                  continue

               # Varying parameter ...
               if __debug__:
                  if self.is_logging_enabled(): 
                     self.logger.log_info('Varying Parameter: {0}...'.format(param.name)) 

               param_index = rank-1  # parameter zero-based index                       
               if not param_fitting_data.has_records():
                  step = param.step_size
               else:
                  # Last step for target parameter
                  step = param_fitting_data.get_last_fitting_data('steps')

               # Stop varying if the parameter has been declared as constant
               if param.is_constant:
                  param_fitting_data.stop_varying = True # will stop varying at the next iteration
                  new_param_value = prev_param_values[param_index]
               else:
                  # TODO: to redesign: pass a user-defined parameter scaling function
                  # Parameter scaling: Location or Scale
                  new_param_value = prev_param_values[param_index] + step

#               isScaleParam = (param.scaling == ParamScaling.Scale)

#               # Vary target parameter value
#               if not isScaleParam:
#                  newScaledParamValue = prev_param_values[param_index] + step
#                  new_param_value = newScaledParamValue
#               else:
#                  newScaledParamValue = math.log10(prev_param_values[param_index]) + step
#                  new_param_value = 10**newScaledParamValue
#                  #self.logger.log_info('Scaled: %11.6e  Unscaled: %11.6e' %(newScaledParamValue, new_param_value))

               # Check if the move if compatible with specified parameter range
               #if (not param.param_range[0] is None) and new_param_value < param.param_range[0]:

               if new_param_value < param.param_range[0]:
                  if __debug__:
                     if self.is_logging_enabled(): 
                        self.logger.log_warning('Lower bound of parameter was reached: {0}={1:6.3e} < {2:6.3e}. Value reset to {3:6.3e}'.format(
                                                 param.name, new_param_value, param.param_range[0], param.param_range[0]))

                  [status, old_param_value, new_param_value] = param.out_of_bound_handler(iiter, param, new_param_value, 0)
                  if status < 0:
                     # Give up and abort fitting
                     self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Lower bound {2:6.3e} of parameter was reached {3:6.3e} => aborting fitting'.format(
                                                       iiter, param.name, param.param_range[0], new_param_value))
                     raise Minimizer.ParamOutOfBoundException(iiter, param, old_param_value, 0, '')
                  else:
                     # Attempt to fix parameter
                     #new_param_value = param.param_range[0] 
                     param_fitting_data.lowerBoundaryReached = True
                     self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Lower bound {2:6.3e} of parameter was reached {3:6.3e}'.format(
                                                          iiter, param.name, param.param_range[0], new_param_value))

               #if (not param.param_range[1] is None) and new_param_value > param.param_range[1]:
               if new_param_value > param.param_range[1]:
                  if __debug__:
                     if self.is_logging_enabled(): 
                        self.logger.log_warning('Upper bound of parameter was reached: {0}={1:6.3e} > {2:6.3e}. Value reset to {3:6.3e}'.format(
                                                 param.name, new_param_value, param.param_range[1], param.param_range[1]))
                  [status, old_param_value, new_param_value]  = param.out_of_bound_handler(iiter, param, new_param_value, 1)
                  if status < 0:
                     # Give up and abort fitting
                     self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Upper bound {2:6.3e} of parameter was reached {3:6.3e} => aborting fitting'.format(
                                                          iiter, param.name, param.param_range[1], new_param_value))
                     raise Minimizer.ParamOutOfBoundException(iiter, param, old_param_value, 1, '')
                  else:
                     # Attempt to fix parameter
                     #new_param_value = param.param_range[1]
                     param_fitting_data.upperBoundaryReached = True
                     self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Upper bound {2:6.3e} of parameter was reached {3:6.3e}'.format(
                                                          iiter, param.name, param.param_range[1], new_param_value))

               # Actual nb of iterations (some parameter variations may be skipped for efficiency)
               aiter = param_fitting_data.inc_iter()
    
               # TODO: optimize this
               new_param_values = prev_param_values[:param_index] + [new_param_value]
               new_param_values.extend(prev_param_values[param_index+1:]) 
               if __debug__:
                  if self.is_logging_enabled(): 
                     self.logger.log_info('Last parameter set: {0}'.format(prev_param_values))
                     self.logger.log_info('New parameter set:  {0}'.format(new_param_values))

               # If no fitting data for target parameter, create an initial one
               if not param_fitting_data.has_records():
                  param_fitting_data.record_fitting_data(pval=param.guess_value, feval=last_func_value, 
                                                         deval=math.fabs(last_func_value), step=step, 
                                                         delta_func=last_func_value, dcount=0, iiter=aiter)
                  if __debug__:
                     if self.is_logging_enabled(): 
                        self.logger.log_info('First recorded data - Params: {0} Step: {1:6.3e} delta_func: {2}'.format(
                                              param_fitting_data.get_last_fitting_data('params'), 
                                              param_fitting_data.get_last_fitting_data('steps'), 
                                              param_fitting_data.get_last_fitting_data('delta_func')))

               # Evaluate function at new location if not already done previously
               if not cached_data.has_record(new_param_values, 'feval'):
                  if has_custom_data:
                     [new_func_value, custom_object] = self.func(new_param_values, *self.args)
                  else: 
                     [new_func_value, custom_object] = [self.func(new_param_values, *self.args), None]

                  if np.isnan(new_func_value):
                     if param.numerical_error_handler(iiter, param, new_func_value, 'feval') < 0:
                        # Abort fitting
                        self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Function evaluation returned NaN (last value: {2}'.format(
                                                             iiter, param.name, last_func_value))

                        raise Minimizer.ParamNumericalException(iiter, param, new_func_value, 'feval', '')
                     else:

                        # Attempt to fix the problem
                        new_func_value = np.mean(param_fitting_data.get_last_N_fitting_data('feval', self.fitting_options._stop_delta_func_win_size)) # take previous value assumed valid
                        #new_func_value = param_fitting_data.get_last_fitting_data('feval') # take previous value assumed valid

                        warn_msg = ('Iteration {0} - Parameter {1} - Function evaluation returned NaN, reset to: {2} (last evaluation: {3})'.format(
                                   iiter, param.name, new_func_value, last_func_value))
                        if __debug__:
                           if self.is_logging_enabled():
                              self.logger.log_warning(warn_msg)
                           self.fitting_status.add_warn_msg(warn_msg)

                  if np.isinf(new_func_value):
                     if param.numerical_error_handler(iiter, param, new_func_value, 'feval') < 0:
                        # Abort fitting
                        self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - Function evaluation returned Inf (last value: {2})'.format(
                                                         iiter, param.name, last_func_value))
                        raise Minimizer.ParamNumericalException(iiter, param, new_func_value, 'feval', '')
                     else:
                        new_func_value = np.mean(param_fitting_data.get_last_N_fitting_data('feval', self.fitting_options._stop_delta_func_win_size)) # take previous value assumed valid

                        warn_msg = ('Iteration {0} - Parameter {1} - Function evaluation returned Inf, reset to: {2} (last evaluation: {3})'.format(
                                   iiter, param.name, new_func_value, last_func_value))
                        if __debug__:
                           if self.is_logging_enabled():
                              self.logger.log_warning(warn_msg)
                           self.fitting_status.add_warn_msg(warn_msg)

                  param_fitting_data.inc_nb_func_eval()
                  tot_nb_fev += 1
               else:
                  new_func_value  = cached_data.get_value(new_param_values, 'feval')
                  if __debug__:
                        if self.is_logging_enabled(): 
                           self.logger.log_info('Reusing previously-found function value {0} for parameter set {1} ({2})'.format(
                                                new_func_value, new_param_values, param_fitting_data.nb_func_eval))

               # Successive feval difference decreasing count
               direction_count = param_fitting_data.get_last_fitting_data('dcount')

               # Partially differentiate about the parameter
               delta_func = math.fabs(new_func_value / last_func_value) - 1.0 
               ### delta_func = math.fabs(new_func_value) - math.fabs(last_func_value)
               if np.isnan(delta_func):
                  if param.numerical_error_handler(iiter, param, delta_func, 'delta_func') < 0:
                     if __debug__:
                        if self.is_logging_enabled(): 
                           self.logger.log_warning('DeltaFunc is Nan NewFuncValue: {0} LastFuncValue: {1} => aborting fitting'.format(new_func_value, last_func_value))

                     raise ParamNumericalException(iiter, param, delta_func, 'delta_func', '')
                  else:
                     if __debug__:
                        if self.is_logging_enabled(): 
                           self.logger.log_warning('DeltaFunc and NewFuncValue were NaN, NewFuncValue: {0} LastFuncValue: {1}'.format(new_func_value, last_func_value))
                     delta_func = np.mean(param_fitting_data.get_last_N_fitting_data('delta_func', self.fitting_options._stop_delta_func_win_size))
                     self.fitting_status.add_warn_msg('Iteration {0} - Parameter {1} - DeltaFunc and NewFuncValue were NaN, reset to: ({2}, {3}) (last delta: {4})'.format(
                                                      iiter, param.name, delta_func, new_func_value, param_fitting_data.get_last_fitting_data('delta_func')))

               # --- Check if stopping condition on remaining residuals error has been reached
               last_delta_funcs = fitting_data.get_last_N_fitting_data('delta_func_norm', self.fitting_options._stop_delta_func_win_size)   

               #print "*** SCDM: MinimizationType=", self.fitting_options.mini_type

               ###last_func_value  = param_fitting_data.get_last_fitting_data('feval')
               # Check global function evaluation error if specified
               if self.fitting_options.mini_type & MinimizationType.minimize_func_value == MinimizationType.minimize_func_value: 
                  ###if len(last_delta_funcs) >= self.fitting_options._stop_delta_func_win_size and np.mean(last_delta_funcs) < self.fitting_options.max_func_tol:
                  if len(last_delta_funcs) >= self.fitting_options._stop_delta_func_win_size and np.mean(last_delta_funcs) < self.fitting_options.max_func_tol:
                     # We are done, record final results
                     # TODO: may check sign of second order derivative here (should be > 0) to double-check we have indeed reached a minimum

                     if __debug__:
                        if self.is_logging_enabled():
                           self.logger.log_info('*** Found a minimum at value: {0} for parameters: {1} Total tolerance: {2}'.format(
                                                new_func_value, new_param_values, math.fabs(new_func_value - last_func_value)))
                           self.logger.log_info('    [last_func_value: {0} delta_func: {1}] last_delta_funcs: {2} (mean: {3}'.format(
                                                last_func_value, delta_func, last_delta_funcs, np.mean(last_delta_funcs)))

                     # Update result status                  
                     self.fitting_status.set_fitting_status(1, ('Success - Function evaluation tolerance {0} below required limit {1} - Function value: {2}'.format(
                                                              np.mean(last_delta_funcs), self.fitting_options.max_func_tol, new_func_value)) )  # success

                     stop_fitting = True
                     break

               # --- Check sign of delta_func and decide whether to keep the same direction or choose the other one
               if delta_func < 0.0:
            
                  # DeltaFunc is decreasing => keep varying the parameter value along the same direction and we should come closer to the minimum
                  last_delta_func = param_fitting_data.get_last_fitting_data('delta_func')
                  if __debug__:
                     if self.is_logging_enabled():
                        self.logger.log_info('NEGATIVE Derivative, Residuals Decreasing: [delta_func: {0} last_delta_func: {1}]  step: {2} LastDirectionCount: {3}'.format(
                                             delta_func, last_delta_func, step, direction_count))

                  if iiter > 1: 
                     if last_delta_func < 0.0:

                        # Derivative did not change sign => keep parameter varying in the same direction => increase parameter direction count
                        direction_count += 1

                        # Check last decreasing count, i.e. if we should increase the step...
                        if  direction_count > self.fitting_options._dcount_step_inc_threshold:
                           
                           ## If at least the last two iterations were in the same direction, increase the step [multiply <stepSize> by <step_factor>] 
                           if math.fabs(step) > param.step_range[1]:  # check step size upper limit
                              step = param.step_range[1]
                              if __debug__:
                                 if self.is_logging_enabled():
                                    self.logger.log_info('Direction Count: {0} > {1}  Step reached maximum value => Step Unchanged at: {2} ...'.format(
                                                         direction_count, self.fitting_options._dcount_step_inc_threshold, step))
                           else:
                              step *= param.step_factor
                              if __debug__:
                                 if self.is_logging_enabled():
                                    self.logger.log_info('Direction Count: {0} > {1}   [LastDeltaFunc: {2}  NewDeltaFunc: {3}] =>  *** Step increased to: {4}'.format(
                                                         direction_count, self.fitting_options._dcount_step_inc_threshold, last_delta_func, delta_func, step))
                        else:
                           # Keep current step size
                           if __debug__:
                              if self.is_logging_enabled():
                                 self.logger.log_info('DirectionCount: {0} => Step Unchanged at: {1}...'.format(direction_count, step))
                              #pass

                     else:                  
                        # Derivative changed sign...  keep the same direction
                        # TODO: check second derivative here to find out if we are going downhill from a maximum or just going away from a minimum
                        if __debug__:
                           if self.is_logging_enabled():
                              self.logger.log_info('*** Derivative changed sign (NEG => POS)  [newDeltaFunc: {0} last_delta_func: {1}] Step: {2} => keeping varying parameter along same direction'.format(delta_func, last_delta_func, step))

                        # Reset direction count
                        direction_count = 0
                  else:
                     if __debug__:
                        if self.is_logging_enabled():
                           self.logger.log_info('First Iteration: Derivative Decreasing => Keep same direction => Step: {0} DirectionCount: {1}'.format(step, direction_count))
               else: 
                  # Derivative is increasing => we have to the change parameter's direction
                  if __debug__:
                     if self.is_logging_enabled():
                        self.logger.log_info('POSITIVE Derivative, Residuals Increasing: [delta_func: {0}] Step: {1} LastDirectionCount: {2}'.format(
                                             delta_func, step, direction_count))

                  # Check the sign of the function delta         
                  last_delta_func = param_fitting_data.get_last_fitting_data('delta_func')
                  if iiter > 1:
                     if last_delta_func > 0:

                        # We are going towards a maximum Reverse parameter direction
                        if delta_func > last_delta_func:  # safe test as both derivatives are positive
                           # Derivative is positive and increasing => going towards the wrong direction => reversing step
                           if __debug__:
                              if self.is_logging_enabled():
                                self.logger.log_info('Derivative Increasing => *** Reversing direction => Step: ({0} => {1})  DirectionCount: {2}'.format(
                                                     step, -step, direction_count))
                           step = -step
                              
                           # Reset direction count
                           direction_count = 0 
                        else:
                           # Parameter varying along the right direction, keep it
                           if __debug__:
                              if self.is_logging_enabled():
                                self.logger.log_info('Derivative Increasing, direction already reversed => Keeping same direction  DirectionCount: {0}  Step: {1}'.format(
                                                      direction_count, step))

                           # Derivative did not change sign => keep parameter varying in the same direction => increase parameter direction count
                           direction_count += 1
                           #pass

                        # Check if we should increase the step...
                        if direction_count > self.fitting_options._dcount_step_inc_threshold:
                           
                           ## If at least the last two iterations were self._fitting_results = Minimizer._FittingResult(self._fitting_data)in the same direction, increase the step [multiply <stepSize> by <step_factor>] 
                           # Find a smarter algorithm here
                           if math.fabs(step) > param.step_range[1]:  # check step size upper limit
                              step = param.step_range[1]
                              if __debug__:
                                 if self.is_logging_enabled():
                                    self.logger.log_info('Direction Count: {0} > {1} [delta_func: {2} last_delta_func: {3}] Step reached maximum value => Step Unchanged at: {4} ...'.format(
                                                         direction_count, self.fitting_options._dcount_step_inc_threshold, delta_func, last_delta_func, step))
                           else:
                              step *= param.step_factor
                              if __debug__:
                                 if self.is_logging_enabled():
                                    self.logger.log_info('Direction Count: {0} > {1} [LastDeltaFunc: {2}  NewDeltaFunc: {3}]  =>  *** Step increased to: {4}'.format(
                                                         direction_count, self.fitting_options._dcount_step_inc_threshold, last_delta_func, delta_func, step))
                        
                        #else:
                        #   # Keep current step size
                        #   if self.is_logging_enabled():
                        #      self.logger.log_info('Direction Count: %d  Keeping current step: %11.6e' %(direction_count, step)  

                     else:
                        # The previous derivative was negative => we have crossed a temporary minimum => reverse the step direction and reduce the step size [divide <stepSize> by <step_factor>]
                         
                        step /= -(param.step_factor)
                        if __debug__:
                           if self.is_logging_enabled():
                              self.logger.log_info('*** Derivative changed sign (POS => NEG) => Temporary Minimum crossed => Reversing parameter direction and decreasing step to: {0}'.format(step))
                        if math.fabs(step) < param.step_range[0]:  # check step size lower limit
                           #step = param.step_range[0]
                           if __debug__:
                              if self.is_logging_enabled():
                                 self.logger.log_info('Cannot decrease Minimum step size as lower limit reached at: {0:6.3e} => Stop varying parameter '.format(param.step_range[0]))

                           param_fitting_data.stop_varying = True

                        # Reset direction count
                        direction_count = 0
                  else:
                     if __debug__:
                        if self.is_logging_enabled():
                           self.logger.log_info('First Iteration: Derivative Increasing => *** Reversing direction => step: ({0} => {1})  Direction count: {2}'.format(
                                                step, -step, direction_count))
                     step = -step


#               # --- MinimizationType.step_size_lower_bound: check if stopping condition on the sigma of the last N parameters has been reached
#               if self.fitting_options.mini_type & MinimizationType.step_size_lower_bound == MinimizationType.step_size_lower_bound: 

#                     last_step_param_values = param_fitting_data.get_last_N_fitting_data('steps', self.fitting_options._stop_step_param_win_size)

#                     if __debug__:
#                        if self.is_logging_enabled():
#                           self.logger.log_info('last_step_param_values: {0}  Min steps: {1} '.format(
#                                                last_step_param_values, np.all(np.absolute(last_step_param_values) <= param.max_tol) ))                                             
#                     
#                     if len(last_step_param_values) == self.fitting_options._stop_step_param_win_size and np.all(np.absolute(last_step_param_values) <= param.max_tol):
#                        # Take the mean as final value for parameter and stop the fitting
#                        new_param_value = np.mean(last_step_param_values)
#                        new_param_values[param_index] = new_param_value

#                        if __debug__:
#                           if self.is_logging_enabled():
#                              self.logger.log_info('')                                             
#                              self.logger.log_info('*** Average step size: {0} below limit: {1} ==> Will stop varying {2}...'.format(
#                                                    np.mean(last_step_param_values), param.step_range[0], param.name) )
#                              self.logger.log_info('    [last_func_value: {0} delta_func: {1}]'.format(last_func_value, delta_func))
#                              self.logger.log_info('    last_step_param_values: {0}'.format(last_step_param_values))

#                        # Update result status                  
#                        self.fitting_status.set_fitting_status(4, ('Success - Parameter step below required limit {0:6.3e} - Function value: {1}'.format(
#                                                                  param.max_tol, new_func_value)) )  # success

#                        param_fitting_data.stop_varying = True


               # --- MinimizationType.minimize_param_sigma: check if stopping condition on the sigma of the last N parameters has been reached
               if self.fitting_options.mini_type & MinimizationType.minimize_param_sigma == MinimizationType.minimize_param_sigma: 
#               if self.fitting_options.mini_type & MinimizationType.minimize_param_sigma == MinimizationType.minimize_param_sigma: 

                     last_sigma_param_values = np.asarray(param_fitting_data.get_last_N_fitting_data('params', self.fitting_options._stop_sigma_param_win_size))

                     if __debug__:
                        if self.is_logging_enabled():
                           if len(last_sigma_param_values) > 1 and last_sigma_param_values[-1] > 0:
                              self.logger.log_info('last_sigma_param_values: {0}  sigma: {1}'.format(
                                                   last_sigma_param_values/last_sigma_param_values[-1], 
                                                   np.std(last_sigma_param_values/last_sigma_param_values[-1])) )                                             
                     
                     if len(last_sigma_param_values) == self.fitting_options._stop_sigma_param_win_size and \
                        last_sigma_param_values[-1] > 0 and \
                        np.std(last_sigma_param_values/last_sigma_param_values[-1]) <= param.max_tol:
#                        # Take the mean as final value for parameter and stop the fitting
#                        new_param_value = np.mean(last_sigma_param_values)
#                        new_param_values[param_index] = new_param_value

                        if __debug__:
                           if self.is_logging_enabled():
                              self.logger.log_info('')                                             
                              self.logger.log_info('*** Parameter standard deviation: {0} below limit: {1} ==> Will stop varying {2}...'.format(
                                                    np.std(last_sigma_param_values), param.max_tol, param.name) )
                              self.logger.log_info('    [last_func_value: {0} delta_func: {1}]'.format(last_func_value, delta_func))
                              self.logger.log_info('    normalized last_sigma_param_values: {0}'.format(last_sigma_param_values))

                        # Update result status                  
                        self.fitting_status.set_fitting_status(3, ('Success - Parameter sigma below required limit {0:6.3e} - Function value: {1}'.format(
                                                                  param.max_tol, new_func_value)) )  # success

                        param_fitting_data.stop_varying = True


               # --- Check if it is still worth varying current parameter for reducing the residuals: compare the last N function norms
               if self.fitting_options.mini_type & MinimizationType.minimize_param_delta == MinimizationType.minimize_param_delta: 

                  last_param_deltas = np.asarray(param_fitting_data.get_last_N_fitting_data('deval', self.fitting_options._stop_delta_param_win_size))   

                  if __debug__:
                     if self.is_logging_enabled():
                        if len(last_param_deltas) > 1:
                           self.logger.log_info('last_param_delta_values: {0}  mean: {1}'.format(
                                                last_param_deltas, np.mean(last_param_deltas)))

                  if len(last_param_deltas) >= self.fitting_options._stop_delta_param_win_size and \
                     np.mean(last_param_deltas) <= param.max_tol:
                     if __debug__:
                        if self.is_logging_enabled():
                           ###self.logger.log_info('mean(last_param_deltas): %6.3e < limit: %6.3e => *** Will stop varying this parameter...' %(np.mean(last_param_deltas), param.max_tol))
                           self.logger.log_info('')                                             
                           self.logger.log_info('*** mean(last_param_deltas): {0:6.3e} < limit: {1:6.3e} => *** Will stop varying this parameter...'.format(
                                                 np.mean(last_param_deltas), param.max_tol))
                           self.logger.log_info('    [last_func_value:{0} delta_func: {0}]'.format(last_func_value, delta_func))
                           self.logger.log_info('    last_param_deltas: {0}'.format(last_param_deltas))

                     # Update result status                  
                     self.fitting_status.set_fitting_status(2, ('Success - Function variation per parameter below limit {0}  - Function value: {1} '.format(param.max_tol, new_func_value)) ) 
                                              
                     param_fitting_data.stop_varying = True

               # Update Fitted parameters' cache
               cached_data.add_record(new_param_values, new_func_value)

               # Record global fitting data for this iteration and parameter
               fitting_data.record_fitting_data(index=tot_nb_oper, feval=new_func_value, delta_func=delta_func, custom_obj=custom_object)
               ###print 'recorded param data:', 'params:', param_fitting_data.get_last_fitting_data('params'), 'step:', param_fitting_data.get_last_fitting_data('steps'),\
               ###     'deriv:', param_fitting_data.get_last_fitting_data('delta_func'), 'deval:', param_fitting_data.get_last_fitting_data('deval'), 'dcount:', param_fitting_data.get_last_fitting_data('dcount')

               # Record parameter fitting data for this iteration and parameter
               param_fitting_data.record_fitting_data(pval=new_param_value, feval=new_func_value, deval=math.fabs(delta_func), step=step, delta_func=delta_func, dcount=direction_count, iiter=aiter)
                 
               # The new parameter values become the previous one...
               prev_param_values = new_param_values

               # Update latest value residuals to compare during the next parameter iteration
               last_func_value = new_func_value

               # Check if max number of function evaluation exceeded
               stop_fitting = (tot_nb_fev >= self.fitting_options.max_fev)
               if stop_fitting:
                  if __debug__:
                     if self.is_logging_enabled():
                        self.logger.log_info('*** Fitting stopped: Maximum number {0} of function evaluations reached. Residual tolerance: {1}=> Setting status to: -1'.format(tot_nb_fev, new_func_value))
                  self.fitting_status.set_fitting_status(-1, ('Error - Number of function evaluations reached maximum limit {0}'.format(self.fitting_options.max_fev)) )  # error
               else:
                  # Check if none the target parameters are still varying, or if none of the complete set of parameters is varying
                  #print 'STOPPING COND:', not fitting_data.parametersStillVarying(self._target_params), not fitting_data.allParametersStillVarying()
                  stop_fitting = (not fitting_data.parameters_still_varying(self.target_params)) or (not fitting_data.all_parameters_still_varying())
                  if __debug__:
                     if self.is_logging_enabled():
                        if stop_fitting:
                           self.logger.log_info('*** Fitting will stop: all target parameters in: {0} have been found within tolerance limits'.format(self.target_params.param_names))

                  if not set(self.target_params.param_names) == set(self.params.param_names): 
                     if not stop_fitting:
                           # Check if the only parameters still varying are the target parameters. If so, we consider the fit as complete
                           varying_params = fitting_data.get_varying_parameters()
                           stop_fitting = (len(varying_params) == 1 or set(varying_params) == set(self.target_params.param_names))
                           if __debug__:
                              if self.is_logging_enabled():
                                 #if self.is_logging_enabled():
                                 #   self.logger.log_info('*** Parameters {0} are still varying'.format(fitting_data.get_varying_parameters()))
                                 if stop_fitting:
                                    self.logger.log_info('*** Fitting will stop: only target parameters in: {0} are still varying'.format(self.target_params.param_names))
                     

               # ***Note: even when the entire fit is stopping, we complete the whole iteration loop to have a consistent set of recorded data
               #if stop_fitting:
               #   break

               if __debug__:
                  if self.is_logging_enabled():
                     if rank < nb_params: 
                        self.logger.log_info('-' * 177)

            ### end for ###

            if __debug__:
               if self.is_logging_enabled(): 
                  self.logger.log_info('')

            # Check if should stop the fitting process
            if stop_fitting:
               break   

         #### End of fitting process ####

         if iiter >= self.fitting_options.max_iter:
            if __debug__:
               if self.is_logging_enabled(): 
                  self.logger.log_info('**** Fitting stopped: Maximum number {0} of iterations reached. Residual error: {1} => Setting status to: -2'.format(iiter, new_func_value))
               self.fitting_status.set_fitting_status(-2, ('Error - Number of iterations reached maximum limit {0}'.format(self.fitting_options.max_iter)))

         # Record final number of iteration 
         #self.nb_iter(iiter)
         
         
         self._set_nb_iter(iiter)

         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_info('=' * 175)

         return [new_param_values, self.fitting_results]      # return fitted set of parameters

      except Minimizer.ParamOutOfBoundException:
         self.fitting_status.set_fitting_status(-3, sys.exc_info()[1])
         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_error('Fitting stopped: {0}'.format(sys.exc_info()[1]))
         return [prev_param_values, self.fitting_results]            

      except Minimizer.ParamNumericalException:
         self.fitting_status.set_fitting_status(-4, sys.exc_info()[1])
         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_error('Fitting stopped: {0}'.format(sys.exc_info()[1]))
         return [prev_param_values, self.fitting_results]            

#      except ParamOutOfBoundException as exc:
#         self.fitting_status.set_fitting_status(-3, str(exc))
#         if __debug__:
#            if self.is_logging_enabled(): self._fitting_results = Minimizer._FittingResult(self._fitting_data)
#               self.logger.log_error('Fitting stopped: %s' %str(exc))
#         return [prev_param_values, self.fitting_results]            

#      except ParamNumericalExceptionas exc:
#         self.fitting_status.set_fitting_status(-4, str(exc))
#         if __debug__:
#            if self.is_logging_enabled(): 
#               self.logger.log_error('Fitting stopped: %s' %str(exc))
#         return [prev_param_values, self.fitting_results]            

      except OverflowError:
         self.fitting_status.set_fitting_status(-10, 'Numerical Error - Numerical error ({0})'.format(sys.exc_info()[1]))
         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_error('Fitting stopped: Numerical error ({0})'.format(sys.exc_info()[1]))
         return [prev_param_values, self.fitting_results]            
   
      except (ValueError, KeyError, IndexError):

         exc_type, exc_value, exc_traceback = sys.exc_info()
         #traceback.print_tb(exc_traceback, limit=1, file=self.logger.getFileDescr())

         self.fitting_status.set_fitting_status(-11, 'Value Error - ({0} {1})'.format(exc_type, exc_value))
         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_error('Fitting stopped: error ({0}) encountered'.format(sys.exc_info()[1]))
         return [prev_param_values, self.fitting_results]            

      except MemoryError:
         self.fitting_status.set_fitting_status(-12, 'Memory Error - {0} '.format(sys.exc_info()[1]))
         if __debug__:
            if self.is_logging_enabled(): 
               self.logger.log_error('Fitting stopped: Memory error ({0}) encountered'.format(sys.exc_info()[1]))
         return [prev_param_values, self.fitting_results]            
         
### DISABLED 

#      except Exception:
#         self.fitting_status.set_fitting_status(-13, 'Misc Error - {0}'.format(sys.exc_info()[1]))
#         if __debug__:
#            if self.is_logging_enabled(): 
#               self.logger.log_error('Fitting stopped: error ({0}) encountered'.format(sys.exc_info()[1]))
#         return [prev_param_values, self.fitting_results]            
         

   # ~~~~~~~~~~~~~~~~~~~~~
   # Fit all parameters
   # ~~~~~~~~~~~~~~~~~~~~~
   def minimize(self):

      if __debug__:

         # Logs and diagnostic options
         if self.is_logging_enabled():
            started = time.clock()      
            log_directory = self.get_log_directory()
            if not os.path.exists(log_directory):
               os.mkdir(log_directory)
            self.logger = slogger.FileLogger(logger='', directory=log_directory, filename=self.get_log_filename())

         if self.is_stats_enabled():
            started = time.clock()      
            self._stats_helper = self.StatsHelper(self, self.get_stats_directory())

         if self.is_plotting_enabled():
            self._plot_helper  = self.PlotHelper(self, self.get_plot_directory())

      # Fit each parameter in turn (a parallel version should be better...)
      [fitted_param_values, fitting_results] = self.fit_parameters()

#      print "***SCDM:", [fitted_param_values, fitting_results]

      # Store fitted parameter values
      fitting_results.fitted_param_values = fitted_param_values

      if __debug__:

         if self.is_stats_enabled() or self.is_logging_enabled():
            elapsed = time.clock() - started      

         # Display stats
         if self.is_stats_enabled():

            # Output some statistics
            self._stats_helper.output_parameter_stats()
      
            # Output summary information
            print('\n')
            print('-' * 125)
            print('Status: {0} {1}'.format(self.fitting_status.errno, self.fitting_status.msg))
            for msg in self.fitting_status.warn_msgs:
               print('*** Warning: {0}'.format(msg))
            print('Fitted Parameters: {0}'.format(fitted_param_values))
            print('Total iterations: {0}  Total function evaluations: {1}'.format(self.get_nb_iter(), self.get_nb_func_eval())    )  
            print('Fitting time: {0:6.2f} sec'.format(elapsed))                     
            print('-' * 125)

         # Log stats
         if self.is_logging_enabled():
            self.logger.log_info('Minimization ended. Status: {0} {1}'.format(self.fitting_status.errno, self.fitting_status.msg))

            self.logger.log_blank_line()
            self.logger.log_info('=' * 177)
            self.logger.log_info('Status: {0} {1}'.format(self.fitting_status.errno, self.fitting_status.msg))
            for msg in self.fitting_status.warn_msgs:
               self.logger.log_info('*** Warning: {0}'.format(msg))
            self.logger.log_info('Fitted Parameters: {0}'.format(fitted_param_values))
            self.logger.log_info('Total iterations: {0}  Total function evaluations: {1}'.format(self.get_nb_iter(), self.get_nb_func_eval()))   
            self.logger.log_info('Fitting time {0:6.2f} sec'.format(elapsed))                     
            self.logger.log_info('=' * 177)

         # Plot data      
         if self.is_plotting_enabled():

            # Plot global data (aggregated data for all parameters)   
            self._plot_helper.plot_global_fitting_data()
            
            # Plot individual parameter data
            self._plot_helper.plot_parameter_fitting_data()

         if self.is_logging_enabled() and self.logger.is_open():
            self.logger.flush()
            self.logger.close()
   

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Returns: 
      # - list of fitted parameters
      # - dictionaries with additional fitting information
      # - fitting-related data: residuals, etc.
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      return [fitted_param_values, fitting_results]

   # ~~~~~~~~~~~~~~~
   # Private methods
   # ~~~~~~~~~~~~~~~

   def _setup_options(self):
      """! Setup fitting and diagnostic optiopns """    

      # If diagnostic/tuning required we keep the entire fitting history, otherwise, 
      #  just keep minimum history data to save memory
      if __debug__:

         if (self.is_stats_enabled() or self.is_plotting_enabled()):
            # Keep all fitting data for all parameters
            self._fitting_options.data_window_size = -1
            self._fitting_data  = self._FittingData(self._params, -1)    

         else:         
            # Keep an history of size <data_window_size
            self._fitting_data  = self._FittingData(self._params, 
                                                    self.fitting_options._data_window_size)         

      else:         
         # Keep an history of size <data_window_size> 
         self._fitting_data  = self._FittingData(self._params, 
                                           self.fitting_options._data_window_size)           

      # Base Directory      
      self._base_output_dir = self.diag_options.base_output_dir
      if not os.path.exists(self._base_output_dir):
         os.mkdir(self._base_output_dir)

   def _dummy_func(self, parameters):   
      """! Dummy Evaluation function: used if no eval function provided """
      #print('Please provide an evaluation function for fitting parameters: {0}'.format(parameters.param_names))

      
   def _set_nb_iter(self, nb_iter):
      if self._fitting_results is None:
         self._fitting_results = Minimizer._FittingResult(self._fitting_data)
      self._fitting_results.nb_iter = nb_iter

   # ~~~~~~~~~~~~~~~~~~~
   # Fitting Exceptions
   # ~~~~~~~~~~~~~~~~~~~

   class FittingException(Exception):
      """! Base fitting exception """

      def __init__(self, iiter, msg='FittingException: error: {0}'.format(sys.exc_info()[1])):
        self.iiter = iiter
        self.msg = msg
      def __str__(self):
         return 'FittingException: iteration {0} - error: {1} ({2})'.format(self.iiter, self.msg, sys.exc_info()[1])

   class ParamOutOfBoundException(FittingException):
      """! Exception thrown in case of out of parameter bounds access """

      def __init__(self, iiter, param, param_value, bound_index, msg='ParamOutOfBoundException: error: {0}'.format(sys.exc_info()[1])):
        FittingException.__init__(self, iiter, msg)
        self.param = param 
        self.param_value = param_value
        self.bound_index = bound_index 
      def __str__(self):
         if self.bound_index == 0:  # lower bound reached
            return 'ParamOutOfBoundException: Iteration {0} - Parameter {1} - Lower bound {2:6.3e} of parameter was reached ({3:6.3e}) {4}'.format(
                    self.iiter, self.param.name, self.param.param_range[self.bound_index], self.param_value, self.msg)
         elif self.bound_index == 1: # upper bound reached  
            return 'ParamOutOfBoundException: Iteration {0} - Parameter {1} - Upper bound {2:6.3e} of parameter was reached ({3:6.3e}) {4}'.format(
                    self.iiter, self.param.name, self.param.param_range[self.bound_index], self.param_value, self.msg)

   class ParamNumericalException(FittingException):
      """! Exception thrown in case numerical error """

      def __init__(self, iiter, param, value, value_name, msg='ParamNumericalException: error: {0}'.format(sys.exc_info()[1])):
        Minimizer.FittingException.__init__(self, iiter, msg)
        self.param = param 
        self.value = value
        self.value_name = value_name
      def __str__(self):
         if np.isnan(self.value):
            return 'ParamNumericalException: Iteration {0} - Parameter {1} - {2} is NaN {3}'.format(self.iiter, self.param.name, self.value_name, self.msg)
         elif np.isinf(self.value): 
            return 'ParamNumericalException: Iteration {0} - Parameter {1} - {2} is Inf {3}'.format(self.iiter, self.param.name, self.value_name, self.msg)



   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Statistics Helper class
   # - Public inner class of Minimizer 
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class StatsHelper(object):

      """!
         Output statistics related to Parameter fitting. 
         A inner class of the Minimizer class 
      """
      def __init__(self, minimizer, output_directory='./stats'):

         self._minimizer = minimizer
         if self._minimizer.is_logging_enabled:
            self._logger = self.minimizer.logger         
         self._output_directory = output_directory

         #if not os.path.exists(output_directory):
         #   os.mkdir(output_directory)

      # ~~~~~~~~~~~
      # Properties 
      # ~~~~~~~~~~~

      # --- Getters ---

      @property
      def minimizer(self):
         return self._minimizer

      @property
      def logger(self):
         return self._logger

      @property
      def output_directory(self):
         return self._output_directory

      # ~~~~~~~~~~~~~~~
      # Public Methods 
      # ~~~~~~~~~~~~~~~

      # Output statistics related to Parameter fitting
      def output_parameter_stats(self):

         if self.minimizer.is_diagnose:
            self.logger.log_blank_line()
            self.logger.log_info('=' * 177)

         for param_name in self.minimizer.params.param_names:

            fitting_data = self.minimizer.fitting_results.fitting_data
            if fitting_data.has_data():
               param_fitting_data = fitting_data.get_param_fitting_data_for(param_name)
               param_fitting_data_dico = param_fitting_data.get_fitting_data_dico()

               if self.minimizer.is_diagnose:
                  print('\n*** Statistics for parameter: {0} ***'.format(param_name))
                  self.logger.log_info('*** Statistics for parameter: {0} ***'.format(param_name))

               for key in param_fitting_data_dico.keys():
                  if len(param_fitting_data_dico[key])  > 0:
                     if self.minimizer.is_diagnose:
                        print('{0}\t- Avg:\t{1:9.6e}\t\tmin: {2:9.6e}\t\tMax:\t{3:9.6e}'.format(
                              key, np.mean(param_fitting_data_dico[key]), 
                              min(param_fitting_data_dico[key]), max(param_fitting_data_dico[key]) ))
                        self.logger.log_info('{0} - Avg:\t{1:9.6e}\t\tmin: {2:9.6e}\t\tMax:\t{3:9.6e}'.format(
                              key, np.mean(param_fitting_data_dico[key]), 
                              min(param_fitting_data_dico[key]), max(param_fitting_data_dico[key])))
               if self.minimizer.is_diagnose:
                  self.logger.log_info('-' * 177)

         if self.minimizer.is_diagnose:
            self.logger.log_info('=' * 177)


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #  Plotting helper class  
   # - Public inner class of Minimizer 
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class PlotHelper(object):

      """!
         Helper class for plotting fitting results.
         This is is a inner class of the Minimizer class.  
      """

      def __init__(self, minimizer, output_directory='./plots'):

         self.minimizer = minimizer
         if self.minimizer.is_logging_enabled:
            self.logger = self.minimizer.logger        
         self.output_directory = output_directory

         if self.minimizer.is_plotting_enabled:
            if not os.path.exists(output_directory):
               os.mkdir(output_directory)

      # ~~~~~~~~~~~~~
      #  plot stamp  
      # ~~~~~~~~~~~~~
      def plot_stamp(self, stamp, plot_title='Stamp', has_grid=False, output_dir='.', output_file='stamp_plot', show=False):

         """! plot a stamp in 2D """

         try:
            height, width = stamp.shape

            figure = matplotlib.figure.Figure()
            figtitle = plot_title
            pylab.gcf().text(0.5, 0.93, plot_title, horizontalalignment='center', fontproperties=FontProperties(size=12, weight='bold'))
            pylab.grid(has_grid)
            axes = figure.gca() 
            axes.set_aspect('equal')
            pylab.xlim(1, width)
            pylab.ylim(1, height)
            
            x = np.arange(0, width+1)   
            y = np.arange(0, height+1)   
            plot = pylab.pcolor(x, y, stamp, shading='flat')
            pylab.colorbar(plot, extend='both', shrink=0.8)
            #pylab.pcolor(stamp, edgecolors=None, shading='flat')
            
            pylab.savefig(os.path.join(output_dir, output_file))

            if show:
               pylab.show()

         except (OverflowError, ValueError, IndexError):        
               print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
               #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))# Misc modules

         finally:
            pylab.clf()
            pylab.close()


      # ~~~~~~~~~~~~~~~~~~~
      # plot contour stamp
      # ~~~~~~~~~~~~~~~~~~~
      def plot_contour_stamp(self, stamp, plot_title='Stamp contour', has_grid=False, output_dir='.', output_file='stamp_contour_plot', show=False):

         """!plot stamp contours in 3D """

         try:

            stamp_size, stamp_size = stamp.shape

            figure = pylab.figure()
            figtitle = plot_title
            pylab.gcf().text(0.5, 0.93, plot_title, horizontalalignment='center', fontproperties=FontProperties(size=12, weight='bold'))
            pylab.grid(has_grid)
            axes = figure.gca() 
            axes.set_aspect('equal')
            pylab.xlim(1, stamp_size)
            pylab.ylim(1, stamp_size)
            
            lin_space = np.linspace(0, stamp_size -1, stamp_size)
            x, y = np.meshgrid(lin_space, lin_space)

            # draw contour plot   
            ctplot = pylab.contour(x, y, stamp)
         #   pylab.clabel(ctplot,
         #               inline=1,
         #               fmt='%1.1f',
         #               fontsize=14)
            
            # make a colorbar for the contour lines
            #cbar_contour = pylab.colorbar(ctplot, shrink=0.8, extend='both')
            
            im = pylab.imshow(stamp, interpolation='bilinear', origin='lower')

            # We can still add a colorbar for the image, too.
            cbar_image = pylab.colorbar(im, orientation='vertical', shrink=0.8)


            pylab.savefig(os.path.join(output_dir, output_file))

            if show:
               pylab.show()

         except (OverflowError, ValueError, IndexError):        
               print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
               #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))# Misc modules

         finally:
            pylab.clf()
            pylab.close()

      # ~~~~~~~~~~~~~~~~~~~
      # Plot a Stamp in 3D
      # ~~~~~~~~~~~~~~~~~~~
      def plot_stamp_3D(self, stamp, plot_title='Stamp 3D', x_label=None, y_label=None, z_label=None, has_grid=False, output_dir='.', output_file='stamp_plot_3D', show=False):

         """! plot a stamp in 3D """

         try:
        
            stamp_size, stamp_size = stamp.shape

            figure = pylab.figure()
            figtitle = plot_title
            pylab.gcf().text(0.5, 0.93, plot_title, horizontalalignment='center', fontproperties=FontProperties(size=12, weight='bold'))
            pylab.grid(has_grid)
         #   axes = figure.gca() 
         #   axes.set_aspect('equal')
            
            lin_space = np.linspace(0, stamp_size -1, int(stamp_size))
            x, y = np.meshgrid(lin_space, lin_space)
            z = stamp

            
            #ax = matplotlib.axes3d.Axes3D(figure)
            ax = pylab3.Axes3D(figure)
            #surface = ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='gist_heat')
            surface = ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap=pylab.cm.jet)

            if not x_label is None:
               ax.set_xlabel(x_label, size='medium', weight='heavy', color='red')
            if not y_label is None:
               ax.set_ylabel(y_label, size='medium', weight='heavy', color='red')
            if not z_label is None:
               ax.set_zlabel(z_label, size='medium', weight='heavy', color='red')

            #figure.colorbar(surface, shrink=0.5, aspect=5)

            pylab.savefig(os.path.join(output_dir, output_file))

            if show:
               pylab.show()
      
         except (OverflowError, ValueError, IndexError):        
               print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
               #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))# Misc modules

         finally:
            pylab.clf()
            pylab.close()


      # ~~~~~~~~~~~~~~~~~~~
      # Draw an array plot
      # ~~~~~~~~~~~~~~~~~~~
      def plot_arrays(self, data_dico, value_names,
                      plot_filename, output_directory, plot_title='', pfilter=None, first_index=-1, last_index=-1,
                      x_label=None, y_label=None, has_grid=True,pfigure=None,
                      pcolor='blue', pmarker='', pmarkersize=1.0, pmarkerfacecolor='red', plinestyle='-', plinewidth=1, is_normed=False, has_equalaxes=False,
                      has_lstsqfit=False, pfitdegree=2,
                      pfitcolor='green', pfitmarker='', pfitmarkersize=1.0, pfitmarkerfacecolor='black', pfitlinewidth=2, pfitlinestyle='-.', show=False):
                             
         try:

            if pfigure == None:
               fig = pylab.figure()
            else:
               fig = pfigure
            pylab.grid(has_grid)
            if has_equalaxes:
               axes = fig.gca() 
               axes.set_aspect('equal')

            iplot = 1

            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = 'Line plot of %s: ' %(value_names)

            # Two values to represent
            if value_names[0] in data_dico and value_names[1] in data_dico:
               values_x = np.array(data_dico[value_names[0]])
               values_y = np.array(data_dico[value_names[1]])

               values_x = np.nan_to_num(values_x)
               values_y = np.nan_to_num(values_y)

               if pfilter != None:
                  filtered_values_x = values_x[pfilter]
                  filtered_values_y = values_y[pfilter]
                  if len(filtered_values_x) > 0:
                     values_x = filtered_values_x
                  if len(filtered_values_y) > 0:
                     values_y = filtered_values_y

               min_index = 0
               max_index = len(values_x)
               if first_index > -1 and first_index < max_index:
                  min_index = first_index

               if last_index > -1 and last_index < max_index:
                  max_index = last_index

               values_x = values_x[min_index:max_index]
               values_y = values_y[min_index:max_index]     

               #print 'values_x', values_x
               #print 'values_y', values_y

            #   if first > -1 and len(values_x) > max_len and len(values_y) > max_len :
            #      if not last_values:
            #         # Take first  <max_len> values
            #         values_x = values_x[0:max_len]
            #         values_y = values_y[0:max_len]     
            #      else:
            #         # Take last  <max_len> values
            #         values_x = values_x[-max_len:len(values_x)]
            #         values_y = values_y[-max_len:len(values_y)]     

               # Plot the data
               pylab.plot(values_x, values_y,  linestyle=plinestyle, color=pcolor,\
                          marker=pmarker, markersize=pmarkersize, markerfacecolor=pmarkerfacecolor, linewidth=plinewidth)

               # Also include a least-squares fit if requested
               #Q_factor_calc_freq = string.atof(config.get_param_value('Q_FACTOR_CALC_FREQ'))
               if has_lstsqfit:
                  if len(values_x) > 20:
                     #if has_lstsqfit and len(values_x) > Q_factor_calc_freq:
                     p = np.polyfit(values_x, values_y, pfitdegree)
                     fitted_values_y = np.polyval(p, values_x)
                     
                     values_x = np.extract(fitted_values_y >=0, values_x)
                     fitted_values_y = np.extract(fitted_values_y >=0, fitted_values_y)
                     str_p = str([string.strip('%7.3f' %(e)) for e in p])
                     fit_info = '(degree ' + str(pfitdegree) + ' lstsq fit: ' + str_p + ')'
               
                     # Plot the data
                     pylab.plot(values_x, fitted_values_y, linestyle=pfitlinestyle, color=pfitcolor,\
                                marker=pfitmarker, markersize=pfitmarkersize, markerfacecolor=pfitmarkerfacecolor, linewidth=pfitlinewidth)
                  else:
                     fit_info = ''
                     print 'Warning: Not enough values to fit the data'
               else:
                  fit_info = ''

               if not x_label is None:
                  pylab.xlabel(x_label)
               if not y_label is None:
                  pylab.ylabel(y_label)    

               figtitle += '\n' + fit_info     

               pylab.gcf().text(0.5, 0.93, figtitle, horizontalalignment='center', fontproperties=FontProperties(size=12, weight='bold'))

               pylab.savefig(os.path.join(output_directory, plot_filename))

               if show:
                  pylab.show()

            else:
               print 'Error: %s or %s not found in data dictionary' %(value_names[0], value_names[1])

         except (OverflowError, ValueError, IndexError):        
               print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
               #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))# Misc modules

         finally:
            pylab.clf()
            pylab.close()

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Plot a set of arrays in 3D
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~
      def plot_arrays_3D(self, data_dico, value_names,
                       plot_filename, output_directory, plot_title='', extra_title=None, pfilter=None, first_index=-1, last_index=-1,
                       x_label=None, y_label=None, z_label=None, has_grid=True,pfigure=None,
                       is_normed=False, has_equalaxes=False,
                       show=False):
                             
         try:

            if pfigure == None:
               fig = pylab.figure()
            else:
               fig = pfigure
            pylab.grid(has_grid)
            if has_equalaxes:
               axes = fig.gca() 
               axes.set_aspect('equal')

            iplot = 1

            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = 'Surface plot of %s: ' %(value_names)

            # Two values to represent
            if value_names[0] in data_dico and value_names[1] and value_names[2] in data_dico:
               values_x = np.array(data_dico[value_names[0]])
               values_y = np.array(data_dico[value_names[1]])
               values_z = np.array(data_dico[value_names[2]])

               values_x = np.nan_to_num(values_x)
               values_y = np.nan_to_num(values_y)
               values_z = np.nan_to_num(values_z)

               if pfilter != None:
                  filtered_values_x = values_x[pfilter]
                  filtered_values_y = values_y[pfilter]
                  filtered_values_z = values_y[pfilter]
                  if len(filtered_values_x) > 0:
                     values_x = filtered_values_x
                  if len(filtered_values_y) > 0:
                     values_y = filtered_values_y
                  if len(filtered_values_z) > 0:
                     values_z = filtered_values_z

               min_index = 0
               max_index = len(values_x)
               if first_index > -1 and first_index < max_index:
                  min_index = first_index

               if last_index > -1 and last_index < max_index:
                  max_index = last_index

               values_x = values_x[min_index:max_index]
               values_y = values_y[min_index:max_index]     
               values_z = values_z[min_index:max_index]     


               # Points must be equally spaced... so have to interpolate to fill missing points
               #spline = scipy.interpolate.Rbf(values_x, values_y, values_z,function='thin-plate')
               spline = scipy.interpolate.Rbf(values_x, values_y, values_z,function='linear', smooth=1)
               xi = np.linspace(values_x.min(), values_x.max(), len(values_x))
               yi = np.linspace(values_y.min(), values_y.max(), len(values_y))
               xi, yi = np.meshgrid(xi, yi)
               zi = spline(xi,yi)  # spline interpolation          

               # Plot the data
               ax = pylab3.Axes3D(fig)
               #ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap='Oranges_r')
               ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=pylab.cm.jet)

               ###pylab.setp(ax.get_xticklabels(), family='sans-serif')
               ###pylab.setp(ax.get_xticklabels(), fontname='Lucida Grande')
               ###pylab.setp(ax.get_xticklabels(), size='4' )

               if not x_label is None:
                  ax.set_xlabel('%s' %value_names[0], size='medium', weight='heavy', color='red')
               if not y_label is None:
                  ax.set_ylabel('%s' %value_names[1], size='medium', weight='heavy', color='red')
               if not z_label is None:
                  ax.set_zlabel('%s' %value_names[2], size='medium', weight='heavy', color='red')


               #pylab.xlabel(x_label)
               #pylab.ylabel(y_label)    

               # Title and Extra title
               fig.text(0.5, 0.93, figtitle, horizontalalignment='center', fontproperties=FontProperties(size=12, weight='bold'))
               if not extra_title is None:
                  pylab.gcf().text(0.05, 0.10, extra_title, fontproperties=FontProperties(size=10), rotation='vertical', color='red')
                  #fig.text(0.0, 1.0, extra_title, horizontalalignment='center', fontproperties=FontProperties(size=10, rotation='vertical', facecolor='red'))

               pylab.savefig(os.path.join(output_directory, plot_filename))

               if show:
                  pylab.show()

            else:
               print 'Error: %s, %s or %s not found in data dictionary' %(value_names[0], value_names[1], value_names[2])

         except (OverflowError, ValueError, IndexError):        
               print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
               #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))# Misc modules

         finally:
            pylab.clf()
            pylab.close()


         
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 3D plots of a set of 3 parameters
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      def plot_parameter_set_3D(self, param_names, dataType,
                                plot_filename, output_directory, plot_title='', extra_title=None, pfilter=None, first_index=-1, last_index=-1,
                                x_label=None, y_label=None, z_label=None, has_grid=True, pfigure=None, is_normed=False, has_equalaxes=False, show=False):
         
         # Collect parmeter date for dataType <dataType> in one single dictionary 
         param_data_dico = {}
         for param_name in param_names:
            fitting_data = self.minimizer.fitting_results.fitting_data
            if fitting_data.has_Data():
               param_fitting_data = fitting_data.get_param_fitting_data_for(param_name)
               param_fitting_data_dico = param_fitting_data.fitting_data_dico
               param_data_dico[param_name] = param_fitting_data_dico[dataType]

         # Plot the 3-parameter surface      
         self.plot_arrays_3D(param_data_dico, param_names,
                             plot_filename, output_directory, plot_title, extra_title, pfilter, first_index, last_index,
                             x_label, y_label, z_label, has_grid, pfigure, is_normed, has_equalaxes, show)   

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # 2D Plots for all parameters
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
      def plot_global_fitting_data(self, first_array_index=0, last_array_index=-1):

         fitting_data = self.minimizer.fitting_results.fitting_data
         if fitting_data.has_data():
            fitting_data_dico = fitting_data.fitting_data_dico

            plot_title_prefix = ('All parameters - ')
            plot_filename_prefix = ('global')

            self.plot_arrays(fitting_data_dico, ['index', 'feval'],\
                     plot_filename_prefix + '_residuals_all_parameters', self.output_directory, plot_title_prefix + 'Evolution residuals (all parameters)',\
                     first_index=first_array_index, last_index=last_array_index,\
                     x_label='iterations', y_label='residuals',\
                     pcolor='#6a5acd', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

            self.plot_arrays(fitting_data_dico, ['index', 'delta_func'],\
                     plot_filename_prefix + 'deltafunc_iter_all_parameters', self.output_directory, plot_title_prefix + 'Evolution of function variations over step size (delta_func)',\
                     first_index=first_array_index, last_index=last_array_index,\
                     x_label='iterations', y_label='residuals',\
                     pcolor='#ff8c00', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)


         else:
            print 'Info: no data are available for plotting'


      # Per parameter plots 
      def plot_parameter_fitting_data(self, first_array_index=0, last_array_index=-1):

         for param_name in self.minimizer.params.param_names:

            fitting_data = self.minimizer.fitting_results.fitting_data
            if fitting_data.has_data():


               param_fitting_data = fitting_data.get_param_fitting_data_for(param_name)
               param_fitting_data_dico = param_fitting_data.get_fitting_data_dico()

               print '\n*** Parameter %s ***' %param_name

               # Averages   
               for key in param_fitting_data_dico.keys():
                  if len(param_fitting_data_dico[key]) > 0:
                     print 'parameter %s - Avg:\t%9.6e\t\tmin: %9.6e\t\tMax:\t%9.6e' %(key, np.mean(param_fitting_data_dico[key]), min(param_fitting_data_dico[key]), max(param_fitting_data_dico[key]))

               print '\nMaking plots for parameter %s in directory %s...' %(param_name, self.output_directory)   

               param_plot_title_prefix = ('Parameter %s - ' %(param_name))
               param_plot_filename_prefix = ('param_%s' %(param_name))

               if len(param_fitting_data_dico['params']) > 0:

                  # Plot evolution of estimated parameter value
                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'params'],\
                           param_plot_filename_prefix + '_estimated_value_iter', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of parameter estimates [%6.3e]' %(param_fitting_data_dico['params'][-1]), \
                           first_index=first_array_index, last_index=last_array_index,\
                           x_label='iterations', y_label='parameter estimates',\
                           pcolor='blue', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)
                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'params'],\
                           param_plot_filename_prefix + '_estimated_value_iter_2', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of parameter estimates [%6.3e]' %(param_fitting_data_dico['params'][-1]),
                           first_index=100, last_index=last_array_index,\
                           x_label='iterations', y_label='parameter estimates',\
                           pcolor='blue', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  self.plot_arrays(param_fitting_data_dico, ['params', 'feval'],\
                           param_plot_filename_prefix + '_sum_residuals_param', self.output_directory, 
                           param_plot_title_prefix + 'Evolution residuals with parameter estimates [%6.3e]' %(param_fitting_data_dico['feval'][-1]),\
                           first_index=first_array_index, last_index=last_array_index,\
                           x_label='parameter estimate', y_label='residuals',\
                           pcolor='#8b8b83', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  # Plot evolution of residuals per parameter value
                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'feval'],\
                           param_plot_filename_prefix + '_sum_residuals_iter', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of residuals per parameter [%6.3e]' %(param_fitting_data_dico['feval'][-1]),\
                           first_index=first_array_index, last_index=100,\
                           x_label='iterations', y_label='residuals',\
                           pcolor='#CD3333', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  # Plot evolution of step sizes
                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'steps'],\
                           param_plot_filename_prefix + '_step_iter', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of step values [%6.3e]' %(param_fitting_data_dico['steps'][-1]),\
                           first_index=first_array_index, last_index=last_array_index,\
                           x_label='iterations', y_label='step size',\
                           pcolor='#548B54', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'deval'],\
                           param_plot_filename_prefix + '_deval_iter', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of function variation over step size [%6.3e]' %(param_fitting_data_dico['deval'][-1]),\
                           first_index=1, last_index=last_array_index,\
                           x_label='iterations', y_label='function variation over step size',\
                           pcolor='#8E2323', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'delta_func'],\
                           param_plot_filename_prefix + '_delta_func_iter', self.output_directory, 
                           param_plot_title_prefix + 'Evolution of function variations over step size [%6.3e]' %(param_fitting_data_dico['delta_func'][-1]),\
                           first_index=1, last_index=last_array_index,\
                           x_label='iterations', y_label=' estimate',\
                           pcolor='#ff8c00', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)

                  self.plot_arrays(param_fitting_data_dico, ['iiter', 'dcount'],\
                           param_plot_filename_prefix + '_direction_count', self.output_directory, 
                           param_plot_title_prefix + 'Evolution direction count [%6.3e]' %(param_fitting_data_dico['dcount'][-1]),\
                           first_index=first_array_index, last_index=last_array_index,\
                           x_label='iterations', y_label='direction count',\
                           pcolor='#8a2be2', pmarker='', pmarkersize=2.5, is_normed=False, has_equalaxes=False, has_grid=True, has_lstsqfit=False)
               else:
                  print 'Info: no data are available for plotting for parameter %s' %(param_name)
            else:
               print 'Info: no data are available for plotting for parameter %s' %(param_name)

      def showPlots(self):
         pylab.show()


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Fitting Cache with tuple of parameters [p0, ... pn] as key
   # - Inner class of Minimizer 
   # - A private class that should not be accessed directly
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class _FittingCache(object):
      """! 
          Private class for managing a cache with tuples of parameters [p0,...,pn] as keys
      """   

      def __init__(self):
         self._fitting_cache_dico = {}
         for value_type in ['feval']:
            self._fitting_cache_dico[value_type] = {} 

      # ~~~~~~~~~~~~~~~
      # Public methods
      # ~~~~~~~~~~~~~~~

      def is_empty(self):
         return (self.get_size() == 0)

      def clear(self):
         del self._fitting_cache_dico
         self.__init__()      

      def get_size(self):
         return len(self._fitting_cache_dico['feval'])

      def has_record(self, param_values, value_type):
         return tuple(param_values) in self._fitting_cache_dico[value_type]
         
      def add_record(self, param_values, feval):
         self._fitting_cache_dico['feval'][tuple(param_values)] = feval

      def get_value(self, param_values, value_type):
         try:
            return self._fitting_cache_dico[value_type][tuple(param_values)] 
         except KeyError:
            print('Error: cound not find data for {0} in cache'.format(param_values))


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Hold global fitting data for all parameters
   # - Inner class of Minimizer 
   # - A private class that should not be accessed directly
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class _FittingData(object):
      """! 
         Hold global fitting data for all parameters. 
         This is a inner class of the Minimizer class.
      """

      def __init__(self, params, data_window_size=-1):

         # Global data
         self._fitting_data_dico = {}
         for value_type in ['index', 'feval', 'delta_func', 'delta_func_norm', 'custom_objs']:
            self._fitting_data_dico[value_type] = [] 

         # Parameter data
         self._params = params
         self._param_data_dico = {}  # dictionary of ParamFittingData objects with parameter name as key
         self._param_names = []      # ordered list of parameter names       
         for param in self._params.param_list:
            self._param_data_dico[param.name] = self._ParamFittingData(
                                                                param, 
                                                                data_window_size=data_window_size)
            self._param_names.append(param.name)
         #self._param_names.sort()

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~

      @property
      def fitting_data_dico(self):
         """! Return the global fitting data dictionary """
         return self._fitting_data_dico

      @property
      def param_data_dico(self):
         """! Return the param data dictionary """
         return self._param_data_dico

      @property
      def param_names(self):
         """! Return the list of param names """
         return self._param_names

      # ~~~~~~~~~~~~~~~
      # Public methods
      # ~~~~~~~~~~~~~~~

      def clear_data(self):
         # Clear global data
         self._fitting_data_dico.clear()
         for value_type in ['index', 'feval', 'delta_func', 'delta_func_norm', 'custom_objs']:
            self._fitting_data_dico[value_type] = []
            
         # Clear  parameter data
         for param_name in self._param_names:
            self._param_data_dico[param_name].clear_data()      
         
      def set_param_fitting_data_for(self, param_name, param_data):
         self._param_data_dico[param_name] = param_data

      def get_param_fitting_data_for(self, param_name):
         return self._param_data_dico[param_name]

      def get_param_fitting_data(self):
         return [self._param_data_dico[p.name] for p in self._params.param_list]

      def has_data(self):
         return len(self._param_names) > 0 and len(self._param_data_dico) > 0 and self._param_data_dico[self._param_names[0]].has_records()
         #return len(self._param_names) > 0 and len(self._param_data_dico[self._param_names[0]]) > 0
      
      def get_fitting_data(self, data_type):
         return self._fitting_data_dico[data_type]

      def get_last_fitting_data(self, data_type):
         try:
            return self._fitting_data_dico[data_type][-1]
         except IndexError:
            print('Error: no record found for {0}'.format(data_type))     

      def get_last_N_fitting_data(self, data_type, N):
         try:
            data = self._fitting_data_dico[data_type]
            if len(data) < N:
               return data
            else:
               return self.get_fitting_data(data_type)[len(data)-N:]

         except IndexError:
            print('Warning: only {0} record found for {1}'.format(
                             len(self._fitting_data_dico[self._param.name][data_type]), data_type))
            return self.get_fitting_data(data_type)

      def all_parameters_still_varying(self):
         return not all([pdata.stop_varying for pdata in self.get_param_fitting_data()])

      def parameters_still_varying(self, params):
         return not all([self.get_param_fitting_data_for(pname).stop_varying \
                                                                  for pname in params.param_names])

      def get_varying_parameters(self):
         return [pname for pname in self._params.param_names \
                       if not self.get_param_fitting_data_for(pname).stop_varying]

#      def only_parameters_varying(self, params):
##         other_params = set(self.params) - set(params)
##         return parameters_still_varying(self, params)       
#         pass

      # Return a dictionary of the last fitted data
      def get_last_param_fitting_data(self):
         last_data_dico = {}
         try:
   #         [last_data_dico[p.name] = self._param_data_dico[p.name].get_last_record(p.name) for p in self._params]
            for p in self._params.param_list:   
               last_data_dico[p.name] = self._param_data_dico[p.name].get_last_record(p.name)
         except IndexError:
            print('Error: no record found for parameter')      
         return last_data_dico 

      def has_custom_data(self):
         return len(self._fitting_data_dico['custom_objs']) > 0

      def record_fitting_data(self, index, feval, delta_func, custom_obj=None):
         self._fitting_data_dico['index'].append(index)  
         self._fitting_data_dico['feval'].append(feval)  
         self._fitting_data_dico['delta_func'].append(delta_func)  
         self._fitting_data_dico['delta_func_norm'].append(math.fabs(delta_func))
         if not custom_obj is None:
            self._fitting_data_dico['custom_objs'] = [custom_obj]
#            self._fitting_data_dico['custom_objs'].append(custom_obj)
            #print "***SCDM:  =", len(self._fitting_data_dico['custom_objs'])

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Hold a record of fitting data for a given parameter of class Param
      # *** Notes:
      # - Inner class of Fittingdata 
      # - A private class that should not be accessed directly
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      class _ParamFittingData(object):
         """! 
            Hold a record of fitting data for a given parameter of class Param.
            This is a inner class of class _FittingData
         """
         def __init__(self, param, data_window_size=-1):

            self._param = param                    # target parameter
            self._nb_func_eval = 0                 # nb of objective function evaluations
            self._stop_varying = False             # True if fitting of the parameter must stop
            self._upper_boundary_reached = False   # True if the upper parameter bound was reached
            self._lower_boundary_reached = False   # True if the lower parameter bound was reached
            self._iiter = 0                        # iteration index (incremented for each records)
            self._fitting_data_dico = {param.name:{}}
            self._data_window_size = data_window_size  # size of history data (default is -1: no limit)
            
            # Recording of fitting process
            self._fitting_data_dico[param.name]['steps']  = []     # step history  
            self._fitting_data_dico[param.name]['delta_func'] = [] # algebraic feval difference over a step interval
            self._fitting_data_dico[param.name]['feval']  = []     # function evaluation values
            self._fitting_data_dico[param.name]['deval']  = []     # norm of feval difference over a step interval 
            self._fitting_data_dico[param.name]['params'] = []     # fitted parameter intermediate values
            self._fitting_data_dico[param.name]['dcount'] = []     # count of successive decreasing steps
            self._fitting_data_dico[param.name]['iiter']  = []     # iteration i

         # ~~~~~~~~~~~
         # Properties
         # ~~~~~~~~~~~

         # --- Getters ---

         @property
         def param(self):
            """! Returns target parameter  """
            return self._param

         @property
         def nb_func_eval(self):
            """! Returns number if objective function evaluations  """
            return self._nb_func_eval

         @property
         def upper_boundary_reached(self):
            """! True if upper boundary of parameter range has been reached  """
            return self._upper_boundary_reached

         @property
         def lower_boundary_reached(self):
            """! True if lower boundary of parameter range has been reached  """
            return self._lower_boundary_reached

         @property
         def iiter(self):
            """! Returns current parameter iteration index """
            return self._iiter

         @property
         def data_window_size(self):
            """! Returns size of the data sliding window """
            return self._data_window_size

         @property
         def stop_varying(self):

            """! Returns True one must stop varying the parameter  """
            return self._stop_varying

         # --- Setters ---

         @upper_boundary_reached.setter
         def upper_boundary_reached(self, upper_boundary_reached):
            self._upper_boundary_reached = upper_boundary_reached

         @lower_boundary_reached.setter
         def lower_boundary_reached(self, lower_boundary_reached):
            self._lower_boundary_reached = lower_boundary_reached

         @data_window_size.setter
         def data_window_size(self, data_window_size):
            self._data_window_size = data_window_size

         @stop_varying.setter
         def stop_varying(self, stop_varying):
            self._stop_varying = stop_varying

         # ~~~~~~~~~~~~~~~
         # Public methods
         # ~~~~~~~~~~~~~~~

         def inc_iter(self):
            self._iiter += 1
            return self._iiter

         def inc_nb_func_eval(self):
            self._nb_func_eval += 1

         def clear_data(self):
            """! Clear accumulated data for the parameter """
            self._fitting_data_dico[self._param.name].clear()
            self.__init__(self._param, self._data_window_size)

         def boundary_reached(self):
            return self._upper_boundary_reached or self._lower_boundary_reached 

         # Record one step in fitting process
         def record_fitting_data(self, pval, feval, deval, step, delta_func, dcount, iiter):

            dataDico = self._fitting_data_dico[self._param.name]       
            dataDico['steps'].append(step)     # step values history
            dataDico['delta_func'].append(delta_func)    # calculated feval difference for parameter
            dataDico['feval'].append(feval)    # calculated function values
            dataDico['deval'].append(deval)    # norm of feval difference over a step interval
            dataDico['dcount'].append(dcount)  # count of successive decreasing steps
            dataDico['params'].append(pval)    # calculated parameter      
            dataDico['iiter'].append(iiter)    # iteration i
            if self._data_window_size != -1 and len(dataDico['iiter']) > self._data_window_size:
               dataDico['steps'].pop(0)
               dataDico['delta_func'].pop(0)      
               dataDico['feval'].pop(0)      
               dataDico['deval'].pop(0)      
               dataDico['dcount'].pop(0)      
               dataDico['params'].pop(0)      
               dataDico['iiter'].pop(0)      

         def get_fitting_data_dico(self):
            """! Returns the fitting data dictionary for the target parameter """
            return self._fitting_data_dico[self._param.name]

         def get_fitting_data(self, data_type):
            """! Returns the fitting data dictionary for the target parameter 
                and data type 
            """
            return self._fitting_data_dico[self._param.name][data_type]

         def get_last_fitting_data(self, data_type):
            try:
               return self.get_fitting_data(data_type)[-1]
            except IndexError:
               print('Error: no record found for {0}'.format(data_type))

         def get_last_N_fitting_data(self, data_type, N):
            try:
               data = self._fitting_data_dico[self._param.name][data_type]
               if len(data) < N:
                  return data
               else:
                  return self.get_fitting_data(data_type)[len(data)-N:]
            except IndexError:
               print('Warning: only {0} record found for {1}'.format(len(self._fitting_data_dico[self._param.name][data_type]), data_type))
               return self.get_fitting_data(data_type)

         def get_last_record(self, param_name=None):
            record_dico = {}
            if param_name is None:
               param_name = self._param.name
            for iType in self._fitting_data_dico[param_name].keys():
               record_dico[iType] = self._fitting_data_dico[param_name][iType][-1]         
            record_dico['nfev'] = self._nb_func_eval
            return record_dico    

         def has_records(self):
            return len(self._fitting_data_dico[self._param.name]['params']) > 0


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Hold the fitting results (i.e fitted parameters and fitting data)
   # *** Notes:
   # - Inner class of Minimizer 
   # - A private class that should not be accessed directly
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class _FittingResult(object):
      """! 
         Hold fitting results (i.e fitted parameters and fitting data) 
         This is a inner class of class _FittingData
      """

      def __init__(self, fitting_data):

         self._fitting_data = fitting_data
         self._fitted_param_values = []
         self._nb_iter = 0
         self._errno  = 0    # no error 
         self._status = 0    # status
          
         if self._fitting_data.has_data():
            self._fitting_result_dico = self._fitting_data.get_last_param_fitting_data()
         else:
            self._fitting_result_dico = {}

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~

      # --- Getters ---

      @property
      def fitting_data(self):
         """! Returns the fitting data """
         return self._fitting_data

      @property
      def fitting_result_dico(self):
         """! Returns the fitting result dictionary """
         return self._fitting_result_dico

      @property
      def fitted_param_values(self):
         """! Returns the fitted parameter values """
         return self._fitted_param_values

      @property
      def nb_iter(self):
         """! Returns the final number of iteration spent for fitting """
         return self._nb_iter

      @property
      def errno(self):
         """! Returns the error code after fitting is complete """
         return self._errno

      @property
      def status(self):
         """! Returns the status resulting from fitting """
         return self._status

      # --- Setters ---

      @fitted_param_values.setter
      def fitted_param_values(self, fitted_param_values):
         self._fitted_param_values = fitted_param_values

      @nb_iter.setter
      def nb_iter(self, nb_iter):
         self._nb_iter = nb_iter

      @errno.setter
      def errno(self, errno):
         self._errno = errno

      @status.setter
      def status(self, status):
         self._status = status

      # ~~~~~~~~~~~~~~~
      # Public methods
      # ~~~~~~~~~~~~~~~

      def has_results(self):
         return (len(self._fitting_result_dico) > 0)   

      def get_fitted_value(self, value_type):
         fitted_values = []
         if len(self._fitting_result_dico) > 0:
            for pname in self._fitting_data.param_names:
               if pname in self._fitting_result_dico:
                  fitted_values.append(self._fitting_result_dico[pname][value_type])
            return fitted_values   
         else:
            print('Warning: no result fitted so far')
            return []

      def get_custom_data(self):
         return self._fitting_data.get_fitting_data('custom_objs')

      def get_nb_func_evals(self):
         return sum(self.get_fitted_value('nfev'))


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Hold the status of the fitting operation
   # - Inner class of Minimizer 
   # - A private class that should not be accessed directly
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   class _FittingStatus(object):
      """! 
         Hold the status of the fitting operation.
         This is is a inner class of the Minimizer class.  
      """

      def __init__(self, errno=0, msg='Fitting succeeded', warn_msg=''):
         self._errno = errno    # execution status: 0: success, < 0:failure
         self._msg = msg        # execution status message
         self._warn_msgs = []   # list of warning messages (may be empty)

      # ~~~~~~~~~~~
      # Properties
      # ~~~~~~~~~~~

      # --- Getters ---

      @property
      def errno(self):
         """! Returns the execution error code: 0: success, < 0:failure """
         return self._errno

      @property
      def msg(self):
         """! Returns the execution status message """
         return self._msg

      @property
      def warn_msgs(self):
         """! Returns a list of warning messages (may be empty) """
         return self._warn_msgs

      # --- Setters ---

      @errno.setter
      def errno(self, errno):
         self._errno = errno

      @msg.setter
      def msg(self, msg):
         self._msg = msg

      # ~~~~~~~~~~~~~~~
      # Public methods
      # ~~~~~~~~~~~~~~~

      def add_warn_msg(self, warn_msg):
         """! Add a warning message """
         if not warn_msg in self._warn_msgs:
            self._warn_msgs.append(warn_msg)         

      def set_fitting_status(self, errno, msg):
         self._errno = errno
         self._msg = msg

      def clear(self):
         self._errno = 0       # execution status: 0: success, < 0:failure
         self._msg = 'Fitting succeeded'   # execution status message
         del self._warn_msgs[:] # list of warning messages (may be empty)



# --- EOF scdm.py ---

