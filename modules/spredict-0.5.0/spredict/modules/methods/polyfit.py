"""! 
   spredict - Polynomial Fitting interpolation method - polyfit.py
"""

# --- Python imports
import os, sys
import numpy
import scipy.optimize

# -- External imports

# -- Module-specific imports     
#from sconfig import *
from sp_interp import *    # definition of the BaseInterpolator class


# --------------------------------------------------------------------------------------------------
class Method(BaseMethod):

   """! 
      Inerpolation using Least-Squares Polynomial Fitting (implementation using scipy.leastsq)
      @note extends parent BaseInterpolator class
   """

   def __init__(self, method_config_dir, interp_config_dir):
      """! 
         Construct an Interpolator object 
         @param method_config_dir configuration directory of the method
         @param interp_config_dir directory of extra configuration files for the method
      """      
      BaseMethod.__init__(self, name="polyfit", 
                                method_config_dir=method_config_dir,
                                interp_config_dir=interp_config_dir)

      self.is_local = self.method_config.get_as_boolean("IS_LOCAL", 
                                                        "{0}.METHOD".format(self.name.upper()))

   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

   # ----------------------------------------------------------------------------------------------- 
   def predict(self, model, param_name, input_coords, target_coords, input_values, job, worker):
      """!
         Spatially interpolate input values @c input_values of input coordinates 
         @c input_coords to a target grid @ target_coords by fitting the data to a low-order 
         2D polynomial function with least-squares

         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param input_coords input coordinate tuples (x_input_coords, y_input_coords)
         @param target_coords target coordinate tuples (x_target_coords, y_target_coords)
         @param model parameter values at positions @c input_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
      """

      result_dico = {}

      # --- Method-specific default configuration
      config_dico = self._setup_method_config(model, job, worker)

      # --- Read extra configuration information (from interpolate.cfg) and override previous 
      #     definitions if found     
      self.adjust_method_config(config_dico, model, param_name, job, worker)

      result_dico = self._fit(model, param_name, config_dico, 
                              input_coords, target_coords, input_values, job, worker) 

      return result_dico

   # ~~~~~~~~~~
   # Properties
   # ~~~~~~~~~~
   
   # --- Getters

   # --- Setters



   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   def _fit(self, model, param_name, config_dico, 
                  input_coords, target_coords, input_values, job, worker):


      fitting_result_dico = {}

      # --- Polynomial fitting function
      fitting_funcname = 'polyfit_' + str(config_dico["DEGREE"])
      polyfit_function = eval("self._"+fitting_funcname) 

      # --- Input and Target coordinates
      x_input_coords, y_input_coords = zip(*input_coords)
      x_input_coords, y_input_coords = numpy.asarray(x_input_coords), numpy.asarray(y_input_coords)

      x_target_coords, y_target_coords = zip(*target_coords)
      x_target_coords, y_target_coords = numpy.asarray(x_target_coords), \
                                         numpy.asarray(y_target_coords)

      # --- Input values
      input_values = numpy.asarray(input_values)

      # --- Input arguments for fitting
      #print "###", config_dico, input_values

      guess_values = config_dico['POLY_DEF'][fitting_funcname]
      fitting_args = (x_input_coords, y_input_coords, input_values, polyfit_function)

      #print "### fitting_args:", fitting_args

      # --- Least Squares interpolation
      fitted_coeffs, cov_matrix, infodict, msg, ier = scipy.optimize.leastsq(
                                                      self._least_squares_polyfit, 
                                                      guess_values, args=fitting_args,
                                                      ftol=config_dico["LEASTSQ_FTOL"], 
                                                      xtol=config_dico["LEASTSQ_XTOL"], 
                                                      gtol=config_dico["LEASTSQ_GTOL"], 
                                                      col_deriv=config_dico["LEASTSQ_COL_DERIV"],
                                                      maxfev=config_dico["LEASTSQ_MAX_FEV"], 
                                                      full_output=1) 

      # --- Check results
      fitting_success = ier in (1, 2, 3, 4)
      if fitting_success:
         # --- Fitting successful...

         # --- Evaluate the polynomial values at target coordinates         
         target_values = polyfit_function(x_target_coords, y_target_coords, fitted_coeffs)

         # --- Store the results
         fitting_result_dico["results"] = {} 
         fitting_result_dico["results"]["x_input_coords"]  = x_input_coords
         fitting_result_dico["results"]["y_input_coords"]  = y_input_coords
         fitting_result_dico["results"]["x_target_coords"] = x_target_coords
         fitting_result_dico["results"]["y_target_coords"] = y_target_coords
         fitting_result_dico["results"]["input_values"]  = input_values
         fitting_result_dico["results"]["target_values"] = target_values

         # --- Summary stats
         fitting_result_dico["stats"] = {}
         fitting_result_dico["stats"]["RSS"] = numpy.sum(infodict["fvec"])
 

      else:
         # --- Interpolation error
         err_msg = "{0} - /{1}/image-{2:03d}-{3:1d} - {4} interpolation error on "\
                   "parameter {5}: ier={6}".format(
                    worker.name, job.get_branch_tree(), job.img_no, job.epoch, self.name, 
                    param_name, ier)

         raise SpatialInterpolator.InterpolationError(err_msg)

      return fitting_result_dico


   # -----------------------------------------------------------------------------------------------
   def _setup_method_config(self, model, job, worker):
      """! Read interpolation options from method configuration """   

      section_name = "{0}.{1}".format(self.name.upper(), "OPTIONS")

      config_dico = {}
      config_dico["DEGREE"]   = self.method_config.get_as_int("DEGREE", section_name)
      config_dico["POLY_DEF"] = self.method_config.get_as_dict("POLY_DEF", section_name)

      config_dico["LEASTSQ_MAX_FEV"] = self.method_config.get_as_int("LEASTSQ_MAX_FEV", 
                                                                     section_name)
      config_dico["LEASTSQ_COL_DERIV"] = self.method_config.get_as_int("LEASTSQ_COL_DERIV", 
                                                                        section_name)
      config_dico["LEASTSQ_FTOL"] = self.method_config.get_as_float("LEASTSQ_FTOL", section_name)
      config_dico["LEASTSQ_XTOL"] = self.method_config.get_as_float("LEASTSQ_XTOL", section_name)
      config_dico["LEASTSQ_GTOL"] = self.method_config.get_as_float("LEASTSQ_GTOL", section_name) 

      return config_dico

   # -----------------------------------------------------------------------------------------------
   def _least_squares_polyfit(self, params, x_coords, y_coords, observed_data, fitting_function):
      """! Polynomial fitting with Least Squares """

      #print '***', 'params:', params

      # Estimate semi-variogram values from the model
      estimated_data = fitting_function( x_coords, y_coords, params)
      
   #   print 'observed_data:', observed_data, 'estimated_data:',  estimated_data

      #residuals = numpy.absolute(observed_data - estimated_data)
      residuals = (observed_data - estimated_data)**2

      return residuals

   # -----------------------------------------------------------------------------------------------
   def _polyfit_1(self, x1, y1, coeffs):
      """! Polynomial evaluation function [linear] """

      return coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1

   # -----------------------------------------------------------------------------------------------
   def _polyfit_2(self, x1, y1, coeffs):
      """! Polynomial evaluation function [quadratic] """

      x2 = x1**2
      y2 = y1**2
      x3 = x1**3
      y3 = y1**3

      return coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1 +\
             coeffs[3] * x2 + coeffs[4] * x1*y1 + coeffs[5] * y2

   # -----------------------------------------------------------------------------------------------
   def _polyfit_3(self, x1, y1, coeffs):
      """! Polynomial evaluation function [cubic] """   

      #print x1, y1, coeffs

      x2 = x1**2
      y2 = y1**2
      x3 = x1**3
      y3 = y1**3

      return coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1 +\
             coeffs[3] * x2 + coeffs[4] * x1*y1 + coeffs[5] * y2 +\
             coeffs[6] * x3 + coeffs[7] * x2*y1 + coeffs[8] * x1*y2 + coeffs[9] * y3  

   # -----------------------------------------------------------------------------------------------
   def _polyfit_4(self, x1, y1, coeffs):
      """! Polynomial evaluation function [quadric] """

      x2 = x1**2
      y2 = y1**2
      x3 = x1**3
      y3 = y1**3
      x4 = x1**4
      y4 = y1**4

      return coeffs[0]  + coeffs[1] * x1 + coeffs[2] * y1 +\
             coeffs[3]  * x2 + coeffs[4]  * x1*y1 + coeffs[5]  * y2 +\
             coeffs[6]  * x3 + coeffs[7]  * x2*y1 + coeffs[8]  * x1*y2 + coeffs[9]  * y3 +\
             coeffs[10] * x4 + coeffs[11] * x3*y1 + coeffs[12] * x2*y2 + coeffs[13] * x1*y3 +\
             coeffs[14] * y4 
     
   # -----------------------------------------------------------------------------------------------
   def _polyfit_5(self, x1, y1, coeffs):
      """! Polynomial evaluation function [quintic] """

      x2 = x1**2
      y2 = y1**2
      x3 = x1**3
      y3 = y1**3
      x4 = x1**4
      y4 = y1**4
      x5 = x1**4
      y5 = y1**4

      return coeffs[0] + coeffs[1] * x1 + coeffs[2] * y1 +\
             coeffs[3] * x2 + coeffs[4] * x1*y1 + coeffs[5] * y2 +\
             coeffs[6] * x3 + coeffs[7] * x2*y1 + coeffs[8] * x1*y2 + coeffs[9] * y3 +\
             coeffs[10] * x4 + coeffs[11] * x3*y1 + coeffs[12] * x2*y2 + coeffs[13] * x1*y3 +\
             coeffs[14] * y4 + coeffs[15] * x5 + coeffs[16] * x4*y1 + coeffs[17] * x3*y2 +\
             coeffs[18] * x2*y3 + coeffs[19] * x1*y4 + coeffs[20] * y5


# --- EOF polyfit.py 
