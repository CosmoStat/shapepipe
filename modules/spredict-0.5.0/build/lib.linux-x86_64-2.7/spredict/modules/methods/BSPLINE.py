"""! 
   spredict - Polynomial Fitting interpolation method - BSPLINE.py
"""

# --- Python imports
import os, sys
import numpy
import scipy.interpolate

# -- External imports

# -- Module-specific imports     
from sp_interp import *    # definition of the BaseInterpolator class


# --------------------------------------------------------------------------------------------------
class Method(BaseMethod):

   """! 
      Inerpolation using BSPLINE (smoothing bivariate splines, from scipy.interpolate)
      @note extends parent BaseInterpolator class
   """

   def __init__(self, method_config_dir, interp_config_dir):
      """! 
         Construct an Interpolator object 
         @param method_config_dir configuration directory of the method
         @param interp_config_dir directory of extra configuration files for the method
      """      
      BaseMethod.__init__(self, name="BSPLINE", 
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
         @c input_coords to a target grid @ target_coords using smmoothing bivariate splines

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

      result_dico = self._interpolate(model, param_name, config_dico, 
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

   def _interpolate(self, model, param_name, config_dico, 
                    input_coords, target_coords, input_values, job, worker):


      result_dico = {}

      # --- Input and Target coordinates
      x_input_coords, y_input_coords = zip(*input_coords)
      x_input_coords, y_input_coords = numpy.asarray(x_input_coords), numpy.asarray(y_input_coords)

      x_target_coords, y_target_coords = zip(*target_coords)
      x_target_coords, y_target_coords = numpy.asarray(x_target_coords), \
                                         numpy.asarray(y_target_coords)

      # --- Input values
      input_values = numpy.asarray(input_values)

      if worker.logging_enabled():
         worker.logger.log_info(
           "{0} - {1} - Building Bivariate Splines for parameter: {2} ({3})...".format(
                                                   worker.name, self.name, param_name, model.name))

      # --- Construct the splines based on input values
      bispl = scipy.interpolate.SmoothBivariateSpline(
                                     x_input_coords, y_input_coords, input_values, 
                                     w=None, bbox=[None, None, None, None], 
                                     kx=config_dico["DEGREE"], ky=config_dico["DEGREE"], 
                                     s=config_dico["SMOOTH"], eps=None)

      if worker.logging_enabled():
         worker.logger.log_info(
           "{0} - {1} - Evaluating Bivariate Splines for parameter: {2} ({3}) "\
           "[deg={4} smooth={5}]...".format(worker.name, self.name, param_name, model.name, 
                                            config_dico["DEGREE"], config_dico["SMOOTH"]))

      # --- Compute the spline value at target coordinates
      target_values = bispl.ev(x_target_coords, y_target_coords)

      # --- Store the results
      result_dico["results"] = {} 
      result_dico["results"]["x_input_coords"]  = x_input_coords
      result_dico["results"]["y_input_coords"]  = y_input_coords
      result_dico["results"]["x_target_coords"] = x_target_coords
      result_dico["results"]["y_target_coords"] = y_target_coords
      result_dico["results"]["input_values"]  = input_values
      result_dico["results"]["target_values"] = target_values

      # --- Summary stats
      result_dico["stats"] = {}
      result_dico["stats"]["RSS"] = bispl.get_residual()

      return result_dico


   # -----------------------------------------------------------------------------------------------
   def _setup_method_config(self, model, job, worker):
      """! Read interpolation options from method configuration """   

      section_name = "{0}.{1}".format(self.name.upper(), "OPTIONS")

      config_dico = {}
      config_dico["DEGREE"]   = self.method_config.get_as_int("DEGREE", section_name)

      for var in ["XB", "XE", "YB", "YE", "SMOOTH"]:
         config_dico[var] = self.method_config.get_as_string(var, section_name)
         if eval(config_dico[var]) is not None:
             config_dico[var] = self.method_config.get_as_float(var, section_name)    

      config_dico["FULL_OUTPUT"] = self.method_config.get_as_int("FULL_OUTPUT", section_name)

      return config_dico


# --- EOF BSPLINE.py 
