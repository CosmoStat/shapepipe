"""! 
   spredict - target values set as function of input values - function.py
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
      BaseMethod.__init__(self, name="function", 
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

      result_dico = self._set_values(model, param_name, config_dico, 
                              input_coords, target_coords, input_values, job, worker) 

      return result_dico


   # ~~~~~~~~~~~~~~~
   # Private Methods
   # ~~~~~~~~~~~~~~~

   def _set_values(self, model, param_name, config_dico, 
                         input_coords, target_coords, input_values, job, worker):

      """ Set target pameter values by applying some function of input parameters """

      fitting_result_dico = {}

      # --- Expression for calculating target values based on input values
      try:


         # --- Input and Target coordinates
         x_input_coords, y_input_coords = zip(*input_coords)
         x_input_coords, y_input_coords = numpy.asarray(x_input_coords), numpy.asarray(y_input_coords)

         x_target_coords, y_target_coords = zip(*target_coords)
         x_target_coords, y_target_coords = numpy.asarray(x_target_coords), \
                                         numpy.asarray(y_target_coords)

         # --- Target values
         target_values = numpy.repeat(eval(config_dico["PREDICTED_VALUES"]), len(target_coords))

         # --- Store the results
         fitting_result_dico["results"] = {} 
         fitting_result_dico["results"]["x_input_coords"]  = x_input_coords
         fitting_result_dico["results"]["y_input_coords"]  = y_input_coords
         fitting_result_dico["results"]["x_target_coords"] = x_target_coords
         fitting_result_dico["results"]["y_target_coords"] = y_target_coords
         fitting_result_dico["results"]["input_values"]  = input_values
         fitting_result_dico["results"]["target_values"] = target_values

      except Exception as detail:

         # --- Interpolation error
         err_msg = "{0} - /{1}/image-{2:03d}-{3:1d} - {4} interpolation error on "\
                   "parameter {5} - {6}".format(
                    worker.name, job.get_branch_tree(), job.img_no, job.epoch, self.name, 
                    param_name, detail)
         raise SpatialInterpolator.InterpolationError(err_msg)
 
      return fitting_result_dico


   # -----------------------------------------------------------------------------------------------
   def _setup_method_config(self, model, job, worker):
      """! Read interpolation options from method configuration """   

      section_name = "{0}.{1}".format(self.name.upper(), "OPTIONS")

      config_dico = {}
      config_dico["PREDICTED_VALUES"]  = self.method_config.get_as_string("PREDICTED_VALUES", 
                                                                          section_name)

      return config_dico


# --- EOF function.py 
