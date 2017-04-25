"""! 
   MultiFit - Moffat model using GalSim
"""

# -- Python imports
import math
import numpy

# --- External imports
import galsim
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models', config_filename="gmoffatsh.cfg"):
      """! 
         Construct a @c gmoffatsh FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="gmoffatsh", 
                                      config_dir=config_dir, config_filename=config_filename)
     
      # --- Set parameter definitions from configuration
      self.set_params(self.get_params_from_config() )

      # --- Set module parametric function
      if self.config.has_section("MODEL.FUNC"):

         # --- Use the profile function specified in the configuration 
         func_name = self.name
         if self.config.has_key("func_name", "MODEL.FUNC"):
            func_name = self.config.get_as_string("func_name", "MODEL.FUNC")
         self.func = eval("self."+func_name)

      else:
         # --- Look for a profile function with same name as the model 

         self.func = self.gmoffatsh


   # ~~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~~   

   # --- Getters     

   @property
   def out_of_bound_handler(self):
      """! 
         Return the callback handler function to notify a parameter value came outside the 
                 allowed range of varation specified by Param.guess_value
         @return the callback handler function
      """
      #print "***MFIT: gmoffatsh, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: gmoffatsh, getting  numerical_error_handler..."

      return self._numerical_error_handler

   # -----------------------------------------------------------------------------------------------
   def get_model_func_from_config(self):
      """!
         Return a pointer to the model function for evaluating the models given parameters
         @return model function
         @note this method is called by multifit.BaseFittingModel and must be implemented
      """
 
      if self.config.has_section("MODEL.FUNC"):

         # --- Use the profile function specified in the configuration 
         model_func_name = self.name

         if self.config.has_key("func_name", "MODEL.FUNC"):
            model_func_name = self.config.get_as_string("func_name", "MODEL.FUNC")
         model_func = eval("self."+model_func_name)

      else:
         # --- Look for a profile function with same name as the model 

         model_func = self.gmoffatsh

      return model_func 


   # -----------------------------------------------------------------------------------------------
   def get_params_from_config(self):
      """! 
         Obtain the list of parameter definitions from the configuration file. 
         In this particular case, the model is generic and parameters must be set programmatically 
      """

      return BaseFittingModel.get_params_from_config(self)


   # -----------------------------------------------------------------------------------------------
   def get_best_guess_params(self, coords, actual_galaxy_image, se_galaxy_info_dico):
      """
         Attempt to provide estimates for parameter guess values based on input galaxy data
         @param coords coordinates (not in numpy convention)
         @param 
         @return list of ModelParam objects with estimated guess values
         @todo possibly estimate (e1,e2) using quadrupole moments
      """   

#      # --- Compute guess (e1, e2)
#      stamp_size = actual_galaxy_image.shape[0]
#      weights = self.compute_gaussian_weights(stamp_size, stamp_size/4.0)
##      weights = numpy.ones((stamp_size, stamp_size))
#      moments = self.compute_moments(actual_galaxy_image, weights)
#      [quad_e1, quad_e2, fit_quad_e] = self.compute_ellipticity(moments)

      for param in self.params:

         if param.guess_value is None:
            if param.name == "I0":   
               param.guess_value = se_galaxy_info_dico["SE_FLUX_BEST"][coords]
            elif param.name == "FWHM":   
               param.guess_value = se_galaxy_info_dico["SE_FWHM_IMAGE"][coords]
            elif param.name == "beta":
               param.guess_value = 3.0
            elif param.name == "xc": 
               param.guess_value = (actual_galaxy_image.shape[1]/2.0)-1.0 
            elif param.name == "yc": 
               param.guess_value = (actual_galaxy_image.shape[0]/2.0)-1.0
            elif param.name == "e1" or param.name == "e2":
               param.guess_value = 0.0
#            elif param.name == "e1":
#               param.guess_value = quad_e1
#            elif param.name == "e2":
#               param.guess_value = quad_e2

      return self.params

      
   # -----------------------------------------------------------------------------------------------
   def check_fitted_params(self, param_values):

      """!
         Check that the parameter values  @c param_values are within the allowed bounds
         @param param_values values of the parameters to check against the bounds
      """

      is_valid = BaseFittingModel.check_fitted_params(self, param_values)
      if is_valid:
         # Check global distortion < 1
         e1 = param_values[self.get_param_names().index("e1")]
         e2 = param_values[self.get_param_names().index("e2")]
         is_valid = (math.hypot(e1, e2) < 1.0)
   
         #print "***MFIT: is_valid", is_valid, e1, e2, math.hypot(e1, e2)
      
      return is_valid

   # -----------------------------------------------------------------------------------------------
   def out_of_bound_handler(self, iiter, param, param_value, bound_index):
         
   #   if param.name == 'e1':
   #      status = 0  # continue fitting => may reach < +-1 later on
   #
   #   elif param.name == 'e2':
   #      status = 0  # continue fitting => may reach < +-1 later on

   #   elif param.name == 'dre' and bound_index == 0:
   #      status = -1
   #   elif param.name == 'dre' and bound_index == 1:
   #      status = -1
   #   elif param.name == 'bre' and bound_index == 0:
   #      status = -1
   #   elif param.name == 'bre' and bound_index == 1:
   #      status = -1
   #   else:
   #      status = 0  # continue fitting

      #print "***MFIT: in gmoffatsh out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      if self.config.has_section("MODEL.GALSIM") and \
         self.config.get_as_boolean("abort_fit_uppon_error", "MODEL.GALSIM"):
         return -1   # always abort fitting
      else:
         return 0    # attempt to fix 



  # -----------------------------------------------------------------------------------------------
   def gmoffatsh(self, model, stamp_size, subscale, I0, beta, FWHM, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print I0, beta, FWHM, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # --- Adjust parameters based on resolution and setup model pixel grid
         if subscale > 1:
            stamp_size *= subscale
            FWHM *= subscale
            corr = (subscale - 1)/ (2.0 * subscale)    # account for centroid shift
            xc = (xc + corr) * subscale
            yc = (yc + corr) * subscale

         #print "***gmoffatsh:",  stamp_size, subscale, I0, beta, FWHM, e1, e2, xc, yc

         # Accuracy for real-space approximations
         gsparams = galsim.GSParams(xvalue_accuracy=1.0e-12)

         # Star modeled using a Moffat profile 
         star = galsim.Moffat(beta=beta, fwhm=FWHM, flux=I0, trunc=0, gsparams=gsparams)    

         # Apply shear to both components
         star.applyShear(g1=e1, g2=e2)

         # The center of the object is normally placed at the center of the postage stamp image.
         # You can change that with applyShift:
         half = stamp_size / 2.0
         star.applyShift(dx=xc-half, dy=yc-half)

         # Render the image
         star_image = galsim.ImageF(stamp_size, stamp_size, scale=1.0)
         star.draw(star_image, dx=1.0)

         # Image data
         star_array = star_image.array
   
         if subscale > 1:
            star_array = self._create_low_resolution_stamp(star_array, stamp_size/subscale)

      except Exception as detail:

         # --- Abort the fit
         star_array = numpy.nan

         print("*** gmoffatsh model *** An error occurred: {0}".format(detail))
         #self.helper.print_error("gmoffatsh model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  star_array

      
   # -----------------------------------------------------------------------------------------------
   def _create_low_resolution_stamp(self, highres_stamp, new_smaller_size):

      """!
         Build low resolution stamp, with a size multiple of the original stamp size
         (no interpolation, just binning). Total flux is conserved.
      """

      # --- Rescale stamp values
      scale_factor  = math.floor(highres_stamp.shape[0] / new_smaller_size)** 2

      # --- Build low resolution stamp (no interpolation, just binning)
      return self._rebin_reduce(highres_stamp, new_smaller_size, new_smaller_size) * scale_factor


   # -----------------------------------------------------------------------------------------------
   def _rebin_reduce(self, stamp, *args):
      """!
         Rebin a stamp with a smaller size, multiple of the original size
         [ same syntax as the IDL equivalent: bin(array/matrix, size, size, sample=0) ]
      """

      shape = stamp.shape
      lenShape = len(shape)
      factor = numpy.asarray(shape)/numpy.asarray(args)
      evList = ['stamp.reshape('] + \
                ['args[%d],factor[%d],'%(i,i) for i in xrange(lenShape)] + \
                [')'] + ['.mean(%d)'%(i+1) for i in xrange(lenShape)]

      return eval(''.join(evList))

### EOF gmoffatsh.py

