"""! 
   MultiFit - Moffat model
"""

# -- Python imports
import math
import numpy

# --- External imports
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models', config_filename="moffat.cfg"):
      """! 
         Construct a @c moffat FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="moffat", 
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
         self.func = self.moffat

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
      #print "***MFIT: moffat, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: moffat, getting  numerical_error_handler..."

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

         model_func = self.moffat

      return model_func 


   # -----------------------------------------------------------------------------------------------
   def get_params_from_config(self):
      """! 
         Obtain the list of parameter definitions from the configuration file. 
         In this particular case, the model is generic and parameters must be set programmatically 
      """

      return BaseFittingModel.get_params_from_config(self)


   # -----------------------------------------------------------------------------------------------
   def get_best_guess_params(self, coords, actual_galaxy_image, 
                                          se_galaxy_info_dico):
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

      #print "***MFIT: in moffat out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #return -1   # always abort fitting

      #print "***MFIT: in moffat numerical_error_handler", iiter, param, value, value_name

      return 0    # attempt to fix 


  # -----------------------------------------------------------------------------------------------
   def moffat(self, model, stamp_size, subscale, I0, beta, FWHM, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print I0, beta, FWHM, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # The alpha parameter must be scaled according to <subscale>
         scaled_stamp_size = stamp_size * subscale
         scaled_alpha = 0.5 * FWHM * subscale / math.sqrt(2** (1/beta) - 1)
         corr = (subscale - 1) / (2.0 * subscale)    # account for centroid shift for subscale > 1
         scaled_xc = (xc + corr) * subscale
         scaled_yc = (yc + corr) * subscale

         # Make a coordinate system of size <lin_space_size> with origin at (<xc>, <yc>)
         lin_space = numpy.linspace(0, scaled_stamp_size-1, int(scaled_stamp_size))
         [x_space, y_space] = numpy.meshgrid(lin_space, lin_space)
         [x_space, y_space] = (x_space - scaled_xc), (y_space - scaled_yc)

         # Shear coordinate system itself (inverse magnification transformation)
         [sheared_x_space, sheared_y_space] = self._shear_coordinates(x_space, y_space, e1, e2)
         sheared_distances_from_centre = numpy.ravel(numpy.hypot(sheared_x_space, sheared_y_space))

         # Flux at distance R from centre
         moffat_flux_array = ((1.0 + (sheared_distances_from_centre / scaled_alpha) ** 2)** (-beta))   
         moffat_flux_sum = numpy.sum(moffat_flux_array) / subscale**2

         # Reduced flux at distance R from centre
         psf_reduced_flux_array = moffat_flux_array / moffat_flux_sum   

         # Total PSF flux
         psf_total_reduced_flux = sum(psf_reduced_flux_array)
         psf_reduced_flux_array /= psf_total_reduced_flux    # Normalise flux in stamps
         psf_reduced_flux_array *= I0

         # Build PSF (normalised) stamp
         psf_stamp = numpy.reshape(psf_reduced_flux_array, (scaled_stamp_size, scaled_stamp_size))

         # Go back to lower resolution
         if subscale > 1:
            psf_stamp = self._create_low_resolution_stamp(psf_stamp, stamp_size)

      except Exception as detail:

         # --- Abort the fit

         psf_stamp = numpy.nan

         print("*** moffat model *** An error occurred: {0}".format(detail))
         #self.helper.print_error("moffat model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  psf_stamp

      
   # -----------------------------------------------------------------------------------------------
   def _shear_coordinates(self, x_space, y_space, g1, g2):

      xform_mat = numpy.matrix([ [1 - g1, -g2], 
                                 [-g2, 1 + g1] ])  # transform matrix

      return numpy.dot(xform_mat, numpy.asmatrix([numpy.ravel(x_space), numpy.ravel(y_space)]))  


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

