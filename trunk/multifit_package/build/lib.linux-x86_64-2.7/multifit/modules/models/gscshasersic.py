"""! 
   MultiFit - Disk Sersic model using GalSim
"""

# -- Python imports
import math
import numpy
import scipy.special

# --- External imports
import galsim
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models', config_filename="gschasersic.cfg"):
      """! 
         Construct a @c gscshasersic FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="gscshasersic", 
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

         self.func = self.gscshasersic

      if self.config.has_section("MODEL.GALSIM"):

         # If True, abort the fit when an error occurs; otherwise tries to find another best fit
         self._abort_fit_uppon_error = self.config.get_as_boolean("abort_fit_uppon_error", 
                                                                  "MODEL.GALSIM")

         # Smallest allowed bulge effective radius; beyond, assumes single disk component galaxy
         self._min_bulge_radius = self.config.get_as_float("min_bulge_sersic_radius", 
                                                           "MODEL.GALSIM")

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
      #print "***MFIT: gscshasersic, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: gscshasersic, getting  numerical_error_handler..."

      return self._numerical_error_handler

   @property
   def abort_fit_uppon_error(self):
      """! 
         Tells whether a fit should be aborted when an error occurs
         @return whethed a fit should be aborted when an error occurs
      """
      return self._abort_fit_uppon_error

   @property
   def min_bulge_radius(self):
      """! 
         Return the smallest allowed galaxy bulge effective radius
         @return the smallest allowed galaxy bulge effective radius
      """
      return self._min_bulge_radius


   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

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

         model_func = self.gscshasersic

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
#               dn  = self.config.get_as_float("guess", "MODEL.PARAMS.dn")
#               dre = se_galaxy_info_dico["SE_FLUX_RADIUS"][coords]
#               F = se_galaxy_info_dico["SE_FLUX_BEST"][coords]
#               gal = galsim.Sersic(n=dn, half_light_radius=dre, flux=F)
#               stamp_size = actual_galaxy_image.shape[0]
#               gal_image = galsim.ImageF(stamp_size, stamp_size, scale=1.0)
#               gal.draw(gal_image, dx=1.0)
#               param.guess_value = numpy.max(gal_image.array)  
#       
#               print "IO:", param.guess_value, "F:", F, "dn:", dn, "dre:", dre

#               xc = yc = (actual_galaxy_image.shape[0]-1)/2.0 
##               guess_I0 = numpy.max(actual_galaxy_image)
#               guess_I0 = numpy.max( actual_galaxy_image[xc-5+1:xc+5+1, yc-5+1:yc+5+1] )
#               print "***", (xc, yc), guess_I0, numpy.sum(actual_galaxy_image)
#               param.guess_value = guess_I0

#               #print "*** Begin Reading guess I0..."
               dn  = self.config.get_as_float("guess", "MODEL.PARAMS.dn")
               dre = se_galaxy_info_dico["SE_FLUX_RADIUS"][coords]
               F = se_galaxy_info_dico["SE_FLUX_BEST"][coords]
               bn = 1.9992 * dn - 0.3271
               param.guess_value = 2.0 * (bn**dn / dre)**2 * F / (2* math.pi * scipy.special.gamma(2 * dn)) 
               #print "*** End Reading guess I0", "dn:", dn, "dre:", dre, "F:", F, "bn:", bn, "Guess I0:", guess_I0, scipy.special.gamma(2 * dn)
            elif param.name == "dre":   
               param.guess_value = se_galaxy_info_dico["SE_FLUX_RADIUS"][coords]
            elif param.name == "xc": 
               param.guess_value = (actual_galaxy_image.shape[1]-1.0)/2.0
            elif param.name == "yc": 
               param.guess_value = (actual_galaxy_image.shape[0]-1.0)/2.0
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

      #print "***MFIT: in gscshasersic out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #print "***MFIT: in gbdsersic numerical_error_handler", iiter, param, value, value_name

      if self.abort_fit_uppon_error:
         return -1   # always abort fitting
      else:
         return 0    # attempt to fix 

  # -----------------------------------------------------------------------------------------------
   def gscshasersic(self, model, stamp_size, subscale, I0, dn, dre, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print I0, dn, dre, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # Galaxy is modelled as a bulge + disk profile
         pixel_scale = 1.0
         gal = galsim.Sersic(n=dn, half_light_radius=dre, flux=pixel_scale)

         gal_image = galsim.ImageF(stamp_size, stamp_size, scale=pixel_scale)
         gal.draw(gal_image, dx=pixel_scale)

         #print "init mode:", numpy.max(gal_image.array), "init flux:", numpy.sum(gal_image.array),
 
         # --- Normalize on profile peak value
         gal = gal / numpy.max(gal_image.array)          

         # Shear transform
         gal = gal.shear(g1=e1, g2=e2)

         # Shift the object based on the estimated centroid
         half = stamp_size / 2.0
         gal = gal.shift(dx=xc-half, dy=yc-half)

         # Redraw
         gal.draw(gal_image, dx=pixel_scale)

         # Image data as a numpy array
         gal_array = gal_image.array * I0
   
         #print "I0:", I0, "new mode:", numpy.max(gal_array), "new Tot flux:", numpy.sum(gal_array)

         #gal_image.write("bd_{0}_{1}_{2}_{3}_{4}_{5}_{6}.fits".format(I0, dn, dre, e1, e2, xc, yc))

      except Exception as detail:

         gal_array = None

         # --- Abort the fit
         print("*** gscshasersic model *** An error occurred: {0}".format(detail))
         #self.helper.print_error("gscshasersic model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  gal_array

      
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
