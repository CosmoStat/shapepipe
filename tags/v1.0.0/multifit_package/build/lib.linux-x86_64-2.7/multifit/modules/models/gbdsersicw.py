"""! 
   MultiFit - Bulge-Disk Sersic model using GalSim
"""

# -- Python imports
import math
import numpy
import pyfits
import tempfile

# --- External imports
import galsim
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models'):
      """! 
         Construct a @c gbdsersicw FittingModel object 
         @param config_dir configuration directory of the method
      """  

      BaseFittingModel.__init__(self, name="gbdsersicw", config_dir=config_dir)
      
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

         self.func = self.gbdsersicw

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
      #print "***MFIT: gbdsersicw, getting  out_of_bound_handler..."

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      #print "***MFIT: gbdsersicw, getting  numerical_error_handler..."

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

         model_func = self.gbdsersicw

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
            elif param.name == "dre":   
               param.guess_value = se_galaxy_info_dico["SE_FLUX_RADIUS"][coords]
            elif param.name == "bre":   
               param.guess_value = se_galaxy_info_dico["SE_FLUX_RADIUS"][coords]
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

      #print "***MFIT: in gbdsersicw out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #print "***MFIT: in gbdsersicw numerical_error_handler", iiter, param, value, value_name

      if self.config.has_section("MODEL.GALSIM") and \
         self.config.get_as_boolean("abort_fit_uppon_error", "MODEL.GALSIM"):
         return -1   # always abort fitting
      else:
         return 0    # attempt to fix 


  # -----------------------------------------------------------------------------------------------
   def gbdsersicw(self, model, stamp_size, subscale, I0, dre, bre, df, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """

      #print I0, dn, dre, bre, df, e1, e2, xc, yc, "subscale:",  subscale, "stamp_size:", stamp_size

      try:

         # --- Adjust parameters based on resolution and setup model pixel grid
         if subscale > 1:
            stamp_size *= subscale
            dre *= subscale
            bre *= subscale
            corr = (subscale - 1)/ (2.0 * subscale)    # account for centroid shift
            xc = (xc + corr) * subscale
            yc = (yc + corr) * subscale


         # Galaxy is modelled as a bulge + disk profile
         disk = galsim.Exponential(flux=df, half_light_radius=dre)
         #disk = galsim.Sersic(n=dn, half_light_radius=dre, flux=df)
         disk.applyShear(e1=e1, e2=e2)

         bulge = galsim.DeVaucouleurs(flux=1.0-df, half_light_radius=bre)
         bulge.applyShear(e1=e1, e2=e2)

         # The flux of an Add object is the sum of the component fluxes
         gal = galsim.Add([disk, bulge])

         # Apply shear to both components
         gal.applyShear(e1=e1, e2=e2)

         # The total flux of bulge+disk is set to I0
         gal.setFlux(I0)

         # Shift the object based on the estimated centroid
         half = stamp_size / 2.0
         gal.applyShift(dx=xc-half, dy=yc-half)

         # Render the image
         gal_image = galsim.ImageF(stamp_size, stamp_size, scale=1.0)
         gal.draw(gal_image, dx=1.0)

#         # Render the image
#         pixel_scale = 1.0
#         pix = galsim.Pixel(pixel_scale)
#         gal = galsim.Convolve([pix, gal])

#         gal_image = galsim.ImageF(stamp_size, stamp_size, pixel_scale)
#         gal.draw(gal_image, dx=pixel_scale)

         # Image data as a numpy array
         gal_array = gal_image.array
   
         #gal_image.write("bd.fits")

         # Subscaling if requested
         if subscale > 1:
            gal_array = self._create_low_resolution_stamp(gal_array, stamp_size)

         # Invoke mr_transform
         gal_array = self.mr_transform(model, gal_array)      


      except Exception as detail:

         gal_array = None

         # --- Abort the fit
         print("*** gbdsersicw model *** An error occurred: {0}".format(detail))
         #self.helper.print_error("gbdsersicw model: an error occurred: {0}".format(sys.exc_info()[1]))

      # --- Return model Sersic profile with flux adjusted
      return  gal_array


#   # -----------------------------------------------------------------------------------------------
#   def mr_transform(self, model, input_image):

#      output_image = None
#      mrt_options = model.config.get_as_string("options", "MODEL.MR_TRANSFORM")
#      mrt_path = model.config.get_as_string("exec_path", "MODEL.MR_TRANSFORM")
#      if eval(mrt_options) is None:
#         mrt_options = "-t14 -L"
#      
#      temp = tempfile.NamedTemporaryFile()

#      tmp_dirname, tmp_filename = os.path.split(temp.name)
#      input_filepath  = os.path.join(tmp_dirname, "gal_{0}.fits".format(tmp_filename))
#      output_filepath = os.path.join(tmp_dirname, "gal_{0}.mr".format(tmp_filename))
#      self._write_as_fits(input_image, input_filepath)

#      if len(mrt_options) >  0:
#         cmd_line = "{0} {1} {2} {3}".format(
#            mrt_path, mrt_options, input_filepath, output_filepath)
#      else:
#         cmd_line = "{0} {1} {2}".format(mrt_path, input_filepath, output_filepath)

#      result = os.system(cmd_line)


#      mr_output_filepath = output_filepath

#      if result == 0:
#         output_image = pyfits.getdata(mr_output_filepath).astype(numpy.float32)

#         os.remove(input_filepath)
#         os.remove(mr_output_filepath)

#      else:
#         print("ERROR: mr_transform failed with error {0}".format(result))   

#      return output_image


   # -----------------------------------------------------------------------------------------------
   def _file_exists(self, filepath):
      """! 
         Safe way of checking if a file (not a directory) does exist. 
         @param filepath file path to chech for existence
         @retval True if file exists
         @retval False otherwise
      """
      try:
        with open(filepath) as file:
            return True
      except IOError as e:
        return False

   # -----------------------------------------------------------------------------------------------
   def _write_as_fits(self, data, output_filepath, header=None):
      """! 
         Write a two-dimensional data numpy array as a .fits file to some given path
         @param data data to write
         @param output_filepath full path of the -FITS file
         @param header optional FITS header data 
      """

      if os.path.exists(output_filepath):
         os.remove(output_filepath)

      if data is not None and type(data) == numpy.ndarray:
         pyfits.writeto(output_filepath, data, header=header)       


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
