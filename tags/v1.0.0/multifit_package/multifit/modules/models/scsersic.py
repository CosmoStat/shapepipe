"""! 
   MultiFit - Single-component Sersic model
"""

# -- Python imports
import math
import numpy
import scipy.special

# --- External imports
from multifit import *


# --------------------------------------------------------------------------------------------------
class FittingModel(BaseFittingModel):

   """! 
      Single-Component Sersic galaxy profile 
      @note extends parent BaseFittingModel class
   """

   def __init__(self, config_dir='./config/models', config_filename="scsersic.cfg"):
      """! 
         Construct a @c scsersic FittingModel object 
         @param config_dir configuration directory of the method
         @param config_filename name of the configuration file
      """  

      BaseFittingModel.__init__(self, name="scsersic", 
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

         self.func = self.scsersic


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
      #print "***MFIT: in scsersic, getting out_of_bound_handler()"

      return self._out_of_bound_handler

   @property
   def numerical_error_handler(self):
      """! 
         Return the callback function to notify a numerical error that occurred while computing 
                 the parameter value
         @return the callback function to notify a numerical error
      """
      return self._numerical_error_handler


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

         model_func = self.scsersic

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
                                   galaxy_info_dico, se_galaxy_info_dico):
      """
         Attempt to provide estimates for parameter guess values based on input galaxy data
         @param coords coordinates of left-bottom corner of postage stamp (not in numpy convention)
         @param actual_galaxy_image galaxy postage stamp to fit
         @param se_galaxy_info_dico galaxy data dictionary (depends on dataset)
         @param se_galaxy_info_dico SExtractor data dictionary
         @return list of ModelParam objects with estimated guess values
         @todo possibly estimate (e1,e2) using quadrupole moments
      """   

      for param in self.params:

         if param.guess_value is None:

            if param.name == "flux":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_flux", 
                                                         se_galaxy_info_dico["SE_GAL_FLUX"])

            elif param.name == "dn":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_dn", 1.0)

            elif param.name == "dre":   
               param.guess_value = self.get_guess_value(coords, 
                                                         galaxy_info_dico, se_galaxy_info_dico, 
                                                         param.name, "guess_dre", 
                                                         se_galaxy_info_dico["SE_GAL_FLUX_RADIUS"])
            elif param.name == "xc": 
               if "guess_xc" in self.col_mapping_dico:
                  param.guess_value = self.get_guess_value(coords, 
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            param.name, "guess_xc",
                                                            (actual_galaxy_image.shape[1]-1.0)/2.0)
               else:
                  param.guess_value = (actual_galaxy_image.shape[1]-1.0)/2.0

            elif param.name == "yc": 
               if "guess_yc" in self.col_mapping_dico:
                  param.guess_value = self.get_guess_value(coords, 
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            param.name, "guess_yc", 
                                                            (actual_galaxy_image.shape[0]-1.0)/2.0)
               else:
                  param.guess_value = (actual_galaxy_image.shape[0]-1.0)/2.0

            elif param.name == "e1" or param.name == "e2":
               param.guess_value = 0.0

      return self.params


   #----------------------------------------------------------------------
   def check_fitted_params(self, (x, y), actual_galaxy_image, 
                                         galaxy_info_dico, se_galaxy_info_dico, 
                                         guess_values, param_values):
      """!
         Check that the parameter values  @c param_values are within the allowed bounds
         @param coords coordinates of left-bottom corner of postage stamp (not in numpy convention)
         @param actual_galaxy_image galaxy postage stamp to fit
         @param se_galaxy_info_dico galaxy data dictionary (depends on dataset)
         @param se_galaxy_info_dico SExtractor data dictionaryds
         @return True if the valiity check is successful 
      """

      is_valid = BaseFittingModel.check_fitted_params(self, (x, y), actual_galaxy_image,
                                                            galaxy_info_dico, se_galaxy_info_dico, 
                                                            guess_values, param_values)

      if is_valid:
         # Check global distortion < 1
         e1 = param_values[self.get_param_names().index("e1")]
         e2 = param_values[self.get_param_names().index("e2")]
         is_valid = math.hypot(e1, e2) < 1.0

      #print "***MFIT: is_valid", is_valid, guess_values, param_values
      
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

      #print "***MFIT: in scsersic out_of_bound handler", iiter, param, param_value, bound_index

      status = 0

      return [status, param_value, param.param_range[bound_index]]                              

   # -----------------------------------------------------------------------------------------------
   def numerical_error_handler(self, iiter, param, value, value_name):
         
      #return -1   # always abort fitting
      return 0    # attempt to fix 


   # -----------------------------------------------------------------------------------------------
   def scsersic(self, model, stamp_size, subscale, flux, n, re, e1, e2, xc, yc):
      """!
         Compute the model using specified parameters and arguments.
         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
      """
   
      #print flux, n, re, e1, e2, xc, yc

      scaled_stamp_size = stamp_size * subscale
      lin_space = numpy.linspace(0, stamp_size * subscale, stamp_size * subscale)
      X, Y = numpy.meshgrid(lin_space, lin_space)
      X -= xc * subscale  
      Y -= yc * subscale

      # --- Shear coordinate system itself (inverse magnification transformation)
      bn = self._find_k(n)
      e = math.hypot(e1, e2)
      q = (1.0 - e)/ (1.0 + e)
      psi = 0.5 * math.atan2(e2, e1)

      thetas = numpy.arctan2(Y, X)
      radii  = numpy.hypot(X, Y) * numpy.sqrt(1.0 + (1.0/q**2 - 1.0) * numpy.sin(thetas-psi)**2) 

      scaled_re = re * subscale

      try:

         sersic_flux_array = numpy.expm1(-bn * (radii / scaled_re)**(1.0 / n) ) + 1.0
         #sersic_flux_array = I0 * numpy.exp(-bn * (radii / scaled_re)**(1.0 / n)) 

      except:
         # --- Abort the fit
         #self.helper.print_error("An error occurred: {0}".format(sys.exc_info()[1]))
         print("*** scsersic model *** An error occurred: {0}".format(sys.exc_info()))

         n = 5.0e-1
         sersic_flux_array = numpy.expm1(-bn * (radii/ scaled_re)**(1.0 / n) ) + 1.0

      sersic_flux_array = numpy.reshape(sersic_flux_array, (scaled_stamp_size, scaled_stamp_size) )
      if subscale > 1:
         sersic_flux_array = self._create_low_resolution_stamp(sersic_flux_array, stamp_size)

      # --- Return model Sersic profile with flux adjusted
      return  flux * sersic_flux_array / numpy.sum(sersic_flux_array)  # Normalize on total flux
      


#   # -----------------------------------------------------------------------------------------------
#   def scsersic(self, model, stamp_size, subscale, flux, n, re, e1, e2, xc, yc):
#      """!
#         Compute the model using specified parameters and arguments.
#         In this case, compute the profile @c profile_func on a square grid of size @c grid_size
#      """
#   
#      #print flux, n, re, e1, e2, xc, yc

#      lin_space = numpy.linspace(0, stamp_size-1, int(stamp_size))
#      X, Y = numpy.meshgrid(lin_space, lin_space)
#      X -= xc   
#      Y -= yc 

#      # --- Shear coordinate system itself (inverse magnification transformation)
#      bn = self._find_k(n)

#      [sheared_x_space, sheared_y_space] = self._shear_coordinates(X, Y, e1, e2)
#      distances_from_centre = numpy.ravel(numpy.hypot(sheared_x_space, sheared_y_space))

#      try:
#         sersic_flux_array = numpy.expm1(-bn * (distances_from_centre / re)**(1.0 / n) ) + 1.0

#      except:
#         # --- Abort the fit
#         #self.helper.print_error("An error occurred: {0}".format(sys.exc_info()[1]))
#         print("*** scsersic model *** An error occurred: {0}".format(sys.exc_info()))

#         n = 5.0e-1
#         sersic_flux_array = numpy.expm1(-bn * (distances_from_centre / re)**(1.0 / n) ) + 1.0

#      sersic_flux_array = numpy.reshape(sersic_flux_array, (stamp_size, stamp_size) )
#      if subscale > 1:
#         sersic_flux_array = self._create_low_resolution_stamp(sersic_flux_array, stamp_size)

#      # --- Return model Sersic profile with flux adjusted
#      return  sersic_flux_array * flux / numpy.sum(sersic_flux_array)  # Normalize on total flux

#      return  flux * sersic_flux_array
#      

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # --------------------------------------------------------------------------------------------------
   def _solve_k(self, k, n):

      return scipy.special.gammainc(2*n, k) - 0.5

   # --------------------------------------------------------------------------------------------------
   def _find_k(self, n):

      guess = 2 * n - 1.0/3.0 + 4.0/405.0/n
      true_k = scipy.optimize.fsolve(self._solve_k, guess, args=(n))   
      return true_k[0]

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
