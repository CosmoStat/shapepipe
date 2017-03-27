"""! 
   @package mks.mks_sim Simulation Management
   @author Marc Gentile
   @file mks_sim.py
   Simulation Management
""" 

# -- Python imports
import os, sys
import math
import time
import re
import itertools
import numpy
import scipy.signal
import scipy.special

# -- External imports
from slogger import FileLogger   # logging
import galsim                    # galsim

# --- Module-specific imports
from mks_help import *           # helper utility functions


# -------------------------------------------------------------------------------------------------
class MksSimulator(object):
   
   """! Simulator class """

   # -----------------------------------------------------------------------------------------------
   def __init__(self, master):
      """
         ! Simulator class constructor
         @param master master process instance
      """

      self._start_time = time.clock()     # start time
      self._helper = MksHelper()          # helper utility functions
      #self._sequence_dico = self._setup_sequence_config_dico(master)   # sequence cycle management

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~
   
   @property
   def helper(self):
      """! @return the MksHelper instance. """
      return self._helper
      
#    @property
#    def sequence_dico(self):
#       """! @return the dictionary of sequence cycle objects. """
#       return self._sequence_dico
      
   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def run_simulations(self, file_type, job, worker):

      # --- Setup mosaic images
      gal_image = self._setup_galaxy_mosaic(file_type, job, worker)
      psf_image = self._setup_psf_mosaic(file_type, job, worker) 

      # --- Populate mosaic images
      success = self._populate_mosaics(gal_image, psf_image, file_type, job, worker)

      return success
   

   # -----------------------------------------------------------------------------------------------
   def setup_sequence_config_dico(self, master):
      
      sequence_dico = {}

      # --- Handle potential sequence info in "GALAXY.MODEL" section 
      model_name = master.config.get_as_string("MODEL_NAME", "GALAXY.MODEL")
      section_name = "GALAXY.MODEL.{0}".format(model_name.upper())
      section_dico = master.config.get_section_data(section_name)
      for section_key in section_dico.keys():
         dist_info_dico = eval(section_dico[section_key])
         if dist_info_dico["dist"].lower() == "sequence": 
            sequence_dico[(section_name, section_key)] = itertools.cycle(dist_info_dico["params"])

      # --- Handle potential sequence info in "PSF.MODEL" section 
      model_name = master.config.get_as_string("MODEL_NAME", "PSF.MODEL")
      section_name = "PSF.MODEL.{0}".format(model_name.upper())
      section_dico = master.config.get_section_data(section_name)
      for section_key in section_dico.keys():
         dist_info_dico = eval(section_dico[section_key])
         if dist_info_dico["dist"].lower() == "sequence": 
            sequence_dico[(section_name, section_key)] = itertools.cycle(dist_info_dico["params"])

      # --- Handle potential sequence info in "NOISE.MODEL" section
      model_name = master.config.get_as_string("GALAXY_NOISE_MODEL_NAME", "NOISE.MODEL")
      if model_name != "no_noise":
         section_name = "NOISE.MODEL.{0}".format(model_name.upper())
         section_dico = master.config.get_section_data(section_name)
         for section_key in section_dico.keys():
            dist_info_dico = eval(section_dico[section_key])
            if dist_info_dico["dist"].lower() == "sequence": 
               sequence_dico[(section_name, section_key)] = itertools.cycle(
                                                                        dist_info_dico["params"])

      # --- Handle potential sequence info in "SHEAR" section 
      section_name = "SHEAR"
      section_dico = master.config.get_section_data(section_name)
      for section_key in section_dico.keys():
         dist_info_dico = eval(section_dico[section_key])
         if dist_info_dico["dist"].lower() == "sequence": 
            sequence_dico[(section_name, section_key)] = itertools.cycle(dist_info_dico["params"])

      return sequence_dico
            
#    # -----------------------------------------------------------------------------------------------
#    def setup_power_spectrum_dico(self, stamp_size, pixel_scale, nb_stamps, random_seed,  master):
#       
#       ps_dico = {}
# 
#       # --- Handle potential sequence info in "GALAXY.MODEL" section 
#       model_name = master.config.get_as_string("MODEL_NAME", "GALAXY.MODEL")
#       section_name = "GALAXY.MODEL.{0}".format(model_name.upper())
#       section_dico = master.config.get_section_data(section_name)
#       for section_key in section_dico.keys():
#          dist_info_dico = eval(section_dico[section_key])
#          if dist_info_dico["dist"].lower() == "ps":
#             # power_function which defines the E-mode or B-mode power to use
#             rng = galsim.BaseDeviate(random_seed+nb_stamps)
#             exp = section_dico["params"][0]  
#             ps = galsim.PowerSpectrum(lambda k : k**exp)
#             _, _ = ps.buildGrid(grid_spacing=stamp_size * pixel_scale, ngrid=nb_stamps,
#                                 rng=rng.duplicate())
#             
#             ps_dico[(section_name, section_key)] = galsim.PowerSpectrum(lambda k : k**exp)
            
   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~
   
   # -----------------------------------------------------------------------------------------------
   def _setup_galaxy_mosaic(self, file_type, job, worker):

      """! Create an empty mosaic of galaxies """

      # Compute the number of objects to create
      nb_gal_angles  = worker.config.get_as_int("NB_ANGLES_RING", "SIMULATION")

      # Compute mosaic size
      gal_stamp_size  = worker.config.get_as_int(
                                             "GAL_STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      pixel_scale = worker.config.get_as_float("PIXEL_SCALE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      nb_realizations = worker.config.get_as_int("NB_REALIZATIONS", "SIMULATION")
      nb_objects = nb_gal_angles * nb_realizations

      nb_rows = math.floor(math.sqrt(nb_objects))
      nb_cols = math.ceil(nb_objects / nb_rows)

      mosaic_width  = int(nb_cols * gal_stamp_size)
      mosaic_height = int(nb_rows * gal_stamp_size)  
      gal_image = galsim.ImageF(mosaic_width, mosaic_height, scale=pixel_scale)

      return gal_image

   # -----------------------------------------------------------------------------------------------
   def _setup_psf_mosaic(self, file_type, job, worker):

      """! Create an empty mosaic of PSFs """

      # Compute the number of objects to create
      nb_gal_angles  = worker.config.get_as_int("NB_ANGLES_RING", "SIMULATION")
         
      # --- There should be one PSF per galaxy, so that they can be convolved with each other
      psf_stamp_size  = worker.config.get_as_int(
                                            "PSF_STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      pixel_scale = worker.config.get_as_float("PIXEL_SCALE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      nb_realizations = worker.config.get_as_int("NB_REALIZATIONS", "SIMULATION")
      nb_objects = nb_gal_angles * nb_realizations

      nb_rows = math.floor(math.sqrt(nb_objects))
      nb_cols = math.ceil(nb_objects / nb_rows)
      mosaic_width  = int(nb_cols * psf_stamp_size)
      mosaic_height = int(nb_rows * psf_stamp_size)  

      psf_image = galsim.ImageF(mosaic_width, mosaic_height, scale=pixel_scale)

      return psf_image
   
   # -----------------------------------------------------------------------------------------------
   def _populate_mosaics(self, gal_image, psf_image, file_type, job, worker):

      success = True

      # --- Image properties
      pixel_scale = worker.config.get_as_float("PIXEL_SCALE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      gal_stamp_size  = worker.config.get_as_int(
                                           "GAL_STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
      psf_stamp_size  = worker.config.get_as_int(
                                           "PSF_STAMP_SIZE", "PRIMARY_DATASET.IMAGE_PROPERTIES")
#       pixel_float_size = worker.config.get_as_string("PIXEL_FLOAT_SIZE", 
#                                                      "PRIMARY_DATASET.IMAGE_PROPERTIES")
      gal_image_data = gal_image.array

      nb_x_stamps, nb_y_stamps = gal_image_data.shape[1] / gal_stamp_size, \
                                 gal_image_data.shape[0] / gal_stamp_size

      # --- Output file paths
      file_paths = job.get_file_paths()
      gal_image_filepaths   = job.get_galaxy_image_file_paths(worker)
      psf_image_filepaths   = job.get_star_image_file_paths(worker)
      gal_catalog_filepaths = job.get_galaxy_catalog_file_paths(worker)
      psf_catalog_filepaths = job.get_star_catalog_file_paths(worker)
      
      if len(gal_image_filepaths)   == 0 or len(psf_image_filepaths)   == 0 or \
         len(gal_catalog_filepaths) == 0 or len(psf_catalog_filepaths) == 0:
         if worker.logging_enabled():
            worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Invalid or missing output image or catalog path - Check paths in: {3}".format(
                                               worker.name, job.img_no, job.epoch, file_paths))
         return False

      gal_output_image_filepath   = gal_image_filepaths[0]
      psf_output_image_filepath   = psf_image_filepaths[0]
      gal_output_catalog_filepath = gal_catalog_filepaths[0]
      psf_output_catalog_filepath = psf_catalog_filepaths[0]

      # --- Galsim params
      gsparams = galsim.GSParams(xvalue_accuracy=1.0e-12, maximum_fft_size = 32768)

      # --- Read simulation data
      dist_scope = "field"   # distribution values are calculated once for the whole field
      sim_data_dico = self._read_sim_model_data(job, worker)

      # --- Read model data
      galaxy_model_name = worker.config.get_as_string("MODEL_NAME", "GALAXY.MODEL")
      psf_model_name    = worker.config.get_as_string("MODEL_NAME", "PSF.MODEL")

      galaxy_noise_model_name = worker.config.get_as_string("GALAXY_NOISE_MODEL_NAME", "NOISE.MODEL")

      # --- Check if the galaxy model is implemented
      if not  galaxy_model_name in ["gscsersic", "gbdsersic_r", "gaussian"]:     
         if worker.logging_enabled():
            worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Galaxy model {3} is currently not available".format(
                                    worker.name, job.img_no, job.epoch, galaxy_model_name))                  
         return False   

      # --- Check if the PSF model is implemented
      if not  psf_model_name in ["gaussian", "gmoffat"]:     
         if worker.logging_enabled():
            worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "PSF model {3} is currently not available".format(
                                    worker.name, job.img_no, job.epoch, psf_model_name))                  
         return False  

      # --- Compute model data from configuration
      galaxy_model_dico        = self._read_galaxy_model_data(galaxy_model_name, job, worker)
      psf_model_dico           = self._read_psf_model_data(psf_model_name, job, worker)
      shear_model_dico         = self._read_shear_model_data(job, worker)
      galaxy_noise_model_dico  = self._read_noise_model_data(galaxy_noise_model_name, job, worker)

      #print "*** galaxy_noise_model_dico:", galaxy_noise_model_dico

      # --- Personalize mosaic file names
      gal_output_image_filepath, gal_output_catalog_filepath = \
                             self._adjust_galaxy_filename(
                                   galaxy_model_name, psf_model_name, galaxy_noise_model_name,  
                                   gal_stamp_size,  
                                   galaxy_model_dico, psf_model_dico, 
                                   galaxy_noise_model_dico, shear_model_dico,
                                   gal_output_image_filepath, gal_output_catalog_filepath,
                                   job, worker)
                             
      psf_output_image_filepath, psf_output_catalog_filepath = \
                             self._adjust_psf_filename(
                                   psf_model_name, psf_stamp_size,
                                   psf_model_dico,
                                   psf_output_image_filepath, psf_output_catalog_filepath,
                                   job, worker)

      # --- Simulation properties
      nb_gal_angles   = sim_data_dico["NB_ANGLES_RING"]
      nb_realizations = sim_data_dico["NB_REALIZATIONS"]
      apply_PSF       = sim_data_dico["APPLY_PSF"]
      convol_method   = sim_data_dico["CONVOLUTION_METHOD"]
      nb_objects = nb_gal_angles * nb_realizations       

      # --- Create data dictionaries that will hold the data
      gal_cat_data_dico = self._create_gal_cat_data_dico(galaxy_model_name, job, worker)
      psf_cat_data_dico = self._create_psf_cat_data_dico(psf_model_name, job, worker)

#        print "GAL_KEYS:", gal_cat_data_dico.keys() 

      # --- Populate galaxy and PSF mosaics...
      iobj  = 0   # object index
      irot  = 0   # ring index 
      ireal = -1  # realization index
      random_seed = 3339201
      dist_scope = "object"   # distribution values are recomputed for each object if they have
                              # a scope of 'object'

      for iy in range(nb_y_stamps):

         gal_yc = iy * gal_stamp_size + (gal_stamp_size - 1.0) / 2.0
         psf_yc = iy * psf_stamp_size + (psf_stamp_size - 1.0) / 2.0

         for ix in range(nb_x_stamps):

            # --- Galaxy (geometric) center coordinates   
            gal_xc = ix * gal_stamp_size + (gal_stamp_size - 1.0) / 2.0
            psf_xc = ix * psf_stamp_size + (psf_stamp_size - 1.0) / 2.0


            # --- Local coordinates
            gal_x = math.floor(gal_xc / gal_stamp_size) * gal_stamp_size
            gal_y = math.floor(gal_yc / gal_stamp_size) * gal_stamp_size

            psf_x = math.floor(psf_xc / psf_stamp_size) * psf_stamp_size
            psf_y = math.floor(psf_yc / psf_stamp_size) * psf_stamp_size

            # --- Random number generator
            rng = galsim.BaseDeviate(random_seed + iobj + int(time.time()))

            # --- Update data dictionaries if scope=="object"
            if dist_scope == "object":
               galaxy_model_dico.update(self._read_galaxy_model_data(galaxy_model_name, 
                                                            job, worker, dist_scope=dist_scope))
               psf_model_dico.update(self._read_psf_model_data(psf_model_name, job, worker, 
                                                            dist_scope=dist_scope))
               shear_model_dico.update(self._read_shear_model_data(job, worker, 
                                                            dist_scope=dist_scope))
               if galaxy_noise_model_name != "no_noise":
                  galaxy_noise_model_dico.update(self._read_noise_model_data(
                           galaxy_noise_model_name, job, worker, dist_scope=dist_scope))

               # --- Updated Galaxy properties if scope=="*object
               gal_flux          = galaxy_model_dico["GAL_FLUX"]
               gal_e1_int        = galaxy_model_dico["GAL_E1"]
               gal_e2_int        = galaxy_model_dico["GAL_E2"]
               gal_initial_angle = galaxy_model_dico["GAL_INIT_ANGLE"]
               gal_cshift        = galaxy_model_dico["GAL_CSHIFT"]
          
               # --- Updated PSF properties if scope=="*object
               psf_flux          = psf_model_dico["PSF_FLUX"]
               psf_e1            = psf_model_dico["PSF_E1"]
               psf_e2            = psf_model_dico["PSF_E2"]
               psf_cshift        = psf_model_dico["PSF_CSHIFT"]
      
               # --- Updated Shear properties if scope=="*object
               gal_g1 = shear_model_dico["GAL_G1"]
               gal_g2 = shear_model_dico["GAL_G2"]



            # --------------------------------------- GALAXY --------------------------------------- 

            gal_cat_data_dico["gal_x"].append(gal_x)
            
            gal_cat_data_dico["gal_y"].append(gal_y)

            # --- Galaxy postage stamp
            gal_bounds = galsim.BoundsI(ix*gal_stamp_size+1 , (ix+1)*gal_stamp_size, 
                                        iy*gal_stamp_size+1 , (iy+1)*gal_stamp_size)
            sub_gal_image = gal_image[gal_bounds]

            # --- Compute initial disk position angle
            if irot % nb_gal_angles == 0:

               # Set index of realization
               ireal += 1

               # Set initial position angle in the ring
               gal_angle = gal_initial_angle
               gal_cat_data_dico["gal_pa"].append(gal_angle)
            
            else:

               # Increment position angle
               gal_angle += 180.0 / nb_gal_angles
               gal_cat_data_dico["gal_pa"].append(gal_angle)

            gal_cat_data_dico["iobj"].append(iobj)
            gal_cat_data_dico["irot"].append(irot)
            gal_cat_data_dico["ireal"].append(ireal)

            irot += 1

            # --- Compute whole galaxy centroid shift with respect to the geometric center
            shift_r = gal_cshift * pixel_scale
            gsud = galsim.UniformDeviate(rng)
            theta = gsud() * 2.0 * math.pi
            gdx = shift_r * math.cos(theta) 
            gdy = shift_r * math.sin(theta)
             
            # -- Create galaxy object based on model
            if galaxy_model_name == "gscsersic":   
               
               # --- Create single-component Sersic, GalSim implementation
               
               gal_disk_index = galaxy_model_dico["DISK_INDEX"]
               gal_disk_HLR   = galaxy_model_dico["DISK_HLR"]
               gal = galsim.Sersic(n=gal_disk_index, half_light_radius=gal_disk_HLR, flux=gal_flux,
                                   gsparams=gsparams)
#               print "***Flux:", gal.getFlux()

               # --- Apply reduced shear (1-r) / (1 + r) with additional rotation angle
               disk_intr_angle = 0.5 * math.atan2(gal_e2_int, gal_e1_int)
               gal_angle_rad = disk_intr_angle + gal_angle * math.pi / 180.0
               gal_disk_e_int  = math.hypot(gal_e1_int, gal_e2_int)
               gal_disk_e1_int = gal_disk_e_int * math.cos(2.0 * gal_angle_rad)
               gal_disk_e2_int = gal_disk_e_int * math.sin(2.0 * gal_angle_rad)
               gal = gal.shear(g1=gal_disk_e1_int, g2=gal_disk_e2_int)

#                print("angle:", gal_angle_rad, "angle deg:",  gal_angle_rad * 180.0 / math.pi, 
#                      "e:", gal_disk_e, "(e1, e2):", (gal_disk_e1, gal_disk_e2), 
#                      "intr angle:", disk_intr_angle, "intr angle:", disk_intr_angle * 180.0 / math.pi)
   
               gal_cat_data_dico["gal_flux"].append(gal_flux)
               gal_cat_data_dico["disk_flux"].append(gal_flux)
               gal_cat_data_dico["disk_n"].append(gal_disk_index)
               gal_cat_data_dico["disk_re"].append(gal_disk_HLR / pixel_scale)  # in pixels
               gal_cat_data_dico["gal_e_int"].append(gal_disk_e_int)
               gal_cat_data_dico["gal_e1_int"].append(gal_disk_e1_int)
               gal_cat_data_dico["gal_e2_int"].append(gal_disk_e2_int)

            if galaxy_model_name == "gbdsersic_r":   

               # --- Create bulge + disk Sersic, GalSim implementation

               # --- Disk
               gal_disk_index = galaxy_model_dico["DISK_INDEX"]
               gal_disk_HLR   = galaxy_model_dico["DISK_HLR"]
               gal_disk_flux_frac = galaxy_model_dico["DISK_FLUX_FRAC"]
               disk_flux = gal_disk_flux_frac * gal_flux # disk flux fraction
               disk  = galsim.Sersic(n=gal_disk_index, half_light_radius=gal_disk_HLR, 
                                     flux=disk_flux, gsparams=gsparams)

               gal_cat_data_dico["gal_flux"].append(gal_flux)
               gal_cat_data_dico["disk_flux"].append(disk_flux)
               gal_cat_data_dico["disk_n"].append(gal_disk_index)
               gal_cat_data_dico["disk_re"].append(gal_disk_HLR / pixel_scale) # in pixels

               # --- Bulge
               if gal_disk_flux_frac < 1.0:
                  bulge_flux = gal_flux - disk_flux
                  gal_bulge_index = galaxy_model_dico["BULGE_INDEX"]
                  gal_bulge_HLR   = galaxy_model_dico["BULGE_HLR"]
                  bulge = galsim.Sersic(n=gal_bulge_index, half_light_radius=gal_bulge_HLR, 
                                        flux=bulge_flux, gsparams=gsparams)

                  gal_cat_data_dico["bulge_flux"].append(bulge_flux)
                  gal_cat_data_dico["bulge_n"].append(gal_bulge_index)
                  gal_cat_data_dico["bulge_re"].append(gal_bulge_HLR / pixel_scale)   # in pixels

                  bulge_shift = galaxy_model_dico["BULGE_CSHIFT"] # relative shift for the bulge
                  if bulge_shift > 0:
                     bulge_shift_r = pixel_scale * bulge_shift
                     gbsud = galsim.UniformDeviate(rng)
                     btheta = gbsud() * 2. * math.pi
                     bdx = bulge_shift_r * math.cos(btheta)
                     bdy = bulge_shift_r * math.sin(btheta)
                     bulge = bulge.shift(bdx, bdy) 

                  else:
                     bdx = bdy = 0

                  gal_cat_data_dico["bulge_cshift_x"].append(bdx)
                  gal_cat_data_dico["bulge_cshift_y"].append(bdy)

                  gal_cat_data_dico["bulge_xc"].append(gal_xc + bdx)
                  gal_cat_data_dico["bulge_yc"].append(gal_yc + bdy)

                  # --- Bulge + Disk
                  gal_intr_angle = 0.5 * math.atan2(gal_e2_int, gal_e1_int)
                  gal_angle_rad = gal_intr_angle + gal_angle * math.pi / 180.0
                  gal_e_int = math.hypot(gal_e1_int, gal_e2_int)
                  gal_bd_e1_int = gal_e_int * math.cos(2.0 * gal_angle_rad)
                  gal_bd_e2_int = gal_e_int * math.sin(2.0 * gal_angle_rad)

                  gal = galsim.Add([disk, bulge])
                  gal = gal.withFlux(gal_flux)
                  gal = gal.shear(g1=gal_e1_int, g2=gal_e2_int)

                  gal_cat_data_dico["gal_e_int"].append(gal_e_int)
                  gal_cat_data_dico["gal_e1_int"].append(gal_bd_e1_int)
                  gal_cat_data_dico["gal_e2_int"].append(gal_bd_e2_int)

               else:
                  self.helper.print_error("The galaxy has no bulge, aborting")
                  success = False
                  return success     

            if galaxy_model_name == "gaussian":

               # --- Create single-component Gaussian
               gal_disk_sigma = galaxy_model_dico["DISK_SIGMA"]
               gal = galsim.Gaussian(sigma=gal_disk_sigma, flux=gal_flux, gsparams=gsparams)

               gal_disk_HLR   = gal.getHalfLightRadius() # arcsecs

               # --- Apply reduced shear (1-r) / (1 + r) with additional rotation angle
               disk_intr_angle = 0.5 * math.atan2(gal_e2_int, gal_e1_int)
               gal_angle_rad = disk_intr_angle + gal_angle * math.pi / 180.0
               gal_disk_e_int = math.hypot(gal_e1_int, gal_e2_int)
               gal_disk_e1_int = gal_disk_e_int * math.cos(2.0 * gal_angle_rad)
               gal_disk_e2_int = gal_disk_e_int * math.sin(2.0 * gal_angle_rad)
               gal = gal.shear(g1=gal_disk_e1_int, g2=gal_disk_e2_int)

               #print("angle:", gal_angle_rad, "angle deg:",  gal_angle_rad * 180.0 / math.pi, 
               #      "e:", gal_disk_e, "(e1, e2):", (gal_disk_e1, gal_disk_e2), 
               #      "intr angle:", disk_intr_angle, "intr angle:", disk_intr_angle * 180.0 / math.pi)
   

               gal_cat_data_dico["gal_flux"].append(gal_flux)
               gal_cat_data_dico["disk_flux"].append(gal_flux)
               gal_cat_data_dico["disk_n"].append(0.5)  # Gaussian equivalent to Sersic with n=0.5
               gal_cat_data_dico["disk_sigma"].append(gal_disk_sigma / pixel_scale) # in pixels 
               gal_cat_data_dico["disk_re"].append(gal_disk_HLR / pixel_scale)      # in pixels
               gal_cat_data_dico["gal_e_int"].append(gal_disk_e_int)
               gal_cat_data_dico["gal_e1_int"].append(gal_disk_e1_int)
               gal_cat_data_dico["gal_e2_int"].append(gal_disk_e2_int)


            # --- End if galaxy_model_name

            # --- Apply centroid shift (in arcsecs)
            gal = gal.shift(gdx, gdy)  

            
            # Actual galaxy centroid
            gal_cat_data_dico["gal_cshift_x"].append(gdx / pixel_scale)       # in pixels
            gal_cat_data_dico["gal_cshift_y"].append(gdy / pixel_scale)       # in pixels
            gal_cat_data_dico["gal_xc"].append(gal_xc + gdx / pixel_scale)    # in pixels
            gal_cat_data_dico["gal_yc"].append(gal_yc + gdy / pixel_scale)    # in pixels

            # --- Apply shear
            gal = gal.shear(g1=gal_g1, g2=gal_g2)
            ###print "shear:", (g1, g2)
            gal_cat_data_dico["gal_g1"].append(gal_g1)
            gal_cat_data_dico["gal_g2"].append(gal_g2)

            # --- Sheared (observed) ellipticities
            gal_e1 = gal_cat_data_dico["gal_e1_int"][-1] + gal_g1
            gal_e2 = gal_cat_data_dico["gal_e2_int"][-1] + gal_g2
            gal_cat_data_dico["gal_e1"].append(gal_e1)
            gal_cat_data_dico["gal_e2"].append(gal_e2)
            gal_cat_data_dico["gal_e"].append(math.hypot(gal_e1, gal_e2))

            # ---------------------------------------- PSF -----------------------------------------   

            psf_cat_data_dico["iobj"].append(iobj)
            psf_cat_data_dico["irot"].append(irot)
            psf_cat_data_dico["ireal"].append(ireal)

            # --- Create PSF object
            psf_model = worker.config.get_as_string("MODEL_NAME", "PSF.MODEL")
            if psf_model == "gmoffat":

               # --- Moffat PSF
               psf_FWHM = psf_model_dico["PSF_FWHM"]
               psf_beta = psf_model_dico["PSF_BETA"]
               psf = galsim.Moffat(beta=psf_beta, fwhm=psf_FWHM, flux=psf_flux, 
                                   trunc=0,  gsparams=gsparams)
               psf = psf.shear(e1=psf_e1, e2=psf_e2)

               psf_cat_data_dico["psf_beta"].append(psf_beta)
      
            elif psf_model == "gaussian":

               # --- Gaussian PSF
               psf_FWHM = psf_model_dico["PSF_FWHM"]
               psf = galsim.Gaussian(fwhm=psf_FWHM, flux=psf_flux, gsparams=gsparams)
               psf = psf.shear(e1=psf_e1, e2=psf_e2)
            
            psf_cat_data_dico["psf_flux"].append(psf_flux)
            psf_cat_data_dico["psf_FWHM"].append(psf_FWHM)
            psf_cat_data_dico["psf_e1"].append(psf_e1)
            psf_cat_data_dico["psf_e2"].append(psf_e2)
            psf_cat_data_dico["psf_e"].append(math.hypot(psf_e1, psf_e2))      

            # --- Compute and apply PSF centroid shift
            shift_r = pixel_scale * psf_cshift
            sud = galsim.UniformDeviate(rng)
            theta = sud() * 2.0 * math.pi
            pdx = shift_r * math.cos(theta)
            pdy = shift_r * math.sin(theta)            
            psf = psf.shift(pdx, pdy)

            psf_cat_data_dico["psf_x"].append(psf_x)
            psf_cat_data_dico["psf_y"].append(psf_y)
            
            psf_cat_data_dico["psf_cshift_x"].append(pdx / pixel_scale)
            psf_cat_data_dico["psf_cshift_y"].append(pdy / pixel_scale)
            psf_cat_data_dico["psf_xc"].append(psf_xc + pdx / pixel_scale)
            psf_cat_data_dico["psf_yc"].append(psf_yc + pdy / pixel_scale)

            # --- Create PSF postage stamp in PSF mosaic
            psf_bounds = galsim.BoundsI(ix*psf_stamp_size+1 , (ix+1)*psf_stamp_size, 
                                        iy*psf_stamp_size+1 , (iy+1)*psf_stamp_size)
            sub_psf_image = psf_image[psf_bounds]

            # --- Draw the PSF image:
            psf.drawImage(sub_psf_image, method="real_space")

#             #### TEST ####
#             gal.drawImage(sub_gal_image, method='auto') 
#             cut_gal_array = self.cut_stamp_around_centroid(sub_gal_image.array, 64, (32.0, 32.0))
#             print "truncated gal flux:", numpy.sum(cut_gal_array), "peak:", numpy.max(cut_gal_array)    
#   
#             cut_gal_image = galsim.ImageF(cut_gal_array)
#             cut_gal = galsim.InterpolatedImage(cut_gal_image, 
#                                                       x_interpolant = "quintic",
#                                                       scale=pixel_scale)            
#             cut_psf_array = self.cut_stamp_around_centroid(sub_psf_image.array, 64, (32.0, 32.0))
#             print "truncated psf flux:", numpy.sum(cut_psf_array), "peak:", numpy.max(cut_psf_array)   
#             cut_psf_image = galsim.ImageF(cut_psf_array)
#             cut_psf = galsim.InterpolatedImage(cut_psf_image, 
#                                                       x_interpolant = "quintic",
#                                                       scale=pixel_scale)            
#              
#             cut_conv = galsim.Convolve([cut_psf, cut_gal])
#             cut_conv.drawImage(sub_gal_image, method='auto')
# #             cut_conv.drawImage(sub_psf_image, method="real_space")         
#             print "=> flux convolved:", numpy.sum(sub_gal_image.array), "peak:", numpy.max(sub_gal_image.array)   
#             
#             conv_array = scipy.signal.fftconvolve(cut_psf_array, cut_gal_array, mode='same')            
#             print "=> flux convolved scipy:", numpy.sum(conv_array), "peak:", numpy.max(conv_array)     
#             
#             #### TEST ####


            # ------------------------------------- CONVOLUTION  -----------------------------------   

            # --- Convolve PSF and galaxy
            if apply_PSF:
               if convol_method == "galsim_fft":
                  conv_gal = galsim.Convolve([psf, gal])

               elif convol_method == "scipy_fft":
                  # --- Convolve with scipy FFT convolve() and interpolate to get a galsim image   
                  conv_gal_array = scipy.signal.fftconvolve(sub_psf_image.array.astype(numpy.float),
                                                            sub_gal_image.array.astype(numpy.float), 
                                                            mode='same')
                  conv_gal_image = galsim.ImageF(conv_gal_array)
                  conv_gal = galsim.InterpolatedImage(conv_gal_image, 
                                                      x_interpolant = "quintic",
                                                      scale=pixel_scale)
            else:
               conv_gal = gal

            # --- Draw the convolved galaxy image:
            conv_gal.drawImage(sub_gal_image, method='auto')

            # ---------------------------------------- NOISE  --------------------------------------
            #no_noise_image = sub_gal_image.copy()
            
            if galaxy_noise_model_name != "no_noise":
               sub_gal_image, sigma_noise = self._apply_noise_model(sub_gal_image, 
                                                    galaxy_noise_model_name, galaxy_noise_model_dico, 
                                                    pixel_scale, rng, job, worker)
               gal_cat_data_dico["gal_wanted_SNR"].append(galaxy_noise_model_dico["WANTED_SNR"])
               gal_cat_data_dico["gal_noise_sigma"].append(sigma_noise)
            else:   
               gal_cat_data_dico["gal_wanted_SNR"].append(-1)
               gal_cat_data_dico["gal_noise_sigma"].append(0)

#             rms = numpy.sqrt(numpy.mean((no_noise_image.array - sub_gal_image.array)**2))
#             print "rms:", rms, numpy.max(sub_gal_image.array) / rms, numpy.sum(sub_gal_image.array) / rms   
#             rms = numpy.sqrt(numpy.mean((sub_gal_image.array - numpy.mean(sub_gal_image.array))**2))
#             print rms, numpy.std(sub_gal_image.array)

            #gal_flux_convolved = numpy.sum(sub_gal_image.array)
            #print gal_angle, "estimated flux:",  gal_flux_estimated, "flux convolved:", gal_flux_convolved, "peak estimated: %.9f" %(peak_estimated)

            gal_flux_estimated = numpy.sum(sub_gal_image.array)
            gal_cat_data_dico["gal_flux_estim"].append(gal_flux_estimated)
            gal_peak_estimated = numpy.max(sub_gal_image.array)
            gal_cat_data_dico["gal_peak"].append(gal_peak_estimated)
            gal_cat_data_dico["gal_mean"].append(numpy.mean(sub_gal_image.array))
            gal_cat_data_dico["gal_std"].append(numpy.std(sub_gal_image.array))
            
            psf_flux_estimated = numpy.sum(sub_psf_image.array)
            psf_cat_data_dico["psf_flux_estim"].append(psf_flux_estimated)
            psf_peak_estimated = numpy.max(sub_psf_image.array)
            psf_cat_data_dico["psf_peak"].append(psf_peak_estimated)
            psf_cat_data_dico["psf_mean"].append(numpy.mean(sub_psf_image.array))
            psf_cat_data_dico["psf_std"].append(numpy.std(sub_psf_image.array))
            
            #print "Mode:", numpy.max(sub_gal_image.array)

            #### TEST ####
            #cut_gal = self.cut_stamp_around_centroid(sub_gal_image.array, 22, (32.0, 32.0))
            #print "truncated gal flux:", numpy.sum(cut_gal)   
            #cut_psf = self.cut_stamp_around_centroid(sub_psf_image.array, 22, (32.0, 32.0))
            #print "truncated psf flux:", numpy.sum(cut_psf)   
            #### TEST ####

            gal_cat_data_dico["gal_stamp_size"].append(gal_stamp_size)
            psf_cat_data_dico["psf_stamp_size"].append(psf_stamp_size)

            if nb_objects >= 10 and iobj % 10 == 0:
               if worker.logging_enabled():
                  worker.logger.log_info_p("{0} - {1:03d}-{2:1d} - {3} objects created".format(
                     worker.name, job.img_no, job.epoch, iobj))               

            
            # --- Next object
            iobj += 1

            if iobj >= nb_objects:
               break 
         # --- end for          

         if iobj >= nb_objects:
            break 

      # ---end for

      if worker.logging_enabled():
         worker.logger.log_info_p("{0} - {1:03d}-{2:1d} - {3} objects created".format(
            worker.name, job.img_no, job.epoch, iobj))               

      if worker.logging_enabled():
         worker.logger.log_info_p("Creating galaxy and PSF mosaics...")

      # --- Save galaxy and PSF mosaics to disk
      gal_image.write(gal_output_image_filepath, clobber=True)
      psf_image.write(psf_output_image_filepath, clobber=True)

      # --- Create galaxy and PSF catalogs
      self._create_galaxy_catalog(galaxy_model_name, 
                                  gal_cat_data_dico, gal_output_catalog_filepath,
                                  job, worker)
      self._create_psf_catalog(psf_model_name, 
                                  psf_cat_data_dico, psf_output_catalog_filepath,
                                  job, worker)
      
      
   # -----------------------------------------------------------------------------------------------
   def cut_stamp_around_centroid(self, stamp, stamp_size, (row, col)):
      """!  Cut a postage stamp of a given size around a given centroid relative to the image """
      
      half_stamp_size = stamp_size / 2.0
      if stamp_size % 2 == 0:     
         # --- Even stamp size
         ####print "EVEN stamp size:", stamp_size, "centroid:", (row, col), "half:", half_stamp_size, row-half_stamp_size
         return stamp[(int)(row-half_stamp_size):int(row+half_stamp_size), 
                      int(col-half_stamp_size):int(col+half_stamp_size)]
      else:
         # --- Odd stamp size
         #print "ODD stamp size:", stamp_size, "centroid:", (row, col), half_stamp_size, row-half_stamp_size+0.5
         return stamp[int(row-half_stamp_size+0.5):int(row+half_stamp_size+0.5), 
                      int(col-half_stamp_size+0.5):int(col+half_stamp_size+0.5)]


   # -----------------------------------------------------------------------------------------------
   def _apply_noise_model(self, image, noise_model_name, noise_model_dico, pixel_scale, 
                                rng, job, worker):
      """!
         Apply a noise model to an object
         @param image target GalSim image 
         @param noise_model_name noise model to apply (see configuration)
         @param noise_model_dico dictionary containing the noise model parameters
         @return object with noise added 
      """
      
      noise_sigma = 0.0
      
      if noise_model_name == "no_noise":
         return image, noise_sigma
      
      nud = galsim.UniformDeviate(rng)
         
      # --- Simple Gaussian noise
      if noise_model_name == "simple_gaussian":

         wanted_SNR = noise_model_dico["WANTED_SNR"]            
         additional_sky_pixels  = noise_model_dico["ADDITIONAL_SKY"] * pixel_scale**2  

         sigma_noise = noise_model_dico["SIGMA_NOISE"]

         #print "noise_model_dico:", noise_model_dico
         #print "### additional_sky_pixels:", additional_sky_pixels
         #print "sigma_noise:", sigma_noise


         image += additional_sky_pixels
         noise = galsim.GaussianNoise(nud, sigma=sigma_noise)
         if wanted_SNR != -1:               
            image.addNoiseSNR(noise, snr=wanted_SNR)
         else:   
            image.addNoise(noise)

         noise_sigma = noise.getSigma()

      # --- Simple Poisson noise
      elif noise_model_name == "simple_poisson":

         wanted_SNR = noise_model_dico["WANTED_SNR"]
         substracted_sky_pixels = noise_model_dico["SUBSTRACTED_SKY"] * pixel_scale**2           
         additional_sky_pixels  = noise_model_dico["ADDITIONAL_SKY"] * pixel_scale**2  

         image += additional_sky_pixels
         noise = galsim.PoissonNoise(nud, sky_level=substracted_sky_pixels)               
         if wanted_SNR != -1:               
            image.addNoiseSNR(noise, snr=wanted_SNR)
         else:   
            image.addNoise(noise)

         noise_sigma = math.sqrt(noise.getVariance())

      # --- More realistic CCD noise (Gaussian Read noise + Poisson noise)
      elif noise_model_name == "CCD_noise":

         gain = noise_model_dico["GAIN"]
         read_noise = noise_model_dico["READ_NOISE"]
         wanted_SNR = noise_model_dico["WANTED_SNR"]
         substracted_sky_pixels = noise_model_dico["SUBSTRACTED_SKY"] * pixel_scale**2           
         additional_sky_pixels  = noise_model_dico["ADDITIONAL_SKY"]  * pixel_scale**2  

         image += additional_sky_pixels
         noise = galsim.CCDNoise(nud, gain=gain, 
                                     read_noise=read_noise, sky_level=substracted_sky_pixels)               
         if wanted_SNR != -1:               
            image.addNoiseSNR(noise, snr=wanted_SNR)
         else:   
            image.addNoise(noise)

         noise_sigma = math.sqrt(noise.getVariance())
            
      else:      
         if worker.logging_enabled():
            worker.logger.log_error_p("{0} - {1:03d}-{2:1d} - "\
                  "Unrecognized noise model: {3}".format(
                     worker.name, job.img_no, job.epoch, noise_model_name))                  

      return image, noise_sigma
   

   # -----------------------------------------------------------------------------------------------
   def _read_sim_model_data(self, job, worker, dist_scope=None):
      """! 
         Return a dictionary with the galaxy model data instantiated in the case of a "field" scope
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return dictionary with galaxy model data
      """
      data_dico = {}
      section_name = "SIMULATION"

      data_dico["NB_FIELDS"]          = worker.config.get_as_int("NB_FIELDS", section_name)
      data_dico["NB_REALIZATIONS"]    = worker.config.get_as_int("NB_REALIZATIONS", section_name)
      data_dico["NB_ANGLES_RING"]     = worker.config.get_as_int("NB_ANGLES_RING", section_name)
      data_dico["APPLY_PSF"] = worker.config.get_as_boolean("APPLY_PSF", section_name)
      data_dico["CONVOLUTION_METHOD"] = worker.config.get_as_string("CONVOLUTION_METHOD", 
                                                                  section_name)

      return data_dico        

   # -----------------------------------------------------------------------------------------------
   def _read_galaxy_model_data(self, model_name, job, worker, dist_scope=None):
      """! 
         Return a dictionary with the galaxy model data instantiated in the case of a "field" scope
         @param model_name name of galaxy model
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return dictionary with galaxy model data
      """
      model_dico = {}
      section_name = "GALAXY.MODEL.{0}".format(model_name.upper())
      section_dico = worker.config.get_section_data(section_name)

      for param_name in section_dico.keys():
         
         dist_info_dico = eval(section_dico[param_name])
         self._read_dist_data(dist_info_dico, section_name, param_name, 
                              dist_scope, model_dico, job, worker)     
            
      return model_dico        
      
   # -----------------------------------------------------------------------------------------------
   def _read_psf_model_data(self, model_name, job, worker, dist_scope=None):
      """! 
         Return a dictionary with the galaxy model data instantiated in the case of a "field" scope
         @param model_name name of PSF model
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return dictionary with galaxy model data
      """
      model_dico = {}
      section_name = "PSF.MODEL.{0}".format(model_name.upper())
      section_dico = worker.config.get_section_data(section_name)

      for param_name in section_dico.keys():
         
         dist_info_dico = eval(section_dico[param_name])
         self._read_dist_data(dist_info_dico, section_name, param_name, 
                              dist_scope, model_dico, job, worker)            
      
      return model_dico        
      
   # -----------------------------------------------------------------------------------------------
   def _read_shear_model_data(self, job, worker, dist_scope=None):
      """! 
         Return a dictionary with the shear model data instantiated in the case of a "field" scope
         @param model_name name of shear model
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return dictionary with galaxy model data
      """
      model_dico = {}
      section_name = "SHEAR"
      section_dico = worker.config.get_section_data(section_name)

      for param_name in section_dico.keys():
         
         dist_info_dico = eval(section_dico[param_name])
         self._read_dist_data(dist_info_dico, section_name, param_name, 
                              dist_scope, model_dico, job, worker)     
            
      return model_dico          

   # -----------------------------------------------------------------------------------------------
   def _read_noise_model_data(self, model_name, job, worker, dist_scope=None):
      """! 
         Return a dictionary with the galaxy noise model data instantiated in the case of
         a "field" scope
         @param model_name name of galaxy model
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return dictionary with galaxy model data
      """
      model_dico = {}
      if model_name != "no_noise":
         section_name = "NOISE.MODEL.{0}".format(model_name.upper())
         section_dico = worker.config.get_section_data(section_name)
   
         for param_name in section_dico.keys():
            
            dist_info_dico = eval(section_dico[param_name])
            self._read_dist_data(dist_info_dico, section_name, param_name, 
                                 dist_scope, model_dico, job, worker)     
               
      return model_dico    

   # -----------------------------------------------------------------------------------------------
   def _read_dist_data(self, dist_info_dico, section_name, param_name,  
                             dist_scope, model_dico, job, worker):   

         try:
            
            dist_name = dist_info_dico.get("dist", None)
            if dist_name is None:
               dist_name = "constant"    # if unspecified assume constant
               if worker.logging_enabled():
                  worker.logger.log_warning_p("{0} - Distribution not specified for [{1}].{2}, "\
                                            "assume 'constant' distribution...".format(
                                                         worker.name, section_name, param_name))
                  
            input_dist_scope = dist_info_dico.get("scope", None)

            # --- Compute random value for the given input scope
            rand_value = None
            
            if dist_scope is None:
               rand_value = self._get_random_value(dist_info_dico, section_name, param_name,
                                                   job, worker)
               if not rand_value is None:
                  model_dico[param_name] = rand_value
               else:
                  if worker.logging_enabled():
                     worker.logger.log_error_p("{0} - The random value for [{1}].{2} "\
                                               "could not be obtained".format(
                                                            worker.name, section_name, param_name))
                    
            elif input_dist_scope == dist_scope:

               # --- Update random value using specified scope
               dist = dist_info_dico.get("dist", None)
               value_range  = dist_info_dico.get("range", None)
               
               max_tries = 10
               for itry in range(0, max_tries):
                  rand_value = self._get_random_value(dist_info_dico, section_name, param_name,
                                                      job, worker)
                  # --- Check that computed random value lies within the specified range 
                  if not value_range is None:
                     if rand_value < value_range[0] or rand_value > value_range[1]:
                        if worker.logging_enabled():
                           worker.logger.log_warning_p(
                                       "{0} - Computed random value {1} for distribution' {2}'"\
                                       " must lie in range {3} - Retrying ({4}/{5})...".format(
                                       worker.name, rand_value, dist.lower(), value_range, 
                                       itry+1, max_tries))
                           rand_value = None
                     else:   
                        break
               
               if not rand_value is None:
                  model_dico[param_name] = rand_value
               else:   
                  if worker.logging_enabled():
                     worker.logger.log_error_p("{0} - Invalid or incomplete data in [{1}].{2}"\
                              "...".format(worker.name, section_name, param_name))
            
         except Exception as detail:
            if worker.logging_enabled():
               worker.logger.log_error_p("{0} - Invalid or incomplete data in [{1}].{2} - {3}"\
                        "...".format(worker.name, section_name, param_name, detail))
      
      
   # -----------------------------------------------------------------------------------------------
   def _get_random_value(self, dist_info_dico, section_name, param_name, job, worker):
      """! 
         Return a random variable value drawn from a given probability distribution
         @param dist_info_dico dictionary containing the definition of the random variable
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return computed random value based on the input distribution
      """
      rand_value = None
      
      #seed = int(job.img_no * math.ceil(time.time()))
      numpy.random.seed()
      
      dist = dist_info_dico.get("dist", None)
      params = dist_info_dico.get("params", None)
      if dist is None:
         dist = "constant"
         
      if params is None or not isinstance(params, list) or len(params) == 0:
         if worker.logging_enabled():
            worker.logger.log_warning_p("{0} - No parameter given for distribution '{1}'".format(
                                                   worker.name, dist.lower()))
         return None
      
      if dist.lower() == "constant":           # constant distribution 
         # --- Constant distribution: value remains fixed
         rand_value = float(params[0])
      
      elif dist.lower() == "normal":           # Normal distribution
         mean   = float(params[0])
         stddev = float(params[1])    
         rand_value = numpy.random.randn() * stddev + mean           
      
      elif dist.lower() == "uniform":          # uniform distribution
         dist_min = float(params[0])
         dist_max = float(params[1])    
         rand_value = numpy.random.random() * (dist_max - dist_min) + dist_min

      elif(dist.lower() == 'rayleigh'):
         sigma = float(params[0])
         rand_value =  numpy.random.rayleigh(sigma)
           
      elif dist.lower() == "lognormal":        # lognormal distribution
         mean   = float(params[0])
         stddev = float(params[1])
         rand_value = 10 ** (numpy.random.randn() * stddev + mean)      

#       elif dist.lower() == "lognormal-mean":     # adjusted lognormal distribution
#          # --- Apply correction factor to give the correct mean
#          mean   = float(params[0])
#          stddev = float(params[1])
#          rand_value = 10 ** (numpy.random.randn() * stddev + mean) * \
#                              numpy.exp(-((stddev * numpy.log(10)) ** 2) / 2)
                             
      elif dist.lower() == "sequence":     # adjusted lognormal distribution
         rand_value =  job._sequence_value_dico.get((section_name, param_name), None)
#          print "*** rand_value", job.name, section_name, param_name, rand_value       
                            
#       elif dist.lowwer() == "ps":               # a power spectrum of the form k**n 
                             
                             
      else:
         # --- Distribution not understood: assume constant
         rand_value = float(params[0])
         if worker.logging_enabled():
            worker.logger.log_warning_p(
                                 "{0} - Unknown distribution '{1}', assume 'constant' ".format(
                                                   worker.name, dist.lower()))
 
#       # --- Check that computed random value lies within the specified range 
#       value_range  = dist_info_dico.get("range", None)
#       if not value_range is None:
#          if rand_value < value_range[0] or rand_value > value_range[1]:
#             if worker.logging_enabled():
#                worker.logger.log_warning_p("{0} - Computed random value {1} for distribution' {2}'"\
#                                            " must lie in range {3}".format(
#                                              worker.name, rand_value, dist.lower(), value_range))
#             rand_value = None

      return rand_value   


   # -----------------------------------------------------------------------------------------------
   def _create_gal_cat_data_dico(self, model_name, job, worker):
      """! 
         Create the dictionary that contains the simulated data used for the galaxy simulations; 
         it will be used to generate the simulation catalogs
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return empty galaxy data dictionary
      """
      # --- Possible models: -----------------------------------------------------------------------
      # scsersic        : basic single component Sersic profile implementation in Python
      # gscsersic       : single component Sersic profile implemented using GalSim exponential disk
      # gbdsersic_r     : bulge+disk Sersic profile implemented using GalSim (reduced shear)
      # gbdsersic_d     : bulge+disk Sersic profile implemented using GalSim (distortion)
      # gaussian        : Gaussian profile implemented using GalSim (Sersic with index=0.5)
      # --------------------------------------------------------------------------------------------

      data_dico = OrderedDict()

      var_list = []

      # --- Column names depending on model
      if model_name in ["scsersic", "gscsersic"]:
         var_list = ["iobj", "ireal", "irot", "gal_x", "gal_y", "gal_xc", "gal_yc",\
                      "gal_pa", "gal_flux", "gal_flux_estim", "disk_flux",\
                     "disk_n", "disk_re", "gal_e_int", "gal_e1_int", "gal_e2_int",\
                      "gal_g1", "gal_g2", "gal_e", "gal_e1", "gal_e2",\
                     "gal_cshift_x", "gal_cshift_y",\
                     "gal_wanted_SNR", "gal_noise_sigma", "gal_stamp_size",\
                     "gal_peak", "gal_mean", "gal_std"]

      elif model_name in ["gbdsersic_r", "gbdsersic_d"]:
         var_list = ["iobj", "ireal", "irot", "gal_x", "gal_y",\
                     "gal_xc", "gal_yc", "bulge_xc", "bulge_yc",\
                     "gal_pa", "gal_flux", "gal_flux_estim", "disk_flux", "bulge_flux",\
                     "disk_n", "bulge_n", "disk_re", "bulge_re",\
                     "gal_e_int", "gal_e1_int", "gal_e2_int",\
                     "gal_g1", "gal_g2", "gal_e", "gal_e1", "gal_e2",\
                     "gal_cshift_x", "gal_cshift_y", "bulge_cshift_x", "bulge_cshift_y",\
                     "gal_wanted_SNR", "gal_noise_sigma", "gal_stamp_size",\
                     "gal_peak", "gal_mean", "gal_std"]

      elif model_name == "gaussian":
         var_list = ["iobj", "ireal", "irot", "gal_x", "gal_y", "gal_xc", "gal_yc",\
                     "gal_pa", "gal_flux", "gal_flux_estim", "disk_flux",\
                     "disk_n", "disk_re", "disk_sigma", "gal_e_int", "gal_e1_int", "gal_e2_int",\
                     "gal_g1", "gal_g2", "gal_e", "gal_e1", "gal_e2",\
                     "gal_cshift_x", "gal_cshift_y",\
                     "gal_wanted_SNR", "gal_noise_sigma", "gal_stamp_size",\
                     "gal_peak", "gal_mean", "gal_std"]            

      for var in var_list:
         data_dico[var] = []  

      return data_dico

   # -----------------------------------------------------------------------------------------------
   def _create_psf_cat_data_dico(self, model_name, job, worker):
      """!
         Create the dictionary that contains the PSF simulated data used for the simulations; 
         it will be used to generate the simulation catalogs
         @param job an object of class MksJob to process
         @param worker instance of the worker process         
         @return empty galaxy data dictionary
      """
      
      # --- Profile Parameter Definitions ----------------------------------------------------------
      # PSF_FWHM       : PSF FWHM in arc seconds
      # PSF_BETA       : PSF Moffat beta
      # PSF_E1, PSF_E2 : PSF ellipticity components
      # PSF_INIT_ANGLE : PSF initial position angles in a ring  
      # PSF FLUX       : PSF total flux in ADUs
      # --------------------------------------------------------------------------------------------
         
      data_dico = OrderedDict()
      var_list = []
      # --- Column names depending on model
      
      if model_name == "gmoffat":
         var_list = ["iobj", "ireal", "irot", "psf_x", "psf_y", "psf_xc", "psf_yc", \
                     "psf_flux", "psf_flux_estim", "psf_FWHM", "psf_beta", \
                     "psf_e", "psf_e1", "psf_e2", \
                     "psf_cshift_x", "psf_cshift_y", \
                     "psf_stamp_size", "psf_peak", "psf_mean", "psf_std"]
         
      elif model_name == "gaussian":
         var_list = ["iobj", "ireal", "irot", "psf_x", "psf_y", "psf_xc", "psf_yc", \
                     "psf_flux", "psf_flux_estim", "psf_FWHM", \
                     "psf_e", "psf_e1", "psf_e2", \
                     "psf_cshift_x", "psf_cshift_y", \
                     "psf_stamp_size", "psf_peak", "psf_mean", "psf_std"]

      for var in var_list:
         data_dico[var] = []   

      return data_dico


   # -----------------------------------------------------------------------------------------------
   def _create_galaxy_catalog(self, model_name, data_dico, output_filepath, job, worker):

      col_key_names   = data_dico.keys()   # ordered list of column names
      if "gal_trunc_size" in col_key_names:
         col_key_names.remove("gal_trunc_size")
      
      col_key_indexes = numpy.arange(0, len(col_key_names))
      col_key_map     = dict(zip(data_dico.keys(), col_key_indexes))

      col_fmt_map = {}
      col_integers      = ["ireal", "gal_x", "gal_y", "gal_stamp_size"]
      col_small_floats  = ["irot", "gal_SNR", "disk_n", "bulge_n"]
      col_ellipticities = ["gal_e1_int", "gal_e2_int", "gal_e"]
      col_fmt_map = col_fmt_map.fromkeys(col_integers, "%d")
      col_fmt_map.update(col_fmt_map.fromkeys(col_small_floats, "%.1f"))
      col_fmt_map.update(col_fmt_map.fromkeys(col_ellipticities, "%+.9e"))

      output_directory, output_filename = os.path.split(output_filepath)
      hdu_no   = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.CATALOG_PROPERTIES")
      is_ascii = worker.config.get_as_boolean("IS_ASCII_FORMAT",
                                              "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY")
      is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                   "PRIMARY_DATASET.CATALOG_PROPERTIES.GALAXY")
      
      self.helper.create_from_list_dico(data_dico, output_directory, output_filename,
                                        job, worker,
                                        col_list=col_key_names,   
                                        key_index_map=col_key_map, key_fmt_map=col_fmt_map, 
                                        default_fmt="%.6e", 
                                        is_sextractor=is_sextractor, is_ascii=is_ascii, 
                                        hdu_no=hdu_no)

   # -----------------------------------------------------------------------------------------------
   def _create_psf_catalog(self, model_name, data_dico, output_filepath, job, worker):

      col_key_names   = data_dico.keys()   # ordered list of column names
      if "psf_trunc_size" in col_key_names:
         col_key_names.remove("psf_trunc_size")
      col_key_indexes = numpy.arange(0, len(col_key_names))
      col_key_map     = dict(zip(data_dico.keys(), col_key_indexes))

      col_fmt_map = {}
      col_integers      = ["ireal", "psf_x", "psf_y", "psf_stamp_size"]
      col_small_floats  = ["irot", "psf_SNR", "psf_FWHM", "psf_beta"]
      col_ellipticities = ["psf_e1", "psf_e2", "psf_e"]
      col_fmt_map = col_fmt_map.fromkeys(col_integers, "%d")
      col_fmt_map.update(col_fmt_map.fromkeys(col_small_floats, "%.1f"))
      col_fmt_map.update(col_fmt_map.fromkeys(col_ellipticities, "%+.9e"))

      output_directory, output_filename = os.path.split(output_filepath)
      hdu_no   = worker.config.get_as_int("HDU_NO", "PRIMARY_DATASET.CATALOG_PROPERTIES")
      is_ascii = worker.config.get_as_boolean("IS_ASCII_FORMAT",
                                              "PRIMARY_DATASET.CATALOG_PROPERTIES.STAR")
      is_sextractor = worker.config.get_as_boolean("IS_SEXTRACTOR_FORMAT", 
                                                   "PRIMARY_DATASET.CATALOG_PROPERTIES.STAR")

      self.helper.create_from_list_dico(data_dico, output_directory, output_filename,
                                        job, worker,
                                        col_list=col_key_names,   
                                        key_index_map=col_key_map, key_fmt_map=col_fmt_map, 
                                        default_fmt="%.6e", 
                                        is_sextractor=is_sextractor, is_ascii=is_ascii, 
                                        hdu_no=hdu_no)

   # -----------------------------------------------------------------------------------------------
   def _adjust_galaxy_filename(self, model_name, psf_model_name, noise_model_name, stamp_size, 
                                     model_dico, psf_model_dico, 
                                     noise_model_dico, shear_model_dico,
                                     output_image_filepath, output_catalog_filepath, 
                                     job, worker):
      
      # -- Galaxy Parameter keys that can be adjusted
      relevant_keys = []
      if model_name in ["scsersic", "gscsersic"]:
         relevant_keys.extend([("GAL_FLUX", "F"), ("DISK_INDEX", "dn"), ("DISK_HLR", "dre"),\
                               ("GAL_E1", "e1"), ("GAL_E2", "e2")])
      elif model_name in ["gbdsersic_r", "gbdsersic_d"]:  
         relevant_keys.extend([("GAL_FLUX", "F"), ("DISK_INDEX", "dn"), ("DISK_HLR", "dre"),\
                               ("BULGE_INDEX", "bn"), ("DISK_HLR", "dre"), ("BULGE_HLR", "bre"),\
                               ("GAL_E1", "e1"), ("GAL_E2", "e2")])
      elif model_name in ["gaussian"]:
         relevant_keys.extend([("GAL_FLUX", "F"), ("DISK_INDEX", "dn"), ("DISK_SIGMA", "sig"),\
                               ("GAL_E1", "e1"), ("GAL_E2", "e2")])
         
      psf_relevant_keys = []
      if psf_model_name in ["gaussian", "gmoffat"]:
         psf_relevant_keys.extend([("PSF_FWHM", "PSF_FWHM")])
         
         
      img_file_dir, img_file_name = os.path.split(output_image_filepath)   
      cat_file_dir, cat_file_name = os.path.split(output_catalog_filepath)   

      (img_name, _, _) = img_file_name.split("-")
      (cat_name, _, _) = cat_file_name.split("-")

      # --- Galaxy parameters
      section_name = "GALAXY.MODEL.{0}".format(model_name.upper())
      section_dico = worker.config.get_section_data(section_name)

      name_suffix = "{0}_{1}x{1}".format(model_name, stamp_size)
      for (key, short_key) in relevant_keys:
         if key in section_dico:
            dist_info_dico = eval(section_dico[key])
            if not dist_info_dico is None:
               if  dist_info_dico["dist"].lower() == "constant":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, dist_info_dico["params"][0])
               elif dist_info_dico["dist"].lower() == "sequence":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, 
                                                   job.sequence_value_dico[(section_name, key)])
               else:  
                  name_suffix += "_{0}_{1:.3f}".format(short_key, 
                                                   model_dico[ key])

      # --- PSF parameters
      psf_section_name = "PSF.MODEL.{0}".format(psf_model_name.upper())
      psf_section_dico = worker.config.get_section_data(psf_section_name)

      for (key, short_key) in psf_relevant_keys:
         if key in psf_section_dico:
            psf_dist_info_dico = eval(psf_section_dico[key])
            if not psf_dist_info_dico is None:
               if  psf_dist_info_dico["dist"].lower() == "constant":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, psf_dist_info_dico["params"][0])
               elif psf_dist_info_dico["dist"].lower() == "sequence":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, 
                                                   job.sequence_value_dico[(psf_section_name, key)])
               else:  
                  name_suffix += "_{0}_{1:.3f}".format(short_key, 
                                                   psf_model_dico[ key])

      # -- Shear Parameter
      shear_section_name = "SHEAR"
      shear_section_dico = worker.config.get_section_data(shear_section_name)

      shear_relevant_keys = [("GAL_G1", "g1"), ("GAL_G2", "g2")]
      for (key, short_key) in shear_relevant_keys:
         if key in shear_section_dico:
            dist_info_dico = eval(shear_section_dico[key])
            if not dist_info_dico is None:
               if  dist_info_dico["dist"].lower() == "constant":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, dist_info_dico["params"][0])
               elif dist_info_dico["dist"].lower() == "sequence":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, 
                                                   job.sequence_value_dico[(section_name, key)])
               else:  
                  name_suffix += "_{0}_{1:.3f}".format(short_key, 
                                                   shear_model_dico[ key])
            
      # --- Noise parameters
      if noise_model_name != "no_noise":

         section_name = "NOISE.MODEL.{0}".format(noise_model_name.upper())
         noise_section_dico = worker.config.get_section_data(section_name)      
   
         noise_relevant_keys = []
         if noise_model_name in ["simple_gaussian"]:
            noise_relevant_keys.extend([("SIGMA_NOISE", "sgn"), ("ADDITIONAL_SKY", "addsky"),\
                                        ("WANTED_SNR", "wsnr")])
         elif noise_model_name in ["simple_poisson"]:  
            noise_relevant_keys.extend([("SUBSTRACTED_SKY", "subsky"), \
                                        ("ADDITIONAL_SKY", "addsky"),  ("WANTED_SNR", "wsnr")])
         elif noise_model_name in ["CCD_noise"]:
            noise_relevant_keys.extend([("SUBSTRACTED_SKY", "subsky"), \
                                        ("ADDITIONAL_SKY", "addsky"),  ("WANTED_SNR", "wsnr")])
#                                         ("GAIN", "gain") #, ("READ_NOISE", "readn")])
   
         for (key, short_key) in noise_relevant_keys:
            
            if key in noise_section_dico:
               dist_info_dico = eval(noise_section_dico[key])

               if not dist_info_dico is None:
                  if  dist_info_dico["dist"].lower() == "constant":

                     if key in  ["ADDITIONAL_SKY", "SUBSTRACTED_SKY"]:
                        if dist_info_dico["params"][0] <= 0:
                           continue
                     if key == "WANTED_SNR":
                        if dist_info_dico["params"][0] == -1:
                           continue
                     name_suffix += "_{0}_{1:.2f}".format(short_key, dist_info_dico["params"][0])

                  elif dist_info_dico["dist"].lower() == "sequence":
                     
                     if key in  ["ADDITIONAL_SKY", "SUBSTRACTED_SKY"]:
                        if job.sequence_value_dico[(section_name, key)] <= 0:
                           continue
                     if key == "WANTED_SNR":
                        if job.sequence_value_dico[(section_name, key)] == -1:
                           continue
                     name_suffix += "_{0}_{1:.2f}".format(short_key, 
                                                job.sequence_value_dico[(section_name, key)])
                  else:
                     if dist_info_dico["scope"] == "field":
                        # --- Include the actual randimly-generated value in the filename
                        name_suffix += "_{0}_{1:.2f}".format(short_key, noise_model_dico[key])
                     else:   
                        # --- Include the distribution name and parameters in the filename
                        label = "{0}_{1}".format(
                                    dist_info_dico["dist"], 
                                    "_".join([str(p) for p in dist_info_dico["params"]]))
                        name_suffix += "_{0}_{1}".format(short_key, label) 
                        
      img_file_name = img_file_name.replace(img_name, img_name + "_" + name_suffix)
      output_image_filepath = os.path.join(img_file_dir, img_file_name)

      cat_file_name = cat_file_name.replace(cat_name, cat_name + "_" + name_suffix)
      output_catalog_filepath = os.path.join(cat_file_dir, cat_file_name)
      
      return output_image_filepath, output_catalog_filepath 
   
   
   # -----------------------------------------------------------------------------------------------
   def _adjust_psf_filename(self, model_name, stamp_size, 
                                  model_dico,
                                  output_image_filepath, output_catalog_filepath, 
                                  job, worker):
      
      # -- Keys that can be adjusted
      relevant_keys = []
      if model_name in ["gmoffat"]:
         relevant_keys.extend([("PSF_FLUX", "F"), ("PSF_FWHM", "FWHM"), ("PSF_BETA", "beta"),\
                               ("PSF_E1", "e1"), ("PSF_E2", "e2")])
      elif model_name in ["gaussian"]:  
         relevant_keys.extend([("PSF_FLUX", "F"), ("PSF_FWHM", "FWHM"),\
                               ("PSF_E1", "e1"), ("PSF_E2", "e2")])
         
      img_file_dir, img_file_name = os.path.split(output_image_filepath)   
      cat_file_dir, cat_file_name = os.path.split(output_catalog_filepath)   

      (img_name, _, _) = img_file_name.split("-")
      (cat_name, _, _) = cat_file_name.split("-")

      section_name = "PSF.MODEL.{0}".format(model_name.upper())
      section_dico = worker.config.get_section_data(section_name)

      name_suffix = "{0}_{1}x{1}".format(model_name, stamp_size)
      for (key, short_key) in relevant_keys:
         if key in section_dico:
            dist_info_dico = eval(section_dico[key])
            if not dist_info_dico is None:
               if  dist_info_dico["dist"].lower() == "constant":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, dist_info_dico["params"][0])
               elif dist_info_dico["dist"].lower() == "sequence":
                  name_suffix += "_{0}_{1:.2f}".format(short_key, 
                                                   job.sequence_value_dico[(section_name, key)])
      img_file_name = img_file_name.replace(img_name, img_name + "_" + name_suffix)
      output_image_filepath = os.path.join(img_file_dir, img_file_name)

      cat_file_name = cat_file_name.replace(cat_name, cat_name + "_" + name_suffix)
      output_catalog_filepath = os.path.join(cat_file_dir, cat_file_name)
   
      return output_image_filepath, output_catalog_filepath 