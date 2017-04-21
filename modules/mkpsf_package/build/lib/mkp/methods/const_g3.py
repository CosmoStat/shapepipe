"""! 
   const_g3.py - PSF Reconstruction of a GREAT3 Constant PSF
"""

# -- Python imports
import os, sys
import time
import numpy
import math
import scipy.ndimage
import scipy.ndimage.interpolation
import pyfits

# -- External imports
from sconfig import *      # configuration file management

# --- Module-specific imports
from mkp_plot import *           # plotter 
from mkp_help import *           # helper utility functions

# -------------------------------------------------------------------------------------------------
class PSFBuildMethod(object):
   
   """! 
      PSFBuildMethod: method for reconstructing a PSF field
   """

   def __init__(self, method_name, method_dir):
      """! Construct a PSFBuildMethod class """
   
      self._name = method_name   # reconstruction method name -> xxx.cfg, xxx.py 

      # --- Get relevant configuration info
      self._config = self._load_method_config(method_dir)   # method-specific configuration

      self._helper  = MkpHelper()         # helper utility functions
      self._plotter = MkpPlotter()        # plotter

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   # --- Getters

   @property
   def name(self):
      """! @return the name of the build method """   
      return self._name

   @property
   def config(self):
      """! @return the configuration object of the build method """   
      return self._config

   @property
   def helper(self):
      """! @return the MkpHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the Plotter instance. """
      return self._plotter


   # ~~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def make_PSF_field(self, image_file_type, catalog_file_type, job, worker):
      """! Create a PSF field from the GREAT3 constant PSF field (9 PSFs) """

      result_dico = {}

      # --- Generate a PSF .FITS stamps based on te requested configuration
      build_method = self._config.get_as_string("BUILD_METHOD", "PSF_BUILD") 
      result_dico = self._generate_PSF_stamp(build_method, 
                                             image_file_type, catalog_file_type, job, worker)

      return result_dico
   
   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def _load_method_config(self, base_method_dir):
      """! Load the method configuration file """
      
      # Config name is derived from the method name 
      return SConfig(os.path.join(base_method_dir, self.name) + '.cfg')  

   # -----------------------------------------------------------------------------------------------  
   def _generate_PSF_stamp(self, build_method_name, 
                                 image_file_type, catalog_file_type, job, worker):
      """! Create the PSF postage stamp """

      result_dico = {}  # dictionary to populate

      # -- Source GREAT3 PSF image
      image_filepath  = job.get_file_path(image_file_type)
      _, image_filemain = os.path.split(image_filepath)
      stamp_size = job.get_stamp_size()

      # --- Output directory
      output_dir = os.path.join(worker.result_output_dir, job.get_branch_tree())

      # --- Extract stars as .FITS if requested
      if self.config.get_as_boolean("EXTRACT_STAMPS", "DEBUGGING"):
         float_size = self.config.get_as_string("EXTRACTED_STAMP_FLOAT_TYPE", "DEBUGGING")
         center_star =  self.config.get_as_boolean("EXTRACT_STAMP_CENTERED", "DEBUGGING")
         stamp_output_dir = os.path.join(output_dir, "psf_stamps")
         if not os.path.isdir(stamp_output_dir):
            self.helper.make_dir(stamp_output_dir)
         self.helper._extract_stamps(image_filepath, job.img_no, job.epoch, 
                                     image_file_type, stamp_size, stamp_output_dir, 
                                     float_size, center_star)

      # --- Plot star postage stamps if requested
      if self.config.get_as_boolean("PLOT_STAMPS", "DEBUGGING"):
         stamp_output_dir = os.path.join(worker.plot_output_dir, job.get_branch_tree())
         stamp_output_dir = os.path.join(stamp_output_dir, "psf_stamps")
         if not os.path.isdir(stamp_output_dir):
            self.helper.make_dir(stamp_output_dir)
         self.helper._plot_stamps(self.plotter, image_filepath, job.img_no, job.epoch, 
                                  image_file_type, stamp_size, stamp_output_dir)

      # --- If requested, compute dithering offsets wrt. the 1st star stamp in PSF image
      if self.config.get_as_boolean("DITHERING_OFFSETS", "DEBUGGING"):
         self._compute_dithering_offsets(image_filepath, stamp_size, output_dir, job, worker)

      # --- Create a PSF postage stamp that will populate the entire PSF field
      build_func_ptr = eval("self._make_psf_stamp_using_{0}".format(build_method_name))
      psf_stamp = build_func_ptr(build_method_name, 
                                 image_filepath, stamp_size, image_file_type, job, worker)

      # --- Store the stamp at the appropriate location (preserving the original header)
      psf_output_filename = image_filemain   # we keep the same file name but with a different path
      #psf_output_filename = image_filemain.replace("starfield", "psffield")
      psf_output_filepath = os.path.join(output_dir, psf_output_filename)

      image_header = pyfits.getheader(image_filepath)
      stamp_header = image_header
      self.helper.write_as_fits(psf_stamp, psf_output_filepath, header=stamp_header)

      result_dico["psf_output_filepath"] = psf_output_filepath
      result_dico["psf_image"] = psf_stamp
      result_dico["psf_stamp_size"] = stamp_size 
      result_dico["psf_catalog_file_type"] = catalog_file_type
      result_dico["psf_image_file_type"] = image_file_type

      return result_dico

   # -----------------------------------------------------------------------------------------------
   def _make_psf_stamp_using_STACK(self, build_method_name, image_filepath, stamp_size, file_type, 
                                         job, worker):

      # --- Look for the PSF among the 9 that is centered
      image_data = pyfits.getdata(image_filepath)
      (image_width, image_height) = image_data.shape
      float_size = self.config.get_as_string("EXTRACTED_STAMP_FLOAT_TYPE", "DEBUGGING")

      # --- Find centroid of centered stamp
      ref_x_pixel = ref_y_pixel = (stamp_size / 2.0 - 1.0) 
      (x0, y0) = self._find_centered_star(image_data, stamp_size, (ref_x_pixel, ref_y_pixel))
      centered_stamp = image_data[y0:y0+stamp_size, x0:x0+stamp_size]
      (yc0, xc0) = scipy.ndimage.center_of_mass(centered_stamp)

      # --- Stack all stars at the centered position (xc0, yc0)
      stacked_stamp = numpy.zeros((stamp_size, stamp_size))
      nbstamp = 0
      for x in xrange(0, image_width, stamp_size):
         for y in xrange(0, image_height, stamp_size):
            stamp = image_data[y:y+stamp_size, x:x+stamp_size].astype(float_size) # (row, col) np.

            (yc, xc) = scipy.ndimage.center_of_mass(stamp)
            corr = (yc0 - yc, xc0 - xc)

            stamp = scipy.ndimage.interpolation.shift(stamp, corr)

            stacked_stamp += stamp

            #### CHECK ### 
            #(nyc, nxc) = scipy.ndimage.center_of_mass(stamp)
            #print (x,y), (xc0, yc0), (xc, yc), "=>", (nxc, nyc), "corr:", corr, numpy.sum(stamp)
            #### CHECK ### 

            nbstamp += 1

      stacked_stamp /= nbstamp

      #### CHECK ### 
      #(nyc, nxc) = scipy.ndimage.center_of_mass(stacked_stamp)
      #print "Stacked PSF centroid", (nxc, nyc), numpy.sum(stacked_stamp)
      #### CHECK ### 

      return stacked_stamp



   # -----------------------------------------------------------------------------------------------
   def _make_psf_stamp_using_PICK(self, build_method_name, image_filepath, stamp_size, file_type, 
                                        job, worker):

      # TODO: write the full FITS header 

      # --- Look for the PSF among the 9 that is centered
      image_data = pyfits.getdata(image_filepath)
      (image_width, image_height) = image_data.shape
      centering = self.config.get_as_string("CENTERING", build_method_name)
      ref_x_pixel = ref_y_pixel = (stamp_size / 2.0 - 1.0) 
   
      (x, y) = self._find_centered_star(image_data, stamp_size, (ref_x_pixel, ref_y_pixel))

      stamp = image_data[y:y+stamp_size, x:x+stamp_size]

      if centering in ["IF_NEEDED", "ALWAYS"]:

         (yc, xc) = scipy.ndimage.center_of_mass(stamp)       

         if worker.logging_enabled():
            worker.logger.log_info_p(
               "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - Initial centroid ({5:.3f},{6:.3f})".format(
               worker.name, job.get_branch_tree(), job.img_no, job.epoch,  file_type, xc, yc) )

         if centering == "ALWAYS":
            # --- ALWAYS
            corr = (ref_y_pixel - yc, ref_x_pixel - xc)
            if worker.logging_enabled():
               worker.logger.log_info_p(
                 "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - Applying ALWAYS correction "\
                 "({5:.3f},{6:.3f})".format(
                  worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                               corr[0], corr [1]) )
            stamp = scipy.ndimage.interpolation.shift(stamp, corr)
         else: 
            # --- IF_NEEDED
            (xa, ya) = (math.floor(xc + 0.5), math.floor(yc + 0.5))  # actual alignment
            center_pixels = (ref_x_pixel - xa, ref_y_pixel - ya)

            if xa - ref_x_pixel != 0 or ya - ref_y_pixel != 0:
               corr = [0.0, 0.0]
               if xa - ref_x_pixel != 0:
                  corr[1] = ref_x_pixel - xc
               if ya - ref_y_pixel != 0:
                  corr[0] = ref_y_pixel - yc

               # Apply centroid correction
               stamp = scipy.ndimage.interpolation.shift(stamp, corr)

               if worker.logging_enabled():
                  worker.logger.log_info_p(
                    "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - Applying IF_NEEDED correction "\
                    "({5:.3f},{6:.3f})".format(
                     worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type, 
                                  corr[1], corr [0]) )

               # Checking:
               (new_yc, new_xc) = scipy.ndimage.center_of_mass(stamp)    
               (xa, ya) = (math.floor(new_xc + 0.5), math.floor(new_yc + 0.5))
               if (xa, ya) != (ref_x_pixel, ref_y_pixel):
                  self.helper.print_warning("Misaligned star at {0}".format((xa, ya)))
            else:
               if worker.logging_enabled():
                     worker.logger.log_info_p(
                       "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - NO correction required ".format(
                       worker.name, job.get_branch_tree(), job.img_no, job.epoch, file_type) )                        

         if worker.logging_enabled():
            (new_yc, new_xc) = scipy.ndimage.center_of_mass(stamp)
            worker.logger.log_info_p(
                "{0} - /{1}/image-{2:03d}-{3:1d} - {4} - Final centroid "\
                "({5:.3f},{6:.3f})".format(
                worker.name, 
                job.get_branch_tree(), job.img_no, job.epoch, file_type, new_xc, new_yc) )
         
      return stamp

   # -----------------------------------------------------------------------------------------------
   def _make_psf_stamp_using_PICK_CENTERED(self, build_method_name, image_filepath, 
                                                 stamp_size, file_type, job, worker):

      # TODO: write the full FITS header 

      # --- Look for the PSF among the 9 that is centered
      image_data = pyfits.getdata(image_filepath)
      (image_width, image_height) = image_data.shape
      ref_x_pixel = ref_y_pixel = (stamp_size-1)/2.0

      #print "###", job.get_branch_tree(), "Img {0:03d}".format(job.img_no),

      (x, y) = self._find_best_centered_star(image_data, stamp_size, (ref_x_pixel, ref_y_pixel))

      #print "=>", (x, y)

      return image_data[y:y+stamp_size, x:x+stamp_size]  

   # -----------------------------------------------------------------------------------------------
   def _find_centered_star(self, image_data, stamp_size, ref_pixels):

      found_coords = (0,0)  # by default, return the coordinates of the first star in the image

#      # --- Collect postage stamps in PSF image
#      found_stamp_dico = {}
#      (image_width, image_height) = image_data.shape
#      for x in xrange(0, image_width, stamp_size):
#         for y in xrange(0, image_height, stamp_size):
#            stamp = image_data[y:y+stamp_size, x:x+stamp_size]
#            (yc, xc) = scipy.ndimage.center_of_mass(stamp)
#            (xa, ya) = (math.floor(xc + 0.5), math.floor(yc + 0.5))  # actual alignment

#            # Found candidate stamp
#            if (xa, ya) == ref_pixels:
#               found_stamp_dico[(x,y)] = stamp

#      # --- Among the candidate stamps, find that with the minimum required shift
#      shift_dico = {}
#      for (x, y) in found_stamp_dico.keys():
#         (yc, xc) = scipy.ndimage.center_of_mass(found_stamp_dico[(x,y)])
#         shift_dico[(x,y)] = math.fabs(xc-ref_pixels[0]) + math.fabs(yc-ref_pixels[1]) 

#      if len(shift_dico) > 0:
#         found_coords = sorted(shift_dico.items(), key=itemgetter(1))[0][0]

      return found_coords

   # -----------------------------------------------------------------------------------------------
   def _find_best_centered_star(self, image_data, stamp_size, ref_pixels):
      """! 
         Find the best-centered stamp, i.e. whose coordinates are closest to (23.5, 23.5) or
         (47,5, 47.5) depending on the size of the stamp
      """

      found_coords = (0,0)  # by default, return the coordinates of the first star in the image

#      # --- Among the stamps, find that with the minimum required shift
#      shift_dico = {}
#      (image_width, image_height) = image_data.shape
#      for x in xrange(0, image_width, stamp_size):
#         for y in xrange(0, image_height, stamp_size):
#            stamp = image_data[y:y+stamp_size, x:x+stamp_size]
#   
#            (yc, xc) = scipy.ndimage.center_of_mass(stamp)
#            print "(x,y)", (x,y), "(xc,yc):", (xc, yc),
#            shift_dico[(x,y)] = math.fabs(xc-ref_pixels[0]) + math.fabs(yc-ref_pixels[1])

#            corr = (ref_pixels[1]-yc, ref_pixels[0]-xc)
#            print (corr[1], corr[0]), shift_dico[(x,y)]

      return found_coords

   # -----------------------------------------------------------------------------------------------
   def _compute_dithering_offsets(self, image_filepath, stamp_size, output_dir, job, worker):
      """!  
          Compute dithering offsets compared to the first star stamp in the *starfield-xxx-x.fits
          image, which is ascribed an offset of (dx=0.0, dy=0.0)
      """

      # --- Take first stamp as reference stamp (dx=0, dy=0)
      image_data = pyfits.getdata(image_filepath)
      ref_stamp = image_data[0:0+stamp_size, 0:0+stamp_size]
      ref_yc, ref_xc = scipy.ndimage.center_of_mass(ref_stamp)

      # --- Compute and store dithering offsets compared to the reference stamp
      output_filepath = image_filepath.replace("starfield_image", "dither_offsets")
      output_filepath = output_filepath.replace(".fits", ".txt")
      dirname, filename = os.path.split(output_filepath)
      output_dir = os.path.join(output_dir, "psf_stamps")
      #if os.path.isdir(output_dir):
      self.helper.make_dir(output_dir)

      output_filepath = os.path.join(output_dir, filename)
      fd = open(output_filepath, "w")
      fd.write("# xdither_pixels ydither_pixels\n")

      (image_width, image_height) = image_data.shape
      for x in xrange(0, image_width, stamp_size):
         for y in xrange(0, image_height, stamp_size):
            stamp = image_data[y:y+stamp_size, x:x+stamp_size]
            (yc, xc) = scipy.ndimage.center_of_mass(stamp)
            (dx, dy) = (xc-ref_xc, yc-ref_yc) 
            fd.write("{0:.12f} {1:.12f}\n".format(dx, dy))


      fd.close()

   # -----------------------------------------------------------------------------------------------
   class BuildMethodError(Exception):
      """! 
         Exception thrown related to the PSF build method
      """
      def __init__(self, msg):
         """!
            @param msg error message 
         """
         self._msg = msg
         self._exc_info = sys.exc_info()[1]

      def __str__(self):
         """! String representation of the BuildMethodError class """
         if self._exc_info is not None:
            return "MKPSF *** ERROR ***: {0} ({1})".format(self._msg, self._exc_info)
         else:
            return "MKPSF *** ERROR ***: {0}".format(self._msg)

