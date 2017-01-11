"""! 
   sp_validate.py - Validation of results by resampling
"""

# -- Python imports
import os, sys
import string
import time
import math
import numpy

# --- Module-specific imports
from sp_helper import *       # helper utility functions
from sp_plot import *         # plotter 


# -------------------------------------------------------------------------------------------------
class Validator(object):
   """! Validation of interpolated results by resampling """

   def __init__(self, master):
      """! Valisator constructor """

      self._helper  = SpredictHelper()              # helper utility functions
      self._plotter = SpredictPlotter()             # plotter

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the GfitHelper instance. """
      return self._helper

   @property
   def plotter(self):
      """! @return the GfitPlotter instance. """
      return self._plotter

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~


   # ----------------------------------------------------------------------------------------------- 
   def cross_validate(self, method, model, param_name, 
                            x_tile_no, y_tile_no, input_coords, input_values, job, worker):
      """!
         Remove outliers from input values @c input_values of input coordinates 
         @c input_coords
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param param_name name of the parameter to interpolate
         @param input_coords input coordinate tuples (x_input_coords, y_input_coords)
         @param input_values model parameter values at positions @c input_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instance   
      """

      validation_dico = {}

      # --- Split the dataset into N folds: N-1 for training, 1 for validation with permutations
      fold_data_dico = self._split_dataset_for_CV(model, param_name, 
                                                  input_coords, input_values, job, worker)

      # --- Perform k-fold cross-validation
      validation_dico = self._cross_validate(method, model, param_name, fold_data_dico, 
                                             job, worker)

      # --- Dump validation data for this tile to a catalog
      self._create_validation_catalog(method, model, param_name, 
                                      x_tile_no, y_tile_no, validation_dico, job, worker)

      # --- Make diagnostic validation plots (if requested)
      self._create_validation_scatter_plot(method, model, param_name, 
                                           x_tile_no, y_tile_no, validation_dico, job, worker)
      self._create_validation_histograms(method, model, param_name, 
                                         x_tile_no, y_tile_no, validation_dico, job, worker)

      # --- Create a file with summary cross-validation statistics per tile
      self._dump_validation_stats(method, model, param_name, 
                                  x_tile_no, y_tile_no, validation_dico, job, worker)    

      return validation_dico
      

   # ~~~~~~~~~~~~~~~
   # Private methods 
   # ~~~~~~~~~~~~~~~

   # ----------------------------------------------------------------------------------------------- 
   def _split_dataset_for_CV(self, model, param_name, input_coords, input_values, job, worker):
      """! 
         Split the dataset into N-1 training sets and 1 validation set 
         @param method object of a type inheriting from the BaseMethod class (e.g. RBF)
         @param model object of a type inheriting from the BaseModel class (e.g. Moffat)
         @param param_name name of the parameter to interpolate
         @param input_coords input coordinate tuples (x_input_coords, y_input_coords)
         @param input_values model parameter values at positions @c input_coords
         @param job SpredictJob object instance   
         @param worker SpredictWorkerMPI or SpredictWorkerSMP worker object instanc    
      """

      fold_data_dico = {}      

      # --- Nb of folds
      CV_section_name = "{0}.CV.{1}".format(model.name.upper(), param_name)
      nb_folds = model.config.get_as_int("FOLD_NUMBER", CV_section_name)

      nb_values = len(input_values)    
      fold_size = int(math.ceil(float(nb_values) / nb_folds))

      # --- make sure the data are randomized 
      shuffled_indice = numpy.arange(0, nb_values)
      numpy.random.shuffle(shuffled_indice)

      #print "nb_folds:", nb_folds
      #print "fold_size:", fold_size
      #print "nb_values:", nb_values
      #print "shuffled_indice:", shuffled_indice
      
      # --- populate the folds
      x_coords, y_coords = zip(*input_coords) 
      x_coords = numpy.asarray(x_coords)
      y_coords = numpy.asarray(y_coords)

      for k in xrange(0, nb_folds):
         #print k*fold_size, (k+1)*fold_size
         fold_indice = shuffled_indice[k*fold_size:(k+1)*fold_size]
         #print k, fold_indice
         fold_data_dico[k] = {}
         fold_data_dico[k]["x_coords"] = x_coords[fold_indice]
         fold_data_dico[k]["y_coords"] = y_coords[fold_indice]
         fold_data_dico[k]["z_values"] = input_values[fold_indice]
         #print k, len(fold_data_dico[k]["z_values"])

      return fold_data_dico


   # ----------------------------------------------------------------------------------------------- 
   def _cross_validate(self, method, model, param_name, fold_data_dico, job, worker):

      validation_dico = {}
      validation_dico["calc_values"] = []
      validation_dico["true_values"] = []
      #validation_dico["x_target_coords"] = []
      #validation_dico["y_target_coords"] = []

      fold_indice = numpy.array(sorted(fold_data_dico.keys()))
      nb_folds = len(fold_indice)
      
      validation_per_fold_dico = {}

      # --- Perform k-fold cross-validation
      for k in xrange(0, nb_folds):

         # --- Select the validation fold and the remaining training folds
         valid_fold_index  = fold_indice[0]
         train_fold_indice = fold_indice[1:]

         #print "fold:", k, "valid:", valid_fold_index, "train:",train_fold_indice    
   
         # --- Cross validate given the above fold partition
         validation_per_fold_dico[k] = self._cross_validate_fold(method, model, param_name, 
                                                                 fold_data_dico,
                                                                 valid_fold_index, 
                                                                 train_fold_indice,
                                                                 job, worker)
         #print validation_per_fold_dico[k]["calc_values"][0:10], type(validation_per_fold_dico[k]["calc_values"][0:10])

         #validation_dico["x_target_coords"].extend(validation_per_fold_dico[k]["x_target_coords"])
         #validation_dico["y_target_coords"].extend(validation_per_fold_dico[k]["y_target_coords"])
         validation_dico["calc_values"].extend(validation_per_fold_dico[k]["calc_values"])
         validation_dico["true_values"].extend(validation_per_fold_dico[k]["true_values"])

         # --- Choose another fold combination
         fold_indice = numpy.roll(fold_indice, 1)

      # --- Compute cross-validation statistics and populate the dictionary with summary stats
      self._compute_cross_validation_stats(method, model, param_name, validation_dico, job, worker)

      #print "### stat:", validation_dico["summary_stats"]
      #print "### r:", len(validation_dico["residuals"]), len(numpy.unique(validation_dico["residuals"]))
      #print "### x:", len(validation_dico["x_target_coords"]), len(numpy.unique(validation_dico["x_target_coords"]))
      #print "### y:", len(validation_dico["y_target_coords"]), len(numpy.unique(validation_dico["y_target_coords"]))

      return validation_dico


   # ----------------------------------------------------------------------------------------------- 
   def _cross_validate_fold(self, method, model, param_name, 
                            fold_data_dico, valid_fold_index, train_fold_indice, job, worker):

      result_dico = {}
      result_dico["calc_values"] = []
      result_dico["true_values"] = []
      #result_dico["x_target_coords"] = []
      #result_dico["y_target_coords"] = []

      #print "*** train_fold_indice:", train_fold_indice
      
      # --- Perform k-fold cross-validation on one combination {1 valid_set, N-1 training sets}   
      for k in train_fold_indice:
         
         # --- Coordinates and values of the curent fold
         x_train_coords = fold_data_dico[k]["x_coords"]
         y_train_coords = fold_data_dico[k]["y_coords"]
         z_train_values = fold_data_dico[k]["z_values"]

         # --- Target coordinates and *true* values to check against
         x_target_coords = fold_data_dico[valid_fold_index]["x_coords"]
         y_target_coords = fold_data_dico[valid_fold_index]["y_coords"]
         z_target_values = fold_data_dico[valid_fold_index]["z_values"]

         # --- Interpolate on fold k data
         interp_dico = method.predict(model, param_name, 
                                      zip(x_train_coords, y_train_coords), 
                                      zip(x_target_coords, y_target_coords),  
                                      z_train_values, 
                                      job, worker)   

         #result_dico["x_target_coords"].extend(x_target_coords)
         #result_dico["y_target_coords"].extend(y_target_coords)
         result_dico["calc_values"].extend(interp_dico["results"]["target_values"])
         result_dico["true_values"].extend(z_target_values)

         #print "train no: ", k, len(interp_dico["results"]["target_values"]), len(z_target_values)

      #print "calc_values:", result_dico["calc_values"][0:10], len(result_dico["calc_values"])    
      #print "true_values:", result_dico["true_values"][0:10], len(result_dico["true_values"])    
      #print "x_coords:", result_dico["true_values"][0:10], len(result_dico["x_target_coords"])    
      #print "y_cords:", result_dico["true_values"][0:10], len(result_dico["y_target_coords"])    

      return result_dico

   # ----------------------------------------------------------------------------------------------- 
   def _compute_cross_validation_stats(self, method, model, param_name, 
                                             validation_dico, job, worker):
      """! Compute cross-validation statistics """

      validation_dico["summary_stats"] = {}

      calc_values = numpy.asarray(validation_dico["calc_values"])
      true_values = numpy.asarray(validation_dico["true_values"])    

      residuals = true_values - calc_values
      abs_residuals = numpy.absolute(residuals)
      residuals_sq = residuals**2

      validation_dico["residuals"] = residuals
      validation_dico["abs_residuals"] = abs_residuals

      validation_dico['summary_stats']['min_residuals'] = numpy.min(residuals)
      validation_dico['summary_stats']['max_residuals'] = numpy.max(residuals)
      validation_dico['summary_stats']['mean_residuals'] = numpy.mean(residuals)
      validation_dico['summary_stats']['median_residuals'] = numpy.median(residuals)

      validation_dico['summary_stats']['min_abs_residuals'] = numpy.min(abs_residuals)
      validation_dico['summary_stats']['max_abs_residuals'] = numpy.max(abs_residuals)
      validation_dico['summary_stats']['mean_abs_residuals'] = numpy.mean(abs_residuals)
      validation_dico['summary_stats']['median_abs_residuals'] = numpy.median(abs_residuals)

      validation_dico['summary_stats']['var_residuals'] = numpy.var(residuals)
      validation_dico['summary_stats']['var_abs_residuals'] = numpy.var(abs_residuals)
      validation_dico['summary_stats']['RMSE'] = math.sqrt(numpy.mean(numpy.mean(residuals_sq)))

#      validation_dico['summary_stats']['sum_residuals'] = numpy.sum(validation_dico['residuals'])
#      validation_dico['summary_stats']['sum_abs_residuals'] = numpy.sum(validation_dico['abs_residuals'])
#      validation_dico['summary_stats']['sum_residuals_sq']  = numpy.sum(validation_dico['residuals_sq'])
#      validation_dico['summary_stats']['RMSE'] = math.sqrt(validation_dico['summary_stats']['MSE'])

#      validation_dico['summary_stats']['min_abs_residuals']    = numpy.min(validation_dico['abs_residuals'])
#      validation_dico['summary_stats']['max_abs_residuals']    = numpy.max(validation_dico['abs_residuals'])
#      validation_dico['summary_stats']['mean_abs_residuals']   = numpy.mean(validation_dico['abs_residuals'])
#      validation_dico['summary_stats']['median_abs_residuals'] = numpy.median(validation_dico['abs_residuals'])

      return validation_dico


   # -----------------------------------------------------------------------------------------------
   def _create_validation_scatter_plot(self, method, model, param_name, x_tile_no, y_tile_no,
                                             validation_dico, job, worker):


      # --- Should we make CV validation plots?
      CV_section_name = "{0}.CV.{1}".format(model.name.upper(), param_name)
      create_CV_plots = model.config.get_as_boolean("CREATE_TILE_VALIDATION_PLOTS", CV_section_name)

      if create_CV_plots:

         # --- Where plots will be stored
         branch_tree = job.get_branch_tree()
         plot_output_dir = os.path.join(worker.plot_output_dir, branch_tree)
         plot_output_dir = os.path.join(plot_output_dir, "cross-validation")
         if not os.path.isdir(plot_output_dir):
            worker.helper.make_dir(plot_output_dir)
      
         # --- Setup plot titles and filenames
         RMSE = validation_dico["summary_stats"]["RMSE"]

         plot_title = "{0} ({1}) - Tile: ({2},{3}) - Validation plot - "\
                      "Image {4:03d}-{5} (RMSE={6:.4e})".format(
                        param_name, method.name, x_tile_no, y_tile_no, job.img_no, job.epoch, RMSE)

         plot_filename = "validation_plot_{0}_{1}_{2:03d}_{3}_tile_{4}_{5}_{6:.3f}.png".format(
                        param_name, method.name, job.img_no, job.epoch, x_tile_no, y_tile_no, RMSE)

         # --- Compute plot limits
         min_value = min(numpy.min(validation_dico['calc_values']), 
                         numpy.min(validation_dico['true_values']))
         max_value = max(numpy.max(validation_dico['calc_values']), 
                         numpy.max(validation_dico['true_values']))

         # --- Make the plot 
         plotter = SpredictPlotter()
         plotter.plot_scatter(validation_dico, ['calc_values', 'true_values'],
                              plot_filename, plot_output_dir, plot_title=plot_title, pfilter=None,
                              #pxmin_lim=min_value, pxmax_lim=max_value, 
                              #pymin_lim=min_value, pymax_lim=max_value, 
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                              x_label='Predicted values (red)', y_label='True values (green)',
                              pcolor=['#ff4500','#7fff00'], 
                              pmarker=[3,3], pmarkersize=[8,8], is_normed=False, 
                              has_equalaxes=True, has_grid=True, pfigure=None)


   # -----------------------------------------------------------------------------------------------
   def _create_validation_histograms(self, method, model, param_name, x_tile_no, y_tile_no,
                                           validation_dico, job, worker):


      # --- Should we make CV histogram plots?
      CV_section_name = "{0}.CV.{1}".format(model.name.upper(), param_name)
      create_hist_plots = model.config.get_as_boolean("CREATE_HISTOGRAM_PLOTS", CV_section_name)

      if create_hist_plots:

         # --- Where plots will be stored
         branch_tree = job.get_branch_tree()
         plot_output_dir = os.path.join(worker.plot_output_dir, branch_tree)
         plot_output_dir = os.path.join(plot_output_dir, "cv_residuals")
         if not os.path.isdir(plot_output_dir):
            worker.helper.make_dir(plot_output_dir)
      
         RMSE = validation_dico["summary_stats"]["RMSE"]

         # --- Create residuals histogram plots
         plot_title = "{0} ({1}) - Tile: ({2},{3}) - Residuals - "\
                      "Image {4:03d}-{5} (RMSE={6:.4e})".format(
                        param_name, method.name, x_tile_no, y_tile_no, job.img_no, job.epoch, RMSE)

         plot_filename = "hist_residuals_{0}_{1}_{2:03d}_{3}_tile_{4}_{5}_{6:.3f}.png".format(
                        param_name, method.name, job.img_no, job.epoch, x_tile_no, y_tile_no, RMSE)

#         print validation_dico["residuals"]
#         data_dico = {}
#         data_dico["abs_residuals"] = numpy.absolute(validation_dico["residuals"])

         plotter = SpredictPlotter()
         plotter.plot_histogram(validation_dico, ["residuals"],
                                plot_filename, plot_output_dir, plot_title=plot_title, pfilter=None,
                                pnb_bins=20.0, pmax_xlim=None, pcolor='#6a5acd', is_normed=False, 
                                has_equalaxes=False, has_grid=True, pfigure=None)

   # -----------------------------------------------------------------------------------------------
   def _create_validation_catalog(self, method, model, param_name, 
                                        x_tile_no, y_tile_no, validation_dico, job, worker):

      # --- Should we create the catalogs?
      CV_section_name = "{0}.CV.{1}".format(model.name.upper(), param_name)
      create_catalogs = model.config.get_as_boolean("CREATE_TILE_VALIDATION_CATALOGS", 
                                                    CV_section_name)

      if create_catalogs:

         # --- Create a result catalog, marking failed ellipticity values with a constant >= 10
         output_filename = "validation_{0}_{1}_{2:03d}_{3}_tile_{4}_{5}.txt".format(
                           param_name, method.name, job.img_no, job.epoch, x_tile_no, y_tile_no)

         output_directory = os.path.join(worker.result_output_dir, job.get_branch_tree())
         output_directory = os.path.join(output_directory, "cross-validation")
         if not os.path.isdir(output_directory):
            worker.helper.make_dir(output_directory)

         col_list = ["calc_values", "true_values", "residuals", "abs_residuals"]
         col_key_map = {"calc_values":0, "true_values":1, "residuals":2, "abs_residuals":3}
         col_fmt_map = {"calc_values":"%9.6e", "true_values":"%9.6e", \
                        "residuals":"%9.6f", "abs_residuals":"%9.6f"}

         # --- Create the catalog
         self.helper.save_from_list_dico(validation_dico, output_directory, output_filename, 
                                         col_list, 
                                         key_index_map=col_key_map, key_fmt_map=col_fmt_map)

   # -----------------------------------------------------------------------------------------------
   def _dump_validation_stats(self, method, model, param_name, 
                              x_tile_no, y_tile_no, validation_dico, job, worker):
         
      # --- Should we dump CV statistics ?
      CV_section_name = "{0}.CV.{1}".format(model.name.upper(), param_name)
      dump_stats = model.config.get_as_boolean("DUMP_TILE_VALIDATION_STATS", CV_section_name)

      if dump_stats:

         # --- Dump CV statistics to a file
         output_filename = "stats_{0}_{1}_{2}_{3:03d}_{4}.txt".format(
                           param_name, model.name, method.name, job.img_no, job.epoch)
         output_directory = os.path.join(worker.stat_output_dir, job.get_branch_tree())
         if not os.path.isdir(output_directory):
            worker.helper.make_dir(output_directory)
         output_filepath = os.path.join(output_directory, output_filename)
            
         if worker.helper.file_exists(output_filepath):
            # --- Part of the file has already been generated, append statistics to that file
            fd = open(output_filepath, "a")
         else:
            # --- Stats file not yet created, do it now
            fd = open(output_filepath, "w")
            fd.write('\n')
            title = "CROSS-VALIDATION STATISTICS - /{0}/image-{1:03d}-{2:1d} - "\
                      "{3} ({4})\n".format(
                         job.get_branch_tree(), job.img_no, job.epoch, param_name, method.name)
            fd.write(title)
            fd.write("=" * len(title) + "\n\n")
             
         fd.write('\n--------------- Tile ({0},{1}) ---------------\n'.format(x_tile_no, y_tile_no))

         for stat in sorted(validation_dico['summary_stats'].keys()):
            
            fd.write("{0:<20} = {1:9.6e}\n".format(stat, validation_dico['summary_stats'][stat]))

         fd.close()   

#      validation_dico['summary_stats']['min_residuals'] = numpy.min(residuals)
#      validation_dico['summary_stats']['max_residuals'] = numpy.max(residuals)
#      validation_dico['summary_stats']['mean_residuals'] = numpy.mean(residuals)
#      validation_dico['summary_stats']['median_residuals'] = numpy.median(residuals)

#      validation_dico['summary_stats']['min_abs_residuals'] = numpy.min(abs_residuals)
#      validation_dico['summary_stats']['max_abs_residuals'] = numpy.max(abs_residuals)
#      validation_dico['summary_stats']['mean_abs_residuals'] = numpy.mean(abs_residuals)
#      validation_dico['summary_stats']['median_abs_residuals'] = numpy.median(abs_residuals)

#      validation_dico['summary_stats']['var_residuals'] = numpy.var(residuals)
#      validation_dico['summary_stats']['var_abs_residuals'] = numpy.var(abs_residuals)
#      validation_dico['summary_stats']['RMSE'] = math.sqrt(numpy.mean(numpy.mean(residuals_sq)))


# -- EOF sp_validate.p
