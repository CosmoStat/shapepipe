"""! 
   Plotting functions.
"""

# -- Python imports
import os, sys
import math
import numpy
import imp
from copy import deepcopy
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as pylab3
from matplotlib.font_manager import FontProperties
from matplotlib import rc

# -- External imports

# --- Module-specific imports
from gfit_helper import *        # helper utility functions


# -------------------------------------------------------------------------------------------------
class GfitPlotter(object):
   
   """! 
      Class with a number of convenient plotting methods.   
   """

   def __init__(self):
      """! GfitPlotter constructor """

      self._helper = GfitHelper()      # helper

   # ~~~~~~~~~~
   # Properties 
   # ~~~~~~~~~~

   @property
   def helper(self):
      """! @return the GfitHelper instance. """
      return self._helper

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def plot_stamp(self, stamp, plot_title='Stamp', cmap="jet", 
                        output_dir='.', output_file='stamp_plot', 
                        color_bar=False, show=False, logger=None):

      """! Plot a postage stamp 
         @param stamp postage stamp to plot
         @param plot_title optional plot title
         @param cmap color map
         @param color_bar if True, add a color bar   
         @param output_dir optional output directory  
         @param output_file optional output file name
         @param show True if the plot is shown once drawn, False otherwise  
         @param logger a FileLogger instance or None
      """

      try:
         height, width = stamp.shape

         figure = matplotlib.figure.Figure()
         figtitle = plot_title
         plt.gcf().text(0.45, 0.93, plot_title, 
                           horizontalalignment='center', 
                           fontproperties=FontProperties(size=11, weight='bold'))
         plt.grid(False)
         axes = figure.gca() 
         axes.set_aspect('equal')
         plt.xlim(1, width)
         plt.ylim(1, height)
         
         x = numpy.arange(0, width+1)   
         y = numpy.arange(0, height+1)   
         plot = plt.pcolor(x, y, stamp, shading='flat', cmap=cmap)
         if color_bar:
            plt.colorbar(plot, shrink=1.0)
            #plt.clim(-1.0, 1.0) 
         
         plt.savefig(os.path.join(output_dir, output_file))

         if show:
            plt.show()

         plt.clf()
         plt.close()

      except:    
         if not logger is None:    
            logger.log_error_p('Exception thrown during plotting ({0})'.format(
                                                                               sys.exc_info()[1]))
         else:
            print("Exception thrown during plotting ({0})".format(sys.exc_info()[1]))     

      finally:
         plt.clf()
         plt.close()


   # -----------------------------------------------------------------------------------------------
   def plot_stamp_3D(self, stamp, plot_title='Stamp 3D', 
                           x_label=None, y_label=None, z_label=None, has_grid=False, cmap="jet",
                           rstride=1, cstride=1,   
                           output_dir='.', output_file='stamp_plot_3D', show=False, logger=None):
      """! 
         plot a stamp in 3D 
         @param stamp postage stamp to plot
         @param plot_title optional plot title   
         @param x_label labbel string for the x axis
         @param y_label labbel string for the y axis
         @param z_label labbel string for the z axis
         @param cmap color map
         @param rstride row stride
         @param cstride column stride
         @param has_grid if True, draw a grid 
         @param output_dir optional output directory  
         @param output_file optional output file name
         @param show True if the plot is shown once drawn, False otherwise  
         @param logger a FileLogger instance or None            
      """

      try:
     
         stamp_size, stamp_size = stamp.shape

         figure = plt.figure()
         figtitle = plot_title
         plt.gcf().text(0.5, 0.93, plot_title, horizontalalignment='center', 
                                     fontproperties=FontProperties(size=11, weight='bold'))
         plt.grid(has_grid)
      #   axes = figure.gca() 
      #   axes.set_aspect('equal')
         
         lin_space = numpy.linspace(0, stamp_size -1, int(stamp_size))
         x, y = numpy.meshgrid(lin_space, lin_space)
         z = stamp

         #ax = matplotlib.axes3d.Axes3D(figure)
         ax = pylab3.Axes3D(figure)
         #surface = ax.plot_surface(x,y,z, rstride=1, cstride=1, cmap='gist_heat')
         surface = ax.plot_surface(x,y,z, rstride=rstride, cstride=cstride, cmap=cmap)

         if not x_label is None:
            ax.set_xlabel(x_label, size='medium', weight='heavy', color='red')
         if not y_label is None:
            ax.set_ylabel(y_label, size='medium', weight='heavy', color='red')
         if not z_label is None:
            ax.set_zlabel(z_label, size='medium', weight='heavy', color='red')

         plt.savefig(os.path.join(output_dir, output_file))

         if show:
            plt.show()
   
      except:        
         if not logger is None:    
            logger.log_error_p('Exception thrown during plotting ({0})'.format(
                                                                           sys.exc_info()[1]))
         else:
            print("Exception thrown during plotting ({0})".format(sys.exc_info()[1]))     

      finally:
         plt.clf()
         plt.close()


   # -----------------------------------------------------------------------------------------------
   def plot_histogram(self, data_dico, value_names,
                            plot_filename, output_directory, plot_title, pfilter=None,
                            pnb_bins=10, pmin_xlim=None, pmax_xlim=None, 
                            pcolor='b', pfacecolor="#F5F5F5", is_normed=False, 
                            px_sci=False, py_sci=False, pminor=True,
                            has_equalaxes=False, has_grid=True, pfigure=None, logger=None):

      """! Plot histograms of possibly multiple arrays of floats found in a dictionary """

      try:

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure

         rect = fig.patch
         rect.set_facecolor(pfacecolor)
                            
         axes = fig.gca() 
       
         iplot = 1

         figtitle = plot_title

         if pfigure == None:
            fig.suptitle(figtitle, horizontalalignment='center', position=(0.5, 0.95),
                                   fontproperties=FontProperties(size=11))

         for name in value_names:

            if value_names[iplot-1] in data_dico:

               if len(value_names) > 1:
                  if len(value_names) < 3:         
                     plt.subplot(2, 1, iplot)
                  elif len(value_names) >= 3 and len(value_names) < 5:
                     plt.subplot(2, 2, iplot)
                  else:
                     plt.subplot(3, 2, iplot)

               plt.grid(has_grid)
               if has_equalaxes:
                  axes.set_aspect('equal')

               values = numpy.array(data_dico[value_names[iplot-1]], dtype=numpy.float64)
               values = numpy.nan_to_num(values)

               if len(values) > 0:
                  if pfilter != None:
                     filtered_values = values[pfilter]

                     if len(filtered_values) > 0:
                        values = filtered_values
                  color = pcolor
                  x_label = value_names[iplot-1]
                  y_label = ''

                  plt.xlabel(x_label, fontproperties=FontProperties(size=10, weight='regular'))
                  if pmin_xlim is None:
                     pmin_xlim = numpy.min(values) #FCS ADDED NUMPY.

                  if pmax_xlim is not None:
                     pmax_xlim = numpy.min(pmax_xlim, numpy.max(values))  #FCS ADDED NUMPY.
                  else:
                     pmax_xlim = numpy.max(values)

                  if pmax_xlim > pmin_xlim:
                     plt.xlim(pmin_xlim, pmax_xlim)

                  plt.hist(values, pnb_bins, normed=is_normed, alpha=1.0,
                           edgecolor="k", facecolor=pcolor, histtype='bar')

               iplot += 1
            
         plt.setp(axes.get_xticklabels(), rotation='horizontal', fontsize=10)
         plt.setp(axes.get_yticklabels(), rotation='horizontal', fontsize=10)

         if px_sci:
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0), labelsize=10)
         else:
            plt.ticklabel_format(style='plain', axis='x', scilimits=(0,0), labelsize=10)

         if py_sci:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), labelsize=10)
         else:
            plt.ticklabel_format(style='plain', axis='x', scilimits=(0,0), labelsize=10)

         if pminor:
            plt.minorticks_on()

         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename), 
                                     facecolor=fig.get_facecolor(), edgecolor="k", 
                                     transparent=True)
            plt.clf()
            plt.close()

      except (OverflowError, ValueError, IndexError):       
         if not logger is None:    
            logger.log_error_p('Exception thrown during plotting ({0})'.format(sys.exc_info()[1]))
         else:
            print("Exception thrown during plotting ({0})".format(sys.exc_info()[1])) 
      except Exception as detail:
            print("Exception thrown during plotting ({0})".format(detail) )


   # -----------------------------------------------------------------------------------------------
   def plot_scatter(self, data_dico, value_names,
                    plot_filename, output_directory, plot_title=None, pfilter=None,
                    pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                    x_label=None, y_label=None, pnorm=None, pcmap='jet',
                    pcolor='b', pmarker='o', pmarkersize=5, pfacecolor="#F5F5F5",
                    pmarkerfacecolor=None, pmarkeredgecolor=None, 
                    px_sci=False, py_sci=False, pminor=False, pstats=False,
                    has_lstsqfit=False, pfitdegree=2, pfitcolor='red', pfitmarker='',
                    pfitmarkersize=1.0, pfitmarkerfacecolor='black', 
                    pfitlinewidth=1, pfitlinestyle='-',
                    is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None):

      """! Make a scatter plot """

      try:

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure
         iplot = 1

         axes = fig.gca() 

         rect = fig.patch
         rect.set_facecolor(pfacecolor)

         if not plot_title is None:
            figtitle = plot_title
         else:
            figtitle = 'Scatter plot of %s' %(value_names)

         fig.suptitle(figtitle, horizontalalignment='center', position=(0.5, 0.95),
                                fontproperties=FontProperties(size=11))

         plt.grid(has_grid)
         if has_equalaxes:
            axes.set_aspect('equal')

         # Two values to represent
         if not value_names[0] in data_dico or not value_names[1] in data_dico:
            return
         values_x = numpy.asarray(data_dico[value_names[0]])
         values_y = numpy.asarray(data_dico[value_names[1]])

         if pfilter is not None:

            filtered_values_x = values_x[pfilter]
            filtered_values_y = values_y[pfilter]

            if len(filtered_values_x) > 0:
               values_x = filtered_values_x
            if len(filtered_values_y) > 0:
               values_y = filtered_values_y

         if pxmin_lim is None:
            pxmin_lim = numpy.min(values_x)
         if pxmax_lim is None:
            pxmax_lim = numpy.max(values_x)
         if pymin_lim is None:
            pymin_lim = numpy.min(values_y)
         if pymax_lim is None:
            pymax_lim = numpy.max(values_y)

         if px_sci:
            plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0), labelsize=10)

         if py_sci:
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), labelsize=10)

         if pminor:
            axes.xaxis.grid(True, 'major', ls='solid', lw=0.5, color='gray')
            axes.yaxis.grid(True, 'major', ls='solid', lw=0.5, color='gray')
            axes.xaxis.grid(True, 'minor', ls='solid', lw=0.1, color='gray')
            axes.yaxis.grid(True, 'minor', ls='solid', lw=0.1, color='gray')

            plt.minorticks_on()

         plt.setp(axes.get_xticklabels(), rotation='horizontal', fontsize=10)
         plt.setp(axes.get_yticklabels(), rotation='horizontal', fontsize=10)

         plot = plt.scatter(values_x, values_y, s=pmarkersize, c=pcolor, marker=pmarker,
                            norm=pnorm, cmap=pcmap, facecolor=pmarkerfacecolor, 
                            edgecolor=pmarkeredgecolor)

         if x_label is None:
            x_label = value_names[0]
         if y_label is None:
            y_label = value_names[1]

         plt.xlabel(x_label, fontproperties=FontProperties(size=10, weight='regular'))
         plt.ylabel(y_label, fontproperties=FontProperties(size=10, weight='regular'))

         if not has_equalaxes:
            pxmin = pxmin_lim
            pymin = pymin_lim
            pxmax = pxmax_lim
            pymax = pymax_lim
         else:
            pxmin = pymin = min(pxmin_lim, pymin_lim) 
            pxmax = pymax = max(pxmax_lim, pymax_lim) 

         if pxmax > pxmin:
            plt.xlim(pxmin, pxmax)
         if pxmax > pymin:
           plt.ylim(pymin, pymax)

         # Also include basic stats requested
         if pstats:
            if type(values_y) == list == list:
               values_y = numpy.asarray([values_y])

            if type(values_y) == numpy.ndarray:

               var_min  = numpy.min(values_y)
               var_max  = numpy.max(values_y)
               var_mean = numpy.mean(values_y)
               var_std  = numpy.std(values_y)

               mean_std_label = "Mean / Std: {0:+.3e} / {1:+.3e}".format(var_mean, var_std)
               min_max_label = "Min / Max: {0:+.3e} / {1:+.3e}".format(var_min, var_max)

               plt.annotate(mean_std_label, (0.40, 0.93), xycoords='axes fraction', fontsize=10,  
                            bbox={'facecolor':'#FFFFFF', 'alpha':0.5, 'pad':10})
               plt.annotate(min_max_label, (0.40, 0.86), xycoords='axes fraction', fontsize=10, 
                            bbox={'facecolor':'#FFFFFF', 'alpha':0.5, 'pad':10})


         # Also include a least-squares fit if requested
         if has_lstsqfit: 
            # Must have at least 20 points
            if len(values_x) > 20:
               p = numpy.polyfit(values_x, values_y, pfitdegree)
               fitted_values_y = numpy.polyval(p, values_x)
               
               values_x = numpy.extract(fitted_values_y >=0, values_x)
               fitted_values_y = numpy.extract(fitted_values_y >=0, fitted_values_y)
               str_p = str([string.strip('%7.3f' %e) for e in p])
               fit_info = "Poly deg: " + str(pfitdegree) + " Fit: " + str_p
         
               # Plot the data
               plt.plot(values_x, fitted_values_y, linestyle=pfitlinestyle, color=pfitcolor,\
                          marker=pfitmarker, markersize=pfitmarkersize, 
                          markerfacecolor=pfitmarkerfacecolor, linewidth=pfitlinewidth)
               plt.annotate(fit_info, (0.05, 0.93), xycoords='axes fraction', fontsize=10,  
                            bbox={'facecolor':'#FFFFFF', 'alpha':0.5, 'pad':10})

         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename), 
                                     facecolor=fig.get_facecolor(), edgecolor="k", 
                                     transparent=True)
            plt.clf()
            plt.close()

      except (OverflowError):        
            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
            #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))


   # -----------------------------------------------------------------------------------------------
   def plot_quiver(self, x_coords, y_coords, values_x, values_y, 
                   output_file='quiver_plot', output_dir='.', plot_title='Quiver plot',
                   pfilter=None, pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                   x_label=None, y_label=None,
                   color='#000080', scale=5, headwidth=0, width=0.001, pivot='middle',
                   is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None, show=False):
      """! Make a quiver plot """

      try:

         rc('text', usetex=True)
         rc('mathtext', fontset='custom')
         rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure

         if not plot_title is None:
            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = r'Quiver plot of %s' %(value_names)

            t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center', \
                                 fontproperties=FontProperties(size=14, weight='bold'))

         plt.grid(has_grid)
         if has_equalaxes:
            axes = fig.gca() 
            axes.set_aspect('equal')

         values_x = numpy.asarray(values_x)
         values_y = numpy.asarray(values_y)

         if pfilter != None:
            filtered_values_x = values_x[pfilter]
            filtered_values_y = values_y[pfilter]
            if len(filtered_values_x) > 0:
               values_x = filtered_values_x
            if len(filtered_values_y) > 0:
               values_y = filtered_values_y

         #plt.quiver(x_coords, y_coords, values_x, values_y, scale=0.1, headwidth=0, width=0.01, pivot='middle')
         q_plot = plt.quiver(x_coords, y_coords, values_x, values_y, scale=scale, headwidth=0, width=width, pivot='middle', color=color)

         plt.xlabel(x_label, fontproperties=FontProperties(size=14, weight='bold'))
         plt.ylabel(y_label, fontproperties=FontProperties(size=14, weight='bold'))

         if pxmin_lim is None:
            pxmin_lim = min(x_coords)
         if pxmax_lim is None:
            pxmax_lim = max(y_coords)
         if pymin_lim is None:
            pymin_lim = min(values_y)
         if pymax_lim is None:
            pymax_lim = max(values_y)

         if not has_equalaxes:
            pxmin = pxmin_lim
            pymin = pymin_lim
            pxmax = pxmax_lim
            pymax = pymax_lim
         else:
            pxmin = pymin = min(pxmin_lim, pymin_lim) 
            pxmax = pymax = max(pxmax_lim, pymax_lim) 

         plt.xlim(pxmin, pxmax)
         plt.ylim(pymin, pymax)

         if show:
            plt.show()

         if pfigure is None:
            plt.savefig(os.path.join(output_dir, output_file))
            plt.clf()
            plt.close()

      except (OverflowError):        
            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
            #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))

      return q_plot

   # -----------------------------------------------------------------------------------------------
   def plot_array(self, data_dico, value_names,
                  plot_filename, output_directory, plot_title='', pfilter=None, 
                  first_index=-1, last_index=-1,
                  x_label='', y_label='', x_tick_labels=None, y_tick_labels=None, 
                  has_grid=True, pfigure=None,
                  pcolor='blue', pmarker='', pmarkersize=1.0, pmarkerfacecolor='red', 
                  plinestyle='-', plinewidth=1, is_normed=False, has_equalaxes=False,
                  has_lstsqfit=False, pfitdegree=2,
                  pfitcolor='green', pfitmarker='', pfitmarkersize=1.0, pfitmarkerfacecolor='black',
                  pfitlinewidth=2, pfitlinestyle='-.'):
                          
      try:

         if pfigure is None:
            fig = plt.figure()
         else:
            fig = pfigure
         plt.grid(has_grid)
         axes = fig.gca() 
         if has_equalaxes:
            axes.set_aspect('equal')

         if x_tick_labels is not None:
            axes.set_xticklabels(x_tick_labels)           
         if y_tick_labels is not None:
            axes.set_yticklabels(y_tick_labels)           

         iplot = 1

         if not plot_title is None:
            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = 'Scatter plot of %s: ' %(value_names)

            if pfigure is None:
               plt.gcf().text(0.5, 0.93, figtitle, horizontalalignment='center', 
               fontproperties=FontProperties(size=11, weight='bold'))

         # Two values to represent
         if value_names[0] in data_dico and value_names[1] in data_dico:
            values_x = numpy.array(data_dico[value_names[0]])
            values_y = numpy.array(data_dico[value_names[1]])

            values_x = numpy.nan_to_num(values_x)
            values_y = numpy.nan_to_num(values_y)

            if pfilter != None:
               filtered_values_x = values_x[pfilter]
               filtered_values_y = values_y[pfilter]
               if len(filtered_values_x) > 0:
                  values_x = filtered_values_x
               if len(filtered_values_y) > 0:
                  values_y = filtered_values_y

            min_index = 0
            max_index = len(values_x)
            if first_index > -1 and first_index < max_index:
               min_index = first_index

            if last_index > -1 and last_index < max_index:
               max_index = last_index

            values_x = values_x[min_index:max_index]
            values_y = values_y[min_index:max_index]     

         #   if first > -1 and len(values_x) > max_len and len(values_y) > max_len :
         #      if not last_values:
         #         # Take first  <max_len> values
         #         values_x = values_x[0:max_len]
         #         values_y = values_y[0:max_len]     
         #      else:
         #         # Take last  <max_len> values
         #         values_x = values_x[-max_len:len(values_x)]
         #         values_y = values_y[-max_len:len(values_y)]     

            # Plot the data
            plt.plot(values_x, values_y,  linestyle=plinestyle, color=pcolor,\
                     marker=pmarker, markersize=pmarkersize, markerfacecolor=pmarkerfacecolor, 
                     linewidth=plinewidth)

            # plt.plot(values_x, values_y,  linestyle='-', color='green',  marker='', 
            # markersize=2.5, markerfacecolor='red')


            # Also include a least-squares fit if requested
            if has_lstsqfit: 
               if len(values_x) > 20:
                  p = numpy.polyfit(values_x, values_y, pfitdegree)
                  fitted_values_y = numpy.polyval(p, values_x)
                  
                  values_x = numpy.extract(fitted_values_y >=0, values_x)
                  fitted_values_y = numpy.extract(fitted_values_y >=0, fitted_values_y)
                  str_p = str([string.strip('%7.3f' %e) for e in p])
                  fit_info = '(degree ' + str(pfitdegree) + ' lstsq fit: ' + str_p + ')'
            
                  # Plot the data
                  plt.plot(values_x, fitted_values_y, linestyle=pfitlinestyle, color=pfitcolor,\
                             marker=pfitmarker, markersize=pfitmarkersize, 
                             markerfacecolor=pfitmarkerfacecolor, linewidth=pfitlinewidth)
               else:
                  fit_info = ''
                  print 'Warning: Not enough values to fit the data'

               figtitle += '\n' + fit_info     

            plt.xlabel(x_label, fontproperties=FontProperties(size=11, weight='bold'))
            plt.ylabel(y_label, fontproperties=FontProperties(size=11, weight='bold'))    

            if pfigure is None:
               plt.savefig(os.path.join(output_directory, plot_filename))
               plt.clf()
               plt.close()

         else:
            print 'Error: %s or %s not found in data dictionary' %(value_names[0], value_names[1])

      except Exception as detail:        
            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
            # logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))


   # -----------------------------------------------------------------------------------------------
   def create_image_plots(self, object_per_type_dico, file_type, job, master):
      """!
         Create plots from the data collected during shape measurement 
      """

      # TODO: model-dependent plots should be moved to the model classes (multifit)

      # --- Where plots will be stored
      branch_tree = job.get_branch_tree()
      plot_output_dir = os.path.join(master.plot_output_dir, branch_tree)

      if master.logging_enabled():
         master.logger.log_info_p(
           "{0} - {1}/img-{2:03}-{3:1d} - {4} - Generating plots...".format(
           master, branch_tree, job.img_no, job.epoch,  file_type))

      # ----------------------------------- Results from Fitting -----------------------------------

      result_dico = object_per_type_dico[file_type]["result"]   # shape measurement results
      filter_dico = object_per_type_dico[file_type]["filtered"] # filtered results
      info_dico = object_per_type_dico[file_type]["info"]       # information about the measurement 
      model_dico = object_per_type_dico[file_type]["model"]     # information about the model

      if len(result_dico["GAL_id"]) > 0:

         # --- Result filtering
         failed_value = master.config.get_as_float("FAILED_ELLIPTICITY_VALUE", "SHAPE_MEASUREMENT")   
         valid_ellipticity_filter = invalid_ellipticity_filter = None
         if failed_value != -1.0:   
            #result_dico = object_per_type_dico[file_type]["result"]
            valid_ellipticity_filter = numpy.logical_and(
                                       numpy.asarray(result_dico["e1"]) != failed_value,
                                       numpy.asarray(result_dico["e2"]) != failed_value,
                                       numpy.asarray(result_dico["flag"]) == 0)
            #invalid_ellipticity_filter = numpy.logical_not(valid_ellipticity_filter)

         # --- Histograms
         self.create_histogram_plots(result_dico.copy(), info_dico, model_dico, 
                                           valid_ellipticity_filter, 
                                           file_type, plot_output_dir, job, master, tag="")

         # --- Scatter plots
         self.create_scatter_plots(result_dico, info_dico, model_dico, 
                                   valid_ellipticity_filter,
                                   file_type, plot_output_dir, job, master, tag="")

      #print "filter_dico", filter_dico

      if len(filter_dico["GAL_id"]) > 0:

         # --- Filtered object "filtering"
         valid_ellipticity_filter = numpy.logical_and(
                                                  numpy.asarray(filter_dico["e1"]) != failed_value,
                                                  numpy.asarray(filter_dico["e2"]) != failed_value)
         #invalid_ellipticity_filter = numpy.logical_not(valid_ellipticity_filter)

         # --- Histograms
         self.create_histogram_plots(filter_dico.copy(), info_dico, model_dico, 
                                           valid_ellipticity_filter, 
                                           file_type, plot_output_dir, 
                                           job, master, tag="*filtered*")

         # --- Scatter plots
         self.create_scatter_plots(filter_dico, info_dico, model_dico, 
                                         valid_ellipticity_filter,
                                         file_type, plot_output_dir, 
                                         job, master, tag="*filtered*")
       
      # TODO: also plot data from failed fits (SE_*)
         


   # -----------------------------------------------------------------------------------------------
   def create_histogram_plots(self, data_dico, info_dico, model_dico, pfilter, 
                                    file_type, plot_output_dir, job, master, tag):
      """!
         Create histogram plots 
      """

      # --- Data sources
#      data_dico = object_per_type_dico[file_type]["result"] # shape measurement results
#      info_dico = object_per_type_dico[file_type]["info"]     # information about the measurement 
#      model_dico = object_per_type_dico[file_type]["model"]   # information about the model

      # --- Histograms of Fitted Parameters
      full_var_list = [] 
      except_list = []

      if info_dico["success_count"] > 0:
         param_names  = model_dico["param_names"]
         param_bounds = model_dico["param_bounds"]

         for (param_name, param_bound) in zip(param_names, param_bounds):
            if not param_name in except_list:
               if not pfilter is None:
                  data_dico[param_name] = numpy.asarray(data_dico[param_name])[pfilter]            
                    
               self.create_histogram_plot(data_dico, param_name, param_bound, file_type, 
                                                plot_output_dir, job, master, 
                                                pcolor="#4169e1", tag=tag)

         full_var_list.extend(param_names)

      else:
         param_names = []


      # --- Extra variables produced by fitting (except SExtractor)
      except_list = ["GAL_id", "iobj", "GAL_x", "GAL_y", "GAL_Xc", "GAL_Yc", "GAL_RA", "GAL_DEC",\
                     "gal_stamp_size", "flag", "weight", "true_rel_centroids"]

      var_filter = "not var.startswith('SE_') and not var.endswith('galaxy') "\
                   "and not var in param_names and not var in except_list" 
      var_names  = [var for var in data_dico.keys() if eval(var_filter)]

      # --- Plot the histogram of each variable, applying the input filter
      for var_name in var_names:
         if not pfilter is None:
            data_dico[var_name] = numpy.asarray(data_dico[var_name])[pfilter]            
         self.create_histogram_plot(data_dico, var_name, None, file_type, 
                                          plot_output_dir, job, master, max_xlim=150, 
                                          pcolor="#6495ed", tag=tag)

      full_var_list.extend(var_names)

      # --- SExtractor variables
      se_except_list = ["SE_X_IMAGE", "SE_Y_IMAGE", "SE_ALPHA_J2000", "SE_DELTA_J2000",\
                        "SE_GAL_Xc", "SE_GAL_Yc", "SE_GAL_xc", "SE_GAL_yc",\
                        "SE_PSF_xc", "SE_PSF_yc",
                        "SE_GAL_FLUXERR", "SE_PSF_FLUXERR"]
      se_var_filter = "var.startswith('SE_') and not var in param_names and not var in var_names "\
                      "and not var in se_except_list" 
      se_var_names = [var for var in data_dico.keys() if eval(se_var_filter)]

      for var_name in se_var_names:
         if not pfilter is None:
            data_dico[var_name] = numpy.asarray(data_dico[var_name])[pfilter]            

         self.create_histogram_plot(data_dico, var_name, None, 
                                          file_type, plot_output_dir, job, master,
                                          pcolor="#ffd700", tag=tag)
      full_var_list.extend(se_var_names)

      # --- Derived variables

      #deriv_var_filter = "not var in full_var_list"
      #deriv_var_names = [var for var in data_dico.keys() if eval(deriv_var_filter)]

      #self.create_histogram_plot(data_dico, var_name, None, 
      #                                 file_type, plot_output_dir, job, master,
      #                                 pcolor="#2e8b57", tag=tag)

#      # --- Useful variables
#      except_list = []

#      min_max_residuals_var_filter = "var not in full_var_list and var.endswith('residuals')"
##      min_max_residuals_var_filter = min_max_residuals_var_filter and
##                                     "(var.startswith('min') or var.startswith('max') or var.startswith('sum')))"
##      min_max_residuals_var_filter = "var not in full_var_list and "\
##                                     "(var.endswith('residuals') and "\
##                                     "(var.startswith('min') or var.startswith('max') or var.startswith('sum')))" 
##      min_max_residuals_var_filter = "var not in full_var_list and "\
##                                     "(var.endswith('residuals') or var.endswith('galaxy')) and "\
##                                     "(var.startswith('min') or var.startswith('max'))" 
#      min_max_residuals_var_names =\
#          [var for var in data_dico.keys() if eval(min_max_residuals_var_filter)]


#      for var_name in min_max_residuals_var_names:
#         self.create_histogram_plot(data_dico, var_name, None, tag, file_type, 
#                                          plot_output_dir, job, master)      
#      full_var_list.extend(min_max_residuals_var_names)


   # -----------------------------------------------------------------------------------------------
   def create_histogram_plot(self, data_dico, var_name, var_bound, 
                                         file_type, plot_output_dir, job, master, pcolor="b",
                                         max_xlim=150, tag=""):
      """! Create an histogram plot from the data dictionary of an image """
               
      if var_name in data_dico and len(data_dico[var_name]) > 0:

         # --- Mean_value of main variable
         values = numpy.nan_to_num(data_dico[var_name])
         var_mean = numpy.mean(values)                
         var_std  = numpy.std(values)                
         var_min  = numpy.min(values)                
         var_max  = numpy.max(values)                

         # --- Plot title and filename
         branch_tree = job.get_branch_tree()

         file_main, file_ext = os.path.splitext(file_type)
         plot_title = "{0} {1} - {2:03d}-{3:1d} - "\
                      "Mean/Var: [{4:.3f}, {5:.3f}]  Min/Max: [{6:.3f}, {7:.3f}]".format(
                            var_name, tag, job.img_no, job.epoch, 
                            var_mean, var_std, var_min, var_max)

         plot_filename = "hist_{0}_{1}_img_{2:03d}-{3:1d}_"\
                         "{4}_{5:.3f}_{6:.3f}_{7:.3f}_{8:.3f}_{9}.png".format(
             var_name, "_".join(job.branch_dirs), job.img_no, job.epoch, 
             file_type, var_mean, var_std, var_min, var_max, tag.replace("*", ""))

         plot_title.replace(" -", "-")
         plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

         (pmin_xlim, pmax_xlim) = (None, None)

         self.plot_histogram(data_dico, [var_name],
                             plot_filename, plot_output_dir, plot_title, pfilter=None,
                             pnb_bins=20, pmin_xlim=pmin_xlim, pmax_xlim=pmax_xlim,
                             pcolor=pcolor, is_normed=False, has_grid=False)

  
   # -----------------------------------------------------------------------------------------------
   def create_scatter_plots(self, data_dico, info_dico, model_dico, 
                                        pfilter,
                                        file_type, plot_output_dir, job, master, tag):
      """!
         Create scatter plots 
      """

      branch_tree = job.get_branch_tree()

      # --- Ellipticity (e1, e2) if they exist in the result dico
      #
      # TODO: Refactor, move to the model-dependent code, e.g. gbdsersic.py
      #

      file_main, file_ext = os.path.splitext(file_type)

      if "e1" in data_dico and "e2" in data_dico:      

         e1_label = r"$\mathregular{e_1}$"
         e2_label = r"$\mathregular{e_2}$"
         e1_e2_label = r"$\mathregular{(e_1,\,e_2)}$"

         plot_title = "Ellipticity {0} {1} - {2:03d}-{3:1d}{4}".format(
                          e1_e2_label, tag, job.img_no, job.epoch, file_ext)

         plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{5:1d}{6}_{7}.png".format(
                                           "e1", "e2",
                                           "_".join(job.branch_dirs), file_main,
                                           job.img_no, job.epoch, file_ext, 
                                           tag.replace("*",""))
         plot_title.replace(" -", "-")
         plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

         self.plot_scatter(data_dico, ["e1", "e2"],
                    plot_filename, plot_output_dir, plot_title=plot_title, 
                    pfilter=pfilter,
                    pxmin_lim=-1.0, pxmax_lim=1.0, pymin_lim=-1.0, pymax_lim=1.0, 
                    x_label=e1_label, y_label=e2_label, pnorm=None, pcmap='jet',
                    pcolor='b', pmarker='o', pmarkersize=5, 
                    pmarkerfacecolor=None, pmarkeredgecolor=None, 
                    px_sci=False, py_sci=False, pminor=True, pstats=False,
                    is_normed=False, has_equalaxes=True, has_grid=True, pfigure=None)

      # --- SE_FWHM_IMAGE vs SE_FLUX_RADIUS
      if "SE_GAL_FWHM_IMAGE" in data_dico and "SE_GAL_FLUX_RADIUS" in data_dico:

         plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                         "SE  FWHM", "SE half-light radius",
                          tag, job.img_no, job.epoch, file_ext)
         plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                           "SE_FWHM", "SE_FLUX_RADIUS",
                                           "_".join(job.branch_dirs), file_main,
                                           job.img_no, job.epoch, file_ext, 
                                           tag.replace("*",""))
         plot_title.replace(" -", "-")
         plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

         self.plot_scatter(data_dico, ["SE_GAL_FLUX_RADIUS", "SE_GAL_FWHM_IMAGE"],
                           plot_filename, plot_output_dir, plot_title=plot_title, 
                           pfilter=pfilter,
                           pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                           x_label="SE half-light radius", y_label="SE FWHM",  
                           pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                           pmarkerfacecolor=None, pmarkeredgecolor=None, 
                           px_sci=False, py_sci=False, pminor=True, pstats=False,
                           is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


      # --- SE_FLUX vs fitted flux
      if "SE_GAL_FLUX" in data_dico and "flux" in data_dico:

         plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                         "SE Flux", "Total fitted flux",
                          tag, job.img_no, job.epoch, file_ext)
         plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                           "Total flux ", "SE_GAL_FLUX",
                                           "_".join(job.branch_dirs), file_main,
                                           job.img_no, job.epoch, file_ext, 
                                           tag.replace("*",""))
         plot_title.replace(" -", "-")
         plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

         self.plot_scatter(data_dico, ["SE_GAL_FLUX", "flux"],
                           plot_filename, plot_output_dir, plot_title=plot_title, 
                           pfilter=pfilter,
                           pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                           x_label="SE Flux", y_label="Total fitted flux",  
                           pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                           pmarkerfacecolor=None, pmarkeredgecolor=None, 
                           px_sci=False, py_sci=False, pminor=True, pstats=False,
                           has_lstsqfit=True, pfitdegree=1,
                           is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


      # --- SE GAL_FLUX_RADIUS vs fitted disk half light radius
      if "SE_GAL_FLUX_RADIUS" in data_dico and "dre" in data_dico:

         plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                         "SE half-light radius", "Fitted disk effective radius",
                          tag, job.img_no, job.epoch, file_ext)
         plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                           "dre ", "SE_GAL_FLUX_RADIUS",
                                           "_".join(job.branch_dirs), file_main,
                                           job.img_no, job.epoch, file_ext, 
                                           tag.replace("*",""))
         plot_title.replace(" -", "-")
         plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

         self.plot_scatter(data_dico, ["SE_GAL_FLUX_RADIUS", "dre"],
                           plot_filename, plot_output_dir, plot_title=plot_title, 
                           pfilter=pfilter,
                           pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                           x_label="SE half-light radius", y_label="Fitted disk effective radius",  
                           pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                           pmarkerfacecolor=None, pmarkeredgecolor=None, 
                           px_sci=False, py_sci=False, pminor=True, pstats=False,
                           has_lstsqfit=True, pfitdegree=1,
                           is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


      # --- Residuals versus various features...

      if "reduced_chi2" in data_dico:      

         reduced_chi2 = numpy.nan_to_num(data_dico["reduced_chi2"])

         pfilter = numpy.logical_and(pfilter, numpy.asarray(reduced_chi2) > 0.0)

         # --- Reduced Chi^2 vs SE Flux
         if "SE_GAL_FLUX" in data_dico:      


            plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                            "Reduced chi^2", "SE Flux",
                             tag, job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "reduced_chi2", "SE_GAL_FLUX",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))

            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["SE_GAL_FLUX", "reduced_chi2"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                              x_label="SE Flux", y_label="Reduced Chi^2",  
                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

         # --- Reduced Chi^2 vs SE SNR
         if "SE_GAL_SNR" in data_dico:

            plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                            "Reduced_chi^2", "SE SNR",
                             tag, job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "reduced_chi2", "SE_GAL_SNR",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))
      
            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["SE_GAL_SNR", "reduced_chi2"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                              x_label="SE SNR", y_label="Reduced Chi^2",  
                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

            max_y_zoom_limit = 200
            if numpy.max(data_dico["reduced_chi2"]) > max_y_zoom_limit: 
               plot_filename = "scatter_{0}_{1}_{2}_img_{3:03d}-{4:1d}_{5}_{6}_zoom.png".format(
                                                 "reduced_chi2", "SNR",
                                                 job.branch_dirs, job.img_no, job.epoch, file_type, 
                                                 tag.replace("*",""))
               plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

               self.plot_scatter(data_dico, ["SE_GAL_SNR", "reduced_chi2"],
                                 plot_filename, plot_output_dir, plot_title=plot_title, 
                                 pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                                 x_label="SE SNR", y_label="Reduced Chi^2",  
                                 pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                                 pmarkerfacecolor=None, pmarkeredgecolor=None, 
                                 px_sci=False, py_sci=False, pminor=True, pstats=False,
                                 is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)



#         # --- Residuals vs SE Flux Radius
#         if "SE_FLUX_RADIUS" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Galaxy Flux Radius", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals", "flux_radius",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_FLUX_RADIUS", "sum_norm_residuals"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=0, pymax_lim=None, 
#                              x_label="SE_FLUX_RADIUS", y_label="Normalized Residuals / Pixel", 
#                               pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

##         # --- Residuals vs SE FLUX_RADIUS_ratio
##         if "SE_FLUX_RADIUS_ratio" in data_dico:      
##            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
##                                              "Residuals", "Rg/Rp Flux Radius", tag,
##                                              job.img_no, job.epoch, branch_tree, file_type)
##            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
##                                              "normed_residuals", "flux_radius_ratio",
##                                              job.branch, job.obs_type, job.data_type, 
##                                              job.img_no, job.epoch, file_type,
##                                              tag.replace("*",""))
##            plot_title.replace(" -", "-")
##            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

##            self.plot_scatter(data_dico, ["SE_FLUX_RADIUS_ratio", "sum_norm_residuals"],
##                              plot_filename, plot_output_dir, plot_title=plot_title, 
##                              pfilter=pfilter,
##                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
##                              x_label="Rg/Rp Flux Radius", y_label="Normalized Residuals / Pixel", 
##                               pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
##                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
##                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#         # --- Residuals vs SE FWHM
#         if "SE_FWHM_IMAGE" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "galaxy FWHM", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals", "galaxy_FWHM",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_FWHM_IMAGE", "sum_norm_residuals"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=0, pymax_lim=None, 
#                              x_label="SE_FWHM_IMAGE", y_label="Normalized Residuals / Pixel",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

##         # --- Residuals vs SE FWHM_ratio
##         if "SE_FWHM_ratio" in data_dico:      
##            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
##                                              "Residuals", "Rg/Rp FWHM", tag,
##                                              job.img_no, job.epoch, branch_tree, file_type)
##            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
##                                              "normed_residuals", "FWHM_ratio",
##                                              job.branch, job.obs_type, job.data_type, 
##                                              job.img_no, job.epoch, file_type,
##                                              tag.replace("*",""))
##            plot_title.replace(" -", "-")
##            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

##            self.plot_scatter(data_dico, ["SE_FWHM_ratio", "sum_norm_residuals"],
##                              plot_filename, plot_output_dir, plot_title=plot_title, 
##                              pfilter=pfilter,
##                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
##                              x_label="Rg/Rp FWHM", y_label="Normalized Residuals / Pixel",  
##                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
##                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
##                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#         # --- Residuals vs flux
#         if "flux" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Galaxy flux", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals", "flux",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["flux", "sum_norm_residuals"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=0, pymax_lim=None, 
#                              x_label="flux", y_label="Normalized Residuals / Pixel",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#            max_x_zoom_limit = 200
#            if numpy.max(data_dico["flux"]) > max_x_zoom_limit: 
#               plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}_zoom.png".format(
#                                                 "normed_residuals", "flux",
#                                                 job.branch, job.obs_type, job.data_type, 
#                                                 job.img_no, job.epoch, file_type,
#                                                 tag.replace("*",""))
#               plot_title.replace(" -", "-")
#               plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#               self.plot_scatter(data_dico, ["flux", "sum_norm_residuals"],
#                                 plot_filename, plot_output_dir, plot_title=plot_title, 
#                                 pfilter=pfilter,
#                                 pxmin_lim=0, pxmax_lim=max_x_zoom_limit, 
#                                 pymin_lim=0, pymax_lim=None, 
#                                 x_label="flux", y_label="Normalized Residuals / Pixel",  
#                                 pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                                 pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                                 is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


         # --- Residuals vs dre
         if "dre" in data_dico:      

            plot_title = "{0} vs. {1} {2} - {3:03d}-{4:1d}{5}".format(
                            "Residuals", "Disk effective radius",
                             tag,  job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "normed_residuals", "dre",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))

            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["dre", "sum_norm_residuals"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=0, pymax_lim=None, 
                              x_label="disk effective radius (dre parameter)", 
                              y_label="Normalized Residuals / Pixel",  
                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

         # --- Residuals vs bre
         if "bre" in data_dico:      
            plot_title = "{0} vs. {1} {2} - {3}.{4:03d}-{5:1d}{6}".format(
                            "Residuals", "Bulge effective radius",
                             tag, file_main, job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "normed_residuals", "bre",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))

            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["bre", "sum_norm_residuals"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=0, pymax_lim=None, 
                              x_label="Bulge effective radius (bre parameter)", 
                              y_label="Normalized Residuals / Pixel",  
                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


#      if "rmsd_norm_residuals_2" in data_dico:      

#         pfilter = numpy.logical_and(pfilter,  
#                                     numpy.asarray(data_dico["rmsd_norm_residuals_2"]) > 0.0)

#         # --- Residuals vs SE Flux
#         if "SE_FLUX" in data_dico:      

#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals 2", "SE_FLUX", tag, 
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "SE_FLUX",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type, 
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_FLUX", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="SE_FLUX", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


#         # --- Residuals vs SE SNR
#         if "SE_SNR" in data_dico:      

#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals 2", "SNR", tag, 
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "SNR",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_SNR", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="SE_SNR", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#            max_x_zoom_limit = 200
#            if numpy.max(data_dico["rmsd_norm_residuals_2"]) > max_x_zoom_limit: 
#               plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}_zoom.png".format(
#                                                 "residuals_2", "SNR",
#                                                 job.branch, job.obs_type, job.data_type, 
#                                                 job.img_no, job.epoch, file_type, 
#                                                 tag.replace("*",""))
#               plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#               self.plot_scatter(data_dico, ["SE_SNR", "rmsd_norm_residuals_2"],
#                                 plot_filename, plot_output_dir, plot_title=plot_title, 
#                                 pfilter=pfilter,
#                                 pxmin_lim=0, pxmax_lim=max_x_zoom_limit, 
#                                 pymin_lim=None, pymax_lim=None, 
#                                 x_label="SE_SNR", y_label="Normalized residuals",  
#                                 pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
#                                 pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                                 is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


#         # --- Residuals vs SE Flux Radius
#         if "SE_FLUX_RADIUS" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Galaxy Flux Radius", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "flux_radius",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_FLUX_RADIUS", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="SE_FLUX_RADIUS", y_label="Normalized Residuals / Total Flux", 
#                               pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

##         # --- Residuals vs SE FLUX_RADIUS_ratio
##         if "SE_FLUX_RADIUS_ratio" in data_dico:      
##            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
##                                              "Residuals", "Rg/Rp Flux Radius", tag,
##                                              job.img_no, job.epoch, branch_tree, file_type)
##            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
##                                              "normed_residuals_2", "flux_radius_ratio",
##                                              job.branch, job.obs_type, job.data_type, 
##                                              job.img_no, job.epoch, file_type,
##                                              tag.replace("*",""))
##            plot_title.replace(" -", "-")
##            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

##            self.plot_scatter(data_dico, ["SE_FLUX_RADIUS_ratio", "rmsd_norm_residuals_2"],
##                              plot_filename, plot_output_dir, plot_title=plot_title, 
##                              pfilter=pfilter,
##                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
##                              x_label="Rg/Rp Flux Radius", y_label="Normalized Residuals / Total flux", 
##                               pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
##                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
##                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#         # --- Residuals vs SE FWHM
#         if "SE_FWHM_IMAGE" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "galaxy FWHM", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "galaxy_FWHM",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["SE_FWHM_IMAGE", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="SE_FWHM_IMAGE", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

##         # --- Residuals vs SE FWHM_ratio
##         if "SE_FWHM_ratio" in data_dico:      
##            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
##                                              "Residuals", "Rg/Rp FWHM", tag,
##                                              job.img_no, job.epoch, branch_tree, file_type)
##            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
##                                              "normed_residuals_2", "FWHM_ratio",
##                                              job.branch, job.obs_type, job.data_type, 
##                                              job.img_no, job.epoch, file_type,
##                                              tag.replace("*",""))
##            plot_title.replace(" -", "-")
##            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

##            self.plot_scatter(data_dico, ["SE_FWHM_ratio", "rmsd_norm_residuals_2"],
##                              plot_filename, plot_output_dir, plot_title=plot_title, 
##                              pfilter=pfilter,
##                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
##                              x_label="Rg/Rp FWHM", y_label="Normalized Residuals / Total Flux",  
##                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
##                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
##                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#         # --- Residuals vs flux
#         if "flux" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Galaxy flux", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "flux",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["flux", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="flux", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


#         # --- Residuals vs dre
#         if "dre" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Disk radius", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "dre",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["dre", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="dre", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

#         # --- Residuals vs bre
#         if "bre" in data_dico:      
#            plot_title = "{0} vs {1} {2} - Image {3:03d}-{4:1d} - /{5} - {6}".format(
#                                              "Residuals", "Bulge radius", tag,
#                                              job.img_no, job.epoch, branch_tree, file_type)
#            plot_filename = "scatter_{0}_{1}_{2}_{3}_{4}_img_{5:03d}-{6:1d}_{7}_{8}.png".format(
#                                              "normed_residuals_2", "bre",
#                                              job.branch, job.obs_type, job.data_type, 
#                                              job.img_no, job.epoch, file_type,
#                                              tag.replace("*",""))
#            plot_title.replace(" -", "-")
#            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

#            self.plot_scatter(data_dico, ["bre", "rmsd_norm_residuals_2"],
#                              plot_filename, plot_output_dir, plot_title=plot_title, 
#                              pfilter=pfilter,
#                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                              x_label="bre", y_label="Normalized Residuals / Total Flux",  
#                              pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=4, 
#                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
#                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)


      if "max_norm_residuals" in data_dico:      

         pfilter = numpy.logical_and(pfilter, numpy.asarray(data_dico["max_norm_residuals"]) > 0.0)

         # --- Residuals vs SE Flux
         if "SE_GAL_FLUX" in data_dico:      


            plot_title = "{0} vs. {1} {2} - {3}.{4:03d}-{5:1d}{6}".format(
                            "Max Residuals", "SE Flux",
                             tag, file_main, job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "max_normed_residuals", "SE_GAL_FLUX",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))

            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["SE_GAL_FLUX", "max_norm_residuals"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                              x_label="SE Flux", y_label="Maximum normalized Residuals / Pixel",  
                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

         # --- Residuals vs SE SNR
         if "SE_GAL_SNR" in data_dico:

            plot_title = "{0} vs. {1} {2} - {3}.{4:03d}-{5:1d}{6}".format(
                            "Residuals", "SE SNR",
                             tag, file_main, job.img_no, job.epoch, file_ext)
            plot_filename = "scatter_{0}_{1}_{2}_{3}{4:03d}-{4:1d}{5}_{6}.png".format(
                                              "max_normed_residuals", "SE_GAL_SNR",
                                              "_".join(job.branch_dirs), file_main,
                                              job.img_no, job.epoch, file_ext, 
                                              tag.replace("*",""))
      
            plot_title.replace(" -", "-")
            plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

            self.plot_scatter(data_dico, ["SE_GAL_SNR", "max_norm_residuals"],
                              plot_filename, plot_output_dir, plot_title=plot_title, 
                              pfilter=pfilter,
                              pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                              x_label="SE SNR", y_label="Maximum normalized Residuals / Pixel",  
                              pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                              pmarkerfacecolor=None, pmarkeredgecolor=None, 
                              px_sci=False, py_sci=False, pminor=True, pstats=False,
                              is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)

            max_y_zoom_limit = 200
            if numpy.max(data_dico["max_norm_residuals"]) > max_y_zoom_limit: 
               plot_filename = "scatter_{0}_{1}_{2}_img_{3:03d}-{4:1d}_{5}_{6}_zoom.png".format(
                                                 "max_normed residuals", "SNR",
                                                 job.branch_dirs, job.img_no, job.epoch, file_type, 
                                                 tag.replace("*",""))
               plot_filename = plot_filename.replace("__", "_").replace("_.png", ".png")

               self.plot_scatter(data_dico, ["SE_GAL_SNR", "max_norm_residuals"],
                                 plot_filename, plot_output_dir, plot_title=plot_title, 
                                 pfilter=pfilter,
                                 pxmin_lim=None, pxmax_lim=None, 
                                 pymin_lim=None, pymax_lim=max_y_zoom_limit, 
                                 x_label="SE SNR", y_label="Maximum normalized Residuals / Pixel",  
                                 pnorm=None, pcmap='jet', pcolor='b', pmarker='o', pmarkersize=4, 
                                 pmarkerfacecolor=None, pmarkeredgecolor=None, 
                                 px_sci=False, py_sci=False, pminor=True, pstats=False,
                                 is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None)




#   # -----------------------------------------------------------------------------------------------
#   def create_image_correlation_matrix_plots(self, result_dico, info_dico, model_dico, 
#                                                   file_type, plot_output_dir, 
#                                                   job, master, tag="",
#                                                   color_bar=False, pcmap="jet", pfigure=None):

#      # --- All variables
#      plot_title = "{0} - Image {1:03d}-{2:1d} - /{3}".format(
#                   "Correlation Matrix", job.img_no, job.epoch, job.get_branch_tree())
#      plot_filename = "corr_matrix_{0}_{1}_{2}_img_{3:03d}-{4:1d}.png".format(
#                               job.branch, job.obs_type, job.data_type, job.img_no, job.epoch)

#      result_dico["e"] = numpy.hypot(result_dico["e1"], result_dico["e2"])

#      var_list = ["rmsd_norm_residuals_2", "rmsd_norm_residuals", "dre", "bre", "flux", "e", "SE_FWHM_IMAGE",
#                  "SE_FLUX_RADIUS", "SE_SNR"] 
#      var_labels = ["Rnrm", "Rpix", "dre", "bre", "Itot", "|e|", "FWHM", "HLR", "SNR"]

#      self._create_image_correlation_matrix_plots(result_dico, info_dico, model_dico, file_type,
#                                                  plot_title, plot_filename, plot_output_dir, 
#                                                  job, master, var_list, var_labels, tag,
#                                                  color_bar, pcmap, pfigure)

#      # --- Fitted parameters
#      plot_title = "{0} - Image {1:03d}-{2:1d} - /{3}".format(
#                   "Parameter Correlation Matrix", job.img_no, job.epoch, job.get_branch_tree())
#      plot_filename = "corr_matrix_params_{0}_{1}_{2}_img_{3:03d}-{4:1d}.png".format(
#                               job.branch, job.obs_type, job.data_type, job.img_no, job.epoch)

#      result_dico["e"] = numpy.hypot(result_dico["e1"], result_dico["e2"])

#      var_list = ["dre", "bre", "flux", "e1", "e2", "e"] 
#      var_labels = ["    dre", "    bre", "    Itot", "     e1", "     e2", "     |e|"]

#      self._create_image_correlation_matrix_plots(result_dico, info_dico, model_dico, file_type,
#                                                  plot_title, plot_filename, plot_output_dir, 
#                                                  job, master, var_list, var_labels, tag,
#                                                  color_bar, pcmap, pfigure)



   # -----------------------------------------------------------------------------------------------
   def _import_module(self, module_name):
      """! Find and load a module <module_name> from a directory <module_dir> """

      try:

         file_obj, filename, data = imp.find_module(module_name)
         return imp.load_module(module_name, file_obj, filename, data)
 
      except:
         raise(ModelFitter.FittingError("Could not load module {0}.".format(module_name)))




   # ----------------------------------------------------------------------------------------------- 
   def _create_image_correlation_matrix_plots(self, result_dico, info_dico, model_dico, file_type, 
                                                    plot_title, plot_filename, plot_output_dir, 
                                                    job, worker, var_list, var_labels, tag="", 
                                                    color_bar=False, pcmap="jet", pfigure=None):

      try:  

         #rc('text', usetex=True)
         #rc('mathtext', fontset='custom')
         #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure

         branch_tree = job.get_branch_tree()
         plot_output_dir = os.path.join(worker.plot_output_dir, branch_tree)

         t = plt.gcf().text(0.5, 0.95, plot_title, horizontalalignment='center', \
                            fontproperties=FontProperties(size=11, weight='bold'))


         # --- Preprocess the data
         result_dico = deepcopy(result_dico)
         result_dico = self.helper.remove_nans_from_result_dico(result_dico, worker)

         data_list = [result_dico[k] for k in var_list if k in result_dico]
         data_stack = numpy.vstack(tuple(data_list))  

         corr_matrix = numpy.corrcoef(data_stack)

         plot = plt.pcolor(corr_matrix)
         if color_bar:
            norm1 = matplotlib.colors.Normalize(vmin=0,vmax=9)
            plt.colorbar(plot, cmap=pcmap, norm=norm1)
            plt.clim(-1.0, 1.0) 

#            cb_axes = fig.add_axes([-1, -1, 1, 1])
#            cb = matplotlib.colorbar.ColorbarBase(cb_axes, cmap=pcmap, norm=norm1)
#            plot.colorbar = cb 

         var_labels = ["{0:^7}".format(l) for l in var_labels]

         x_tick_labels = var_labels
         y_tick_labels = var_labels
          

         axes.set_xticklabels(x_tick_labels, horizontalalignment='left', fontsize="11")  
         axes.set_yticklabels(x_tick_labels, verticalalignment='bottom', fontsize="11", 
                                                                         rotation="vertical")



         plt.savefig(os.path.join(plot_output_dir, plot_filename))

      except Exception as detail:
         print("Error: exception thrown during plotting ({0})".format(detail))  


# -- EOF gfit_plot.py
