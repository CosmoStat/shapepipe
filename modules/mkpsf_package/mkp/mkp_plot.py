"""! 
   mkp_plot.py - Plotting functions.
"""

# -- Python imports
import os, sys
import math
import numpy
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as pylab3
from matplotlib.font_manager import FontProperties


# -------------------------------------------------------------------------------------------------
class MkpPlotter(object):
   
   """! 
      Class with a number of convenient plotting methods.   
   """

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
                           fontproperties=FontProperties(size=12, weight='bold'))
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
         
         plt.savefig(os.path.join(output_dir, output_file))

         if show:
            plt.show()

         plt.clf()
         plt.close()

      except:    
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
         @param xlabel labbel string for the x axis
         @param ylabel labbel string for the y axis
         @param zlabel labbel string for the z axis
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
                                     fontproperties=FontProperties(size=12, weight='bold'))
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
         print("Exception thrown during plotting ({0})".format(sys.exc_info()[1]))     

      finally:
         plt.clf()
         plt.close()


   # -----------------------------------------------------------------------------------------------
   def plot_histogram(self, data_dico, value_names,
                            plot_filename, output_directory, plot_title, pfilter=None,
                            pnb_bins=15.0, pmin_xlim=None, pmax_xlim=None, 
                            pcolor='b', is_normed=False, 
                            has_equalaxes=False, has_grid=True, pfigure=None):

      """! Plot histograms of possibly multiple arrays of floats found in a dictionary """

      try:

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure

         iplot = 1

         figtitle = plot_title

         if pfigure == None:
            t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center',
                                 fontproperties=FontProperties(size=12))

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
                  axes = fig.gca() 
                  axes.set_aspect('equal')

               values = numpy.array(data_dico[value_names[iplot-1]], dtype=numpy.float64)
               values = numpy.nan_to_num(values)

               if pfilter != None:
                  filtered_values = values[pfilter]
                  if len(filtered_values) > 0:
                     values = filtered_values
               color = pcolor
               x_label = value_names[iplot-1]
               y_label = ''

               plt.xlabel(x_label)
               plt.ylabel(y_label)

               if pmin_xlim is None:
                  pmin_xlim = min(values)

               if pmax_xlim is None:
                  plt.xlim(pmin_xlim, max(values))
               else:
                  plt.xlim(pmin_xlim, pmax_xlim)

               plt.hist(values, pnb_bins, normed=is_normed, fc=pcolor)
               
               iplot += 1
            
         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename))
            plt.clf()
            plt.close()

      except (OverflowError, ValueError, IndexError):       
         print("Exception thrown during plotting ({0})".format(sys.exc_info()[1]))     


   # -----------------------------------------------------------------------------------------------
   def make_plots(self, job_result, master):
         """! Plots histograms for PSE """

         job = job_result.job   
         object_per_type_dico = job_result.result

         file_types = object_per_type_dico.keys()   # all results for each image types

         for file_type in file_types:

            object_dico = object_per_type_dico[file_type]

            # --- Where plots will be stored
            branch_tree = job.get_branch_tree()
            plot_pathname = os.path.join(master.plot_output_dir, branch_tree)    

            if master.logging_enabled():
               master.logger.log_info_p(
                 "{0} - Img {1:03d} - /{2} - Generating plots for {3}...".format(
                                                  master, job.img_no, branch_tree, file_type))

            max_lim = 150  # max limit threshold for second plot

            # --- Plot Histograms for each column in the catalog
            for col_name in job_result.data_dico[file_type].keys():
               plot_title = "Histogram: {0} - Image {1:03d}-{2:1d} - /{3} - {4}".format(col_name, 
                                            job.img_no, job.epoch, branch_tree, file_type)
               plot_filename = "hist_{0}_{1}_img_{2:03d}-{3:1d}_{4}.png".format(col_name,
                               branch_tree, job.img_no, job.epoch, file_type)
   
               self.plot_histogram(job_result.data_dico[file_type], [col_name],
                                   plot_filename, plot_pathname, plot_title, 
                                   pnb_bins=25.0, pmax_xlim=None, pcolor='b', is_normed=False)

               if math.fabs(numpy.max(job_result.data_dico[file_type][col_name])) >= max_lim:
                  
                  plot_filename = "hist_{0}_{1}_img_{2:03d}-{3:1d}_{4}.png".format(col_name,
                                  branch_tree, job.img_no, job.epoch, file_type)

                  self.plot_histogram(job_result.data_dico[file_type], [col_name],
                                      plot_filename, plot_pathname, plot_title, 
                                      pnb_bins=25.0, pmin_xlim=0, pmax_xlim=max_lim, 
                                      pcolor='g', is_normed=False)


# -- EOF mkp_plot.py
