"""! 
   sp_plot.py - Plotting functions.
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
from matplotlib import rc


# -- External imports

# --- Module-specific imports


# -------------------------------------------------------------------------------------------------
class SpredictPlotter(object):
   
   """! 
      Class with a number of convenient plotting methods.   
   """

   # ~~~~~~~~~~~~~~
   # Public methods 
   # ~~~~~~~~~~~~~~


   # -----------------------------------------------------------------------------------------------
   def create_plots(self, model, method, result_dico, file_type, job, master):
      """!
         Create input plots 
      """

      file_main, file_ext = os.path.splitext(file_type)

      # --- Input plots
      if master.config.get_as_boolean("CREATE_PLOTS", "PLOTS.INPUT"):  
         if master.logging_enabled():
            master.logger.log_info_p(
             "{0} - /{1}/{2}-{3:03d}.{4} - Generating Input plots for parameter(s) "\
             "{5}\t({6}) ...".format(
                 master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
                 job.param_names, job.method_name))
         pass


      # --- Target plots
      if master.config.get_as_boolean("CREATE_PLOTS", "PLOTS.TARGET"):  
         if master.logging_enabled():
            master.logger.log_info_p(
              "{0} - /{1}/{2}-{3:03d}.{4} - Generating Target plots for parameter(s) "\
              "{5}\t({6}) ...".format(
                 master.name, job.get_branch_tree(), file_main, job.img_no, file_ext,
                 job.param_names, job.method_name))
         pass


      # --- Create model-specific plots
      model.create_plots(method, result_dico, file_type, job, master)


   # -----------------------------------------------------------------------------------------------
   def create_image_quiver_plots(self, x_coords, y_coords, x_values, y_values,
                                       x_var_name, y_var_name,
                                       file_type, plot_title, plot_output_dir, plot_filename,
                                       job, master, scale=5, width=0.001, pfilter=None):
      """!
         Create quiver (whisker) plots 
      """

      self.plot_quiver(x_coords, y_coords, x_values, y_values, 
                       output_file=plot_filename, output_dir=plot_output_dir, 
                       plot_title=plot_title, pfilter=pfilter,
                       pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                       x_label=x_var_name, y_label=y_var_name,
                       color='#000080', scale=scale, headwidth=0, width=width, pivot='middle',
                       is_normed=False, has_equalaxes=True, has_grid=True,
                       pfigure=None, show=False)



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
                            pnb_bins=15.0, pmin_xlim=None, pmax_xlim=None, 
                            pcolor='b', is_normed=False, 
                            has_equalaxes=False, has_grid=True, pfigure=None, logger=None):

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
                  plt.xlim(pmin_xlim, min(pmax_xlim, max(values)))

               plt.hist(values, pnb_bins, normed=is_normed, fc=pcolor)
               
               iplot += 1
            
         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename))
            plt.clf()
            plt.close()

      except (OverflowError, ValueError, IndexError):       
         if not logger is None:    
            logger.log_error_p('Exception thrown during plotting ({0})'.format(
                                                                           sys.exc_info()[1]))
         else:
            print("Exception thrown during plotting ({0})".format(sys.exc_info()[1])) 

   # -----------------------------------------------------------------------------------------------
   def plot_scatter(self, data_dico, value_names,
                    plot_filename, output_directory, plot_title='', pfilter=None,
                    pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
                    x_label=None, y_label=None, pnorm=None, pcmap='jet',
                    pcolor='b', pmarker='o', pmarkersize=5, 
                    pmarkerfacecolor=None, pmarkeredgecolor=None, 
                    is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None):

      """! Make a scatter plot """

      try:

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure
         iplot = 1

         if not plot_title is None:
            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = 'Scatter plot of %s' %(value_names)

            t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center',
                                 fontproperties=FontProperties(size=12))

         plt.grid(has_grid)
         if has_equalaxes:
            axes = fig.gca() 
            axes.set_aspect('equal')

         # Two values to represent
         if not value_names[0] in data_dico or not value_names[1] in data_dico:
            return
         values_x = numpy.asarray(data_dico[value_names[0]])
         values_y = numpy.asarray(data_dico[value_names[1]])

         if pfilter != None:

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

         #plot = plt.scatter(values_x, values_y, s=pmarkersize, c=pcolor, marker=pmarker, facecolor=pmarkerfacecolor, edgecolor=pmarkeredgecolor, norm=pnorm, cmap=pcmap)
         #plot = plt.scatter(values_x, values_y, s=pmarkersize,  marker=pmarker, facecolors=pmarkerfacecolor, norm=pnorm, cmap=pcmap)
         plot = plt.scatter(values_x, values_y, s=pmarkersize, c=pcolor, marker=pmarker, norm=pnorm, cmap=pcmap)

         if x_label is None:
            x_label = value_names[0]
         if y_label is None:
            y_label = value_names[1]

         plt.xlabel(x_label, fontproperties=FontProperties(size=11, weight='bold'))
         plt.ylabel(y_label, fontproperties=FontProperties(size=11, weight='bold'))

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

         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename))
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

         #rc('text', usetex=True)
         #rc('mathtext', fontset='custom')
         #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

         if pfigure == None:
            fig = plt.figure()
         else:
            fig = pfigure

         if not plot_title is None:
            if plot_title != '':
               figtitle = plot_title
            else:
               figtitle = r'Quiver plot of %s' %(value_names)

            t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center',
                                 fontproperties=FontProperties(size=12))

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

         #print "*** ", zip(values_x, values_y)

         #plt.quiver(x_coords, y_coords, values_x, values_y, scale=0.1, headwidth=0, width=0.01, pivot='middle')
         q_plot = plt.quiver(x_coords, y_coords, values_x, values_y, scale=scale, 
                             headwidth=0, width=width, pivot='middle', color=color)

         plt.xlabel(x_label, fontproperties=FontProperties(size=11, weight='bold'))
         plt.ylabel(y_label, fontproperties=FontProperties(size=11, weight='bold'))

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

#         plt.xlim(pxmin, pxmax)
#         plt.ylim(pymin, pymax)

         if show:
            plt.show()

         if pfigure is None:
            plt.savefig(os.path.join(output_dir, output_file))
            plt.clf()
            plt.close()

      except (OverflowError):        
            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
            #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))


#   # -----------------------------------------------------------------------------------------------
#   def plot_quiver(self, x_coords, y_coords, values_x, values_y, 
#                   output_file='quiver_plot', output_dir='.', plot_title='Quiver plot',
#                   pfilter=None, pxmin_lim=None, pxmax_lim=None, pymin_lim=None, pymax_lim=None, 
#                   x_label=None, y_label=None,
#                   color='#000080', scale=5, headwidth=0, width=0.001, pivot='middle',
#                   is_normed=False, has_equalaxes=False, has_grid=True, pfigure=None, show=False):
#      """! Make a quiver plot """

#      try:

#         rc('text', usetex=True)
#         rc('mathtext', fontset='custom')
#         rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

#         if pfigure == None:
#            fig = plt.figure()
#         else:
#            fig = pfigure

#         if not plot_title is None:
#            if plot_title != '':
#               figtitle = plot_title
#            else:
#               figtitle = r'Quiver plot of %s' %(value_names)

#            t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center', \
#                                 fontproperties=FontProperties(size=14, weight='bold'))

#         plt.grid(has_grid)
#         if has_equalaxes:
#            axes = fig.gca() 
#            axes.set_aspect('equal')

#         values_x = numpy.asarray(values_x)
#         values_y = numpy.asarray(values_y)

#         if pfilter != None:
#            filtered_values_x = values_x[pfilter]
#            filtered_values_y = values_y[pfilter]
#            if len(filtered_values_x) > 0:
#               values_x = filtered_values_x
#            if len(filtered_values_y) > 0:
#               values_y = filtered_values_y

#         #plt.quiver(x_coords, y_coords, values_x, values_y, scale=0.1, headwidth=0, width=0.01, pivot='middle')
#         q_plot = plt.quiver(x_coords, y_coords, values_x, values_y, scale=scale, headwidth=0, width=width, pivot='middle', color=color)

#         plt.xlabel(x_label, fontproperties=FontProperties(size=14, weight='bold'))
#         plt.ylabel(y_label, fontproperties=FontProperties(size=14, weight='bold'))

#         if pxmin_lim is None:
#            pxmin_lim = min(x_coords)
#         if pxmax_lim is None:
#            pxmax_lim = max(y_coords)
#         if pymin_lim is None:
#            pymin_lim = min(values_y)
#         if pymax_lim is None:
#            pymax_lim = max(values_y)

#         if not has_equalaxes:
#            pxmin = pxmin_lim
#            pymin = pymin_lim
#            pxmax = pxmax_lim
#            pymax = pymax_lim
#         else:
#            pxmin = pymin = min(pxmin_lim, pymin_lim) 
#            pxmax = pymax = max(pxmax_lim, pymax_lim) 

#         plt.xlim(pxmin, pxmax)
#         plt.ylim(pymin, pymax)

#         if show:
#            plt.show()

#         if pfigure is None:
#            plt.savefig(os.path.join(output_dir, output_file))
#            plt.clf()
#            plt.close()

#      except (OverflowError):        
#            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])     
#            #logger.log_error('Exception thrown during plotting (%s)' %(sys.exc_info()[1]))

#      return q_plot





# -- EOF sp_plot.py
