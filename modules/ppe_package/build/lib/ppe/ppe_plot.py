"""!
   ppe_plot.py - Plotting functions.
"""

# -- Python imports
import os, sys
import math
import numpy
import string
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("AGG")
from matplotlib.font_manager import FontProperties


# -------------------------------------------------------------------------------------------------
class PpePlotter(object):

   """!
      Class with a number of convenient plotting methods.
   """

   # ~~~~~~~~~~~~~~
   # Public methods
   # ~~~~~~~~~~~~~~

   # -----------------------------------------------------------------------------------------------
   def plot_stamp(self, stamp, plot_title='Stamp',
                        output_dir='.', output_file='stamp_plot', show=False):

      """! Plot a postage stamp
         @param stamp postage stamp to plot
         @param plot_title optional plot title
         @param output_dir optional output directory
         @param output_file optional output file name
         @param show True if the plot is shown once drawn, False otherwise
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
         plot = plt.pcolor(x, y, stamp, shading='flat')
         plt.colorbar(plot, shrink=1.0)

         plt.savefig(os.path.join(output_dir, output_file))

         if show:
            plt.show()

         plt.clf()
         plt.close()

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

                  #print "VAR:",value_names[iplot-1], "mean:", numpy.mean(values), "std:", numpy.std(values)

                  if pfilter != None:
                     filtered_values = values[pfilter]

                     if len(filtered_values) > 0:
                        values = filtered_values
                  color = pcolor
                  x_label = value_names[iplot-1]
                  y_label = ''

                  plt.xlabel(x_label, fontproperties=FontProperties(size=10, weight='regular'))

                  if pmin_xlim is None:
                     pmin_xlim = min(values)

                  if pmax_xlim is not None:
                     pmax_xlim = min(pmax_xlim, max(values))
                  else:
                     pmax_xlim = max(values)

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


#    # -----------------------------------------------------------------------------------------------
#    def plot_histogram(self, data_dico, value_names,
#                             plot_filename, output_directory, plot_title, pfilter=None,
#                             pnb_bins=15.0, pmin_xlim=None, pmax_xlim=None,
#                             x_labels=None, y_labels=None,
#                             pcolor='b', is_normed=False,
#                             has_equalaxes=False, has_grid=True, pfigure=None):
#
#       """! Plot histograms of possibly multiple arrays of floats found in a dictionary """
#
#       try:
#
#          if pfigure == None:
#             fig = plt.figure()
#          else:
#             fig = pfigure
#
#          iplot = 1
#
#          figtitle = plot_title
#
#          if pfigure == None:
#             t = plt.gcf().text(0.5, 0.95, figtitle, horizontalalignment='center',
#                                  fontproperties=FontProperties(size=12))
#
#          for name in value_names:
#
#             if value_names[iplot-1] in data_dico:
#
#                if len(value_names) > 1:
#                   if len(value_names) < 3:
#                      plt.subplot(2, 1, iplot)
#                   elif len(value_names) >= 3 and len(value_names) < 5:
#                      plt.subplot(2, 2, iplot)
#                   else:
#                      plt.subplot(3, 2, iplot)
#
#                plt.grid(has_grid)
#                if has_equalaxes:
#                   axes = fig.gca()
#                   axes.set_aspect('equal')
#
#                values = numpy.array(data_dico[value_names[iplot-1]], dtype=numpy.float64)
#                values = numpy.nan_to_num(values)
#
#                if pfilter != None:
#                   filtered_values = values[pfilter]
#                   if len(filtered_values) > 0:
#                      values = filtered_values
#                color = pcolor
#                if x_labels is None:
#                   x_label = value_names[iplot-1]
#                else:
#                   x_label = x_labels[iplot-1]
#
#                if y_labels is None:
#                   y_label = ''
#                else:
#                   y_label = y_labels[iplot-1]
#
#                plt.xlabel(x_label)
#                plt.ylabel(y_label)
#
#                if pmin_xlim is None:
#                   pmin_xlim = min(values)
#
#                if pmax_xlim is None:
#                   if max(values) > pmin_xlim:
#                      plt.xlim(pmin_xlim, max(values))
#                else:
#                   plt.xlim(pmin_xlim, pmax_xlim)
#
#                plt.hist(values, bins=pnb_bins, normed=is_normed, fc=pcolor, histtype='stepfilled')
#
#                iplot += 1
#
#          if pfigure is None:
#             plt.savefig(os.path.join(output_directory, plot_filename))
#
#             plt.clf()
#             plt.close()
#
#       finally:
#          plt.clf()
#          plt.close()


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
                                 fontproperties=FontProperties(size=11))

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
         plot = plt.scatter(values_x, values_y, s=pmarkersize, c=pcolor, marker=pmarker, norm=pnorm, cmap=pcmap, facecolor=pmarkerfacecolor, edgecolor=pmarkeredgecolor)

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

         if pxmax > pxmin:
            plt.xlim(pxmin, pxmax)
         if pymax > pymin:
            plt.ylim(pymin, pymax)

         if pfigure is None:
            plt.savefig(os.path.join(output_directory, plot_filename))
            plt.clf()
            plt.close()

      except (OverflowError):
            print 'Error: exception thrown during plotting (%s)' %(sys.exc_info()[1])


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

#   # -----------------------------------------------------------------------------------------------
#   def plot_stamp_mosaic(self, stamp_dico, coord_list, stamp_size, nb_rows, nb_cols,
#                               plot_title='Stamp Mosaic',
#                               output_directory='.', plot_filename='stamp_mosaic',
#                               has_grid=True, pfigure=None, show=False):

#      try:

#         if pfigure == None:
#            w, h = plt.figaspect(2.0)
#            fig = plt.figure(figsize=(w, h))
#         else:
#            fig = pfigure

#         plt.grid(has_grid)

#         iplot = 1

#         for (x, y) in coord_list:
#
##            irow = y/stamp_size+1
##            icol = x/stamp_size+1

##            print iplot, x, y, "irow:", irow, "icol:", icol, "nb_rows:", nb_rows, "nb_cols:", nb_cols,

#            #plt.subplot(irow, icol , iplot)

#            axes = fig.gca()
#            axes.set_aspect('equal')
#            axes.set_axis_off()

#            plt.subplot(nb_rows, nb_cols, iplot)

#            x_values = numpy.arange(0, stamp_size+1)
#            y_values = numpy.arange(0, stamp_size+1)

##            plot = plt.imshow(stamp_dico[(x, y)], aspect='auto', shape=(stamp_size, stamp_size))
#            plot = plt.imshow(stamp_dico[(x, y)], aspect='auto')

##            plot = plt.imshow(stamp_dico[(x, y)], aspect='equal',
##                              interpolation='nearest', cmap=matplotlib.cm.jet, origin='lower')
##            plot = plt.pcolor(x_values, y_values, stamp_dico[(x, y)], shading='flat')

#            iplot += 1

##         plt.xlim(1000)
##         plt.ylim(1000)

#         if pfigure is None:
#            plt.savefig(os.path.join(output_directory, plot_filename))
#            #plt.savefig(os.path.join(output_directory, plot_filename), bbox_inches='tight')

#            plt.clf()
#            plt.close()


#      finally:
#         plt.clf()
#         plt.close()


   # -----------------------------------------------------------------------------------------------
   def make_plots(self, job_result, plot_output_dir, master, plot_prefix=""):
         """! Plots histograms based on generated PSFEx catalog columns """

         job = job_result.job
         object_per_type_dico = job_result.result

         file_types = object_per_type_dico.keys()   # all results for each image types

         for file_type in file_types:

            object_dico = object_per_type_dico[file_type]

            # --- Where plots will be stored
            branch_tree = job.get_branch_tree()

            if master.logging_enabled():
               master.logger.log_info_p(
                 "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Generating plots...".format(
                 master, branch_tree, job.img_no, job.epoch,  file_type))

            max_lim = 150  # max limit threshold for second plot

            if len(job_result.data_dico) > 0:

               # --- Plot Histograms for each column in the catalog
               for col_name in job_result.data_dico[file_type].keys():
                  plot_filename = plot_prefix + \
                                  "hist_{0}_{1}_img_{2:03d}-{3:1d}_{4}.png".format(
                      col_name, "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)

                  if plot_prefix != "":
                     plot_title_prefix = plot_prefix.replace("_", " ").capitalize()
                  else:
                     plot_title_prefix = ""
                  plot_title = "{0}{1} - {2:03d}-{3:1d} - /{4} - {5}".format(
                                              plot_title_prefix, col_name,
                                              job.img_no, job.epoch, branch_tree, file_type)

                  self.plot_histogram(job_result.data_dico[file_type], [col_name],
                                      plot_filename, plot_output_dir, plot_title,
                                      pnb_bins=25.0, pmax_xlim=None, pcolor='b', is_normed=False)

                  if math.fabs(numpy.max(job_result.data_dico[file_type][col_name])) >= max_lim:
                     plot_filename = plot_prefix + \
                                     "hist_{0}_{1}_img_{2:03d}-{3:1d}_{4}_2.png".format(
                          col_name, "_".join(job.branch_dirs),
                          job.img_no, job.epoch, file_type)
                     self.plot_histogram(job_result.data_dico[file_type], [col_name],
                                         plot_filename, plot_output_dir, plot_title,
                                         pnb_bins=25.0, pmin_xlim=0, pmax_xlim=max_lim,
                                         pcolor='g', is_normed=False)


   # -----------------------------------------------------------------------------------------------
   def make_fragment_plots(self, fragment_data_dico, inner_radius, outer_radius,
                                 min_allowed_mag, max_allowed_mag,
                                 file_type, job, worker):
         """! Plot fragment per object related information """

         # --- Where plots will be stored
         branch_tree = job.get_branch_tree()
         plot_pathname = os.path.join(worker.plot_output_dir, branch_tree)

         if worker.logging_enabled():
            worker.logger.log_info_p(
              "{0} - /{1}/img-{2:03}-{3:1d} - {4} - Plotting fragment-related data...".format(
              worker, branch_tree, job.img_no, job.epoch,  file_type))

         # --- Fragment magnitude histogram plot
         plot_title = "Fragment Magnitude - {0:03d}-{1:1d} - /{2} - {3}".format(
                                  job.img_no, job.epoch, branch_tree, file_type)
         plot_filename = "hist_frag_mag_{0}_img_{1:03d}-{2:1d}_{3}.png".format(
             "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)

         self.plot_histogram(fragment_data_dico, ["mag"],
                             plot_filename, plot_pathname, plot_title,
                             x_labels=["magnitude"],
                             y_labels=["nb. fragments"],
                             pnb_bins=25.0, pmax_xlim=30, pcolor='#87cefa', is_normed=False)

         # --- Fragment distance histogram plot
         plot_title = "Fragment Distance/centroid - {0:03d}-{1:1d} - /{2} - {3}".format(
                                  job.img_no, job.epoch, branch_tree, file_type)
         plot_filename = "hist_frag_dist_{0}_img_{1:03d}-{2:1d}_{3}.png".format(
             "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)

         self.plot_histogram(fragment_data_dico, ["dist"],
                             plot_filename, plot_pathname, plot_title,
                             x_labels=["distance from centroid (pixels)"],
                             y_labels=["nb. fragments"],
                             pnb_bins=25.0, pmax_xlim=None, pcolor='#87cefa', is_normed=False)

         # --- Fragment count histogram plot
         plot_title = "Fragment Counts -{0:03d}-{1:1d} - /{2} - {3}".format(
                                  job.img_no, job.epoch, branch_tree, file_type)
         plot_filename = "hist_frag_count_{0}_img_{1:03d}-{2:1d}_{3}.png".format(
             "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)

         self.plot_histogram(fragment_data_dico, ["count"],
                             plot_filename, plot_pathname, plot_title,
                             x_labels=["fragment count"],
                             y_labels=None,
                             pnb_bins=25.0, pmax_xlim=None, pcolor='#87cefa', is_normed=False)

         # --- Fragment distance/magnitude scatter plot
         plot_title = "Fragment Distance/Magnitude - {0:03d}-{1:1d} - /{2} - {3}".format(
                                  job.img_no, job.epoch, branch_tree, file_type)
         plot_filename = "scatter_frag_{0}_img_{1:03d}-{2:1d}_{3}.png".format(
             "_".join(job.branch_dirs), job.img_no, job.epoch, file_type)
         fig = plt.figure()

         # --- Draw scatter plot
         max_mag = min(max(fragment_data_dico["mag"]), 30)
         self.plot_scatter(fragment_data_dico, ["mag", "dist"],
                    plot_filename, plot_pathname, plot_title=plot_title,
                    pfilter=None,
                    pxmin_lim=None, pxmax_lim=max_mag, pymin_lim=None, pymax_lim=None,
                    x_label="magnitude", y_label="distance from centroid (pixels)",
                    pnorm=None, pcmap='jet',  pcolor='b', pmarker='o', pmarkersize=5,
                    pmarkerfacecolor=None, pmarkeredgecolor=None,
                    is_normed=False, has_equalaxes=False, has_grid=True, pfigure=fig)

         # --- Draw distance bundary lines
         plt.plot([min(fragment_data_dico["mag"]), max(fragment_data_dico["mag"])], \
                  [inner_radius, inner_radius], color="r",
                  marker='', markersize=1.0, markerfacecolor='red', linestyle='--', linewidth=1)
         plt.plot([min(fragment_data_dico["mag"]), max(fragment_data_dico["mag"])], \
                  [outer_radius, outer_radius], color="g",
                  marker='', markersize=1.0, markerfacecolor='green', linestyle='--', linewidth=1)
         plt.plot([min_allowed_mag, min_allowed_mag], \
                  [min(fragment_data_dico["dist"]), max(fragment_data_dico["dist"])], color="b",
                  marker='', markersize=1.0, markerfacecolor='blue', linestyle='-', linewidth=1)
         plt.plot([max_allowed_mag, max_allowed_mag], \
                  [min(fragment_data_dico["dist"]), max(fragment_data_dico["dist"])], color="k",
                  marker='', markersize=1.0, markerfacecolor='black', linestyle='-', linewidth=1)

         plt.savefig(os.path.join(plot_pathname, plot_filename))
         plt.clf()
         plt.close()



# -- EOF ppe_plot.py
