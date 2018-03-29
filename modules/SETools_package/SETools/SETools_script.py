# -*- coding: utf-8 -*-
"""SETOOLS SCRIPT

This script contain a class to handle operation on SExtractor output catalog.

:Authors: Axel Guinot

:Date: 16/01/2017

"""

# Compability with python2.x for x>6
from __future__ import print_function


import numpy as np
from math import ceil
import re
import operator
import string

import matplotlib
matplotlib.use("Agg")
import pylab as plt

import os

import scatalog as sc


class SETools(object):
    """SETools class

    Tools to analyse SExtractor catalogs.

    Parameters
    ----------
    cat_filepath: str
        Path to SExtractor catalog (FITS_LDAC format)
    config_filepath: str
        Path to config.setools file
    output_dir: str
        Path to pipeline result directory
    stat_output_dir: str, optional, default=None
        Path to pipeline statistics output directory
    plot_output_dir: str, optiona, default=None
        Path to pipeline plots output directory

    """

    def __init__(self, cat_filepath, config_filepath, output_dir, stat_output_dir=None, plot_output_dir=None):

        self._cat_filepath = cat_filepath
        self._config_filepath = config_filepath
        self._output_dir = output_dir

        # MK: The following dir names could be set to the default names if
        # not given as argument. However, I don't know how to retrieve them here
        # from _worker.
        if stat_output_dir:
            self._stat_output_dir = stat_output_dir
        if plot_output_dir:
            self._plot_output_dir = plot_output_dir

        self._cat_file = sc.FITSCatalog(self._cat_filepath, SEx_catalog=True)
        self._cat_file.open()
        self._cat_size = len(self._cat_file.get_data())

        self._config_file = open(self._config_filepath)


    #################
    # Main function #
    #################

    def process(self):
        """Process

        Main function called to process a SExtractor catalog.

        """

        s=re.split("\-([0-9]{3})\-([0-9]+)\.", self._cat_filepath)
        file_number = '-{0}-{1}'.format(s[1],s[2])

        self.read()

        ### Processing: Create mask = filter input
        if len(self._mask) != 0:
            if not os.path.isdir(self._output_dir + '/mask'):
                try:
                    os.system('mkdir {0}/mask'.format(self._output_dir))
                except:
                    raise Exception('Impossible to create the directory mask in {0}'.format(self._output_dir))
            self._make_mask()
            for i in self.mask.keys():
                if 'NO_SAVE' in self._mask[i]:
                    continue
                file_name = self._output_dir + '/mask/' + i + file_number + '.fits'
                self.save_mask(self.mask[i], file_name)

        if len(self._plot) != 0:
            if not os.path.isdir(self._output_dir + '/plot'):
                try:
                    os.system('mkdir {0}/plot'.format(self._output_dir))
                except:
                    raise Exception('Impossible to create the directory plot in {0}'.format(self._output_dir))
            self._make_plot2()
            for i in self.plot.keys():
                output_path = self._output_dir + '/plot/' + i + file_number
                plot_tmp = SEPlot(self.plot[i], self._cat_file.get_data(), output_path, self.mask)

        if len(self._new_cat) != 0:
            if not os.path.isdir(self._output_dir + '/new_cat'):
                try:
                    os.system('mkdir {0}/new_cat'.format(self._output_dir))
                except:
                    raise Exception('Impossible to create the directory new_cat in {0}'.format(self._output_dir))
            self._make_new_cat()
            for i in self.new_cat.keys():
                file_name = self._output_dir + '/new_cat/' + i + file_number
                self.save_new_cat(self.new_cat[i], file_name)

        if len(self._rand_split) != 0:
            if not os.path.isdir(self._output_dir + '/rand_split'):
                try:
                    os.system('mkdir {0}/rand_split'.format(self._output_dir))
                except:
                    raise Exception('Impossible to create the directory rand_split in {0}'.format(self._output_dir))
            self._make_rand_split()
            for i in self.rand_split.keys():
                output_dir = self._output_dir + '/rand_split/' + i + '_'
                self.save_rand_split(self.rand_split[i], output_dir, file_number)

        if len(self._stat) != 0:
            if not os.path.isdir(self._output_dir + '/stat'):
                try:
                    os.system('mkdir {0}/stat'.format(self._output_dir))
                except:
                    raise Exception('Impossible to create the directory stat in {0}'.format(self._output_dir))
            self._make_stat2()
            for i in self.stat.keys():
                output_path = self._output_dir + '/stat/' + i + file_number + '.txt'
                self.save_stat2(self.stat[i], output_path)


        ### Post-processing of output products
        # Statistics
        # self._make_stat()
        # for mask_type in self.stat.keys():
        #     file_name = self._stat_output_dir + '/' + mask_type + file_number + '.txt'
        #     self.save_stat(mask_type, file_name)

        # Plotting
        # self._make_plot(file_number)


    #####################
    # Reading functions #
    #####################

    def read(self):
        """Read the config file

        This function read the config file and create a dictionary for every task.

        """

        self._config_file.seek(0)

        self._mask={}
        self._plot={}
        self._stat={}
        self._new_cat={}
        self._rand_split={}
        in_section=0
        while True:
            line_tmp = self._config_file.readline()

            if line_tmp == '':
                break

            line_tmp = self._clean_line(line_tmp)

            if line_tmp == None:
                continue

            # Loop over SETools file, look for section headers
            # [SECTION_TYPE:OBJECT_TYPE], e.g.
            # [MASK:star_selection]

            if (in_section !=0) & (re.split('\[',line_tmp)[0] == ''):
                in_section = 0
            if not in_section:
                if (re.split('\[',line_tmp)[0] != ''):
                    raise Exception('No section found')

                sec=re.split('\[|\]', line_tmp)[1]
                if re.split(':', sec)[0] == 'MASK':
                    in_section = 1
                    try:
                        mask_name = re.split(':', sec)[1]
                    except:
                        mask_name = 'mask_{0}'.format(len(self._mask)+1)
                    self._mask[mask_name] = []
                elif re.split(':', sec)[0] == 'PLOT':
                    in_section = 2
                    try:
                        plot_name = re.split(':', sec)[1]
                    except:
                        plot_name = 'plot_{0}'.format(len(self._plot)+1)
                    self._plot[plot_name] = []
                elif re.split(':', sec)[0] == 'STAT':
                    in_section = 3
                    try:
                        stat_name = re.split(':', sec)[1]
                    except:
                        stat_name = 'stat_{0}'.format(len(self._stat)+1)
                    self._stat[stat_name] = []
                elif re.split(':', sec)[0] == 'NEW_CAT':
                    in_section = 4
                    try:
                        new_cat_name = re.split(':', sec)[1]
                    except:
                        new_cat_name = 'new_cat_{0}'.format(len(self._new_cat)+1)
                    self._new_cat[new_cat_name] = []
                elif re.split(':', sec)[0] == 'RAND_SPLIT':
                    in_section = 5
                    try:
                        rand_split_name = re.split(':', sec)[1]
                    except:
                        rand_split_name = 'rand_split_{0}'.format(len(self._rand_split)+1)
                    self._rand_split[rand_split_name] = []
                else:
                    raise Exception("Section has to be in ['MASK','PLOT','STAT','NEW_CAT','RAND_SPLIT']")
            else:
                if in_section == 1:
                    self._mask[mask_name].append(line_tmp)
                elif in_section == 2:
                    self._plot[plot_name].append(line_tmp)
                elif in_section == 3:
                    self._stat[stat_name].append(line_tmp)
                elif in_section == 4:
                    self._new_cat[new_cat_name].append(line_tmp)
                elif in_section == 5:
                    self._rand_split[rand_split_name].append(line_tmp)

    def _clean_line(self, line):
        """Clean Lines

        This function is called during the reading process to clean line
        from spaces, empty lines and ignore comments.

        Returns
        -------
        str
            If the line is not empty or a comment return the contents and
            None otherwise.

        """

        s = re.split('"', line)
        if len(s) == 3:
            line_tmp = s[0].replace(' ', '') + s[1] + s[2].replace(' ', '')
        else:
            line_tmp = line.replace(' ', '')

        if re.split('#',line_tmp)[0] == '':
            return None

        line_tmp = line_tmp.replace('\n','')
        line_tmp = line_tmp.replace('\t','')
        line_tmp = re.split('#', line_tmp)[0]

        if line_tmp != '':
            return line_tmp
        else:
            return None


    ##################
    # Save functions #
    ##################

    def save_mask(self, mask, output_path, ext_name='LDAC_OBJECTS'):
        """Save Mask

        This function will apply a mask on the data and save them into a new
        SExtractor catalog like fits file.

        Parameters
        ----------
        mask : numpy.ndarray
            Numpy array of boolean containing.
        output_path : str
            Path to the general output directory.
        ext_name : str, optional
            Name of the HDU containing masked data (default is 'LDAC_OBJECTS', SExtractor name)

        """

        if mask is None:
            raise ValueError('mask not provided')

        if output_path is None:
            raise ValueError('output path not provided')

        mask_file = sc.FITSCatalog(output_path, open_mode=sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=True)
        mask_file.save_as_fits(data=self._cat_file.get_data()[mask], ext_name=ext_name, sex_cat_path=self._cat_filepath)

    def save_new_cat(self, new_cat, output_path, ext_name='LDAC_OBJECTS'):
        """Save new catalog

        This function create a new catalog with a specific format
        (fits bin table, SExtractor like fits catalog or ASCII).

        Parameters
        ----------
        new_cat : dict
            Dictionary containing the data and "OUTPUT_FORMAT"
        output_path : str
            Path of the output
        ext_name : str
            Name of the extension for fits bin table output

        """

        try:
            output_format = new_cat.pop('OUTPUT_FORMAT')
        except:
            raise ValueError('OUTPUT_FORMAT not provided')

        if output_format == 'fits':
            new_file = sc.FITSCatalog(output_path + '.fits', open_mode= sc.BaseCatalog.OpenMode.ReadWrite)
            new_file.save_as_fits(data= new_cat, ext_name= ext_name)
        elif output_format == 'SEx_cat':
            new_file = sc.FITSCatalog(output_path + '.fits', open_mode= sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog= True)
            new_file.save_as_fits(data= new_cat, ext_name= ext_name, sex_cat_path= self._cat_filepath)
        elif (output_format == 'txt') | (output_format == 'ascii'):
            new_file = open(output_path + '.txt', 'w')
            new_file.write('# HEADER\n')
            new_file.write('# ')
            n_max = -1
            for i in new_cat.keys():
                if len(new_cat[i]) > n_max:
                    n_max = len(new_cat[i])
                new_file.write('{0}\t'.format(i))
            new_file.write('\n')
            for i in range(n_max):
                for j in new_cat.keys():
                    try:
                        new_file.write('{0}\t'.format(new_cat[j][i]))
                    except:
                        new_file.write('\t')
                new_file.write('\n')
            new_file.close()
        else:
            raise ValueError("Format should be in ['fits', 'SEx_cat', 'txt']")

    def save_rand_split(self, rand_split, output_path, file_number, ext_name='LDAC_OBJECTS'):
        """Save random splitted catalogs

        Save two catalogs following the random split specified.

        Parameters
        ----------
        rand_split : dict
            Dictionary containing the indices for the split and mask to apply
        output_path : str
            Path of the output dir
        file_number : str
            Numbering of the pipeline
        ext_name : str
            Name of the extension where data are stored

        """

        if rand_split is None:
            raise ValueError('rand_split not provided')

        if output_path is None:
            raise ValueError('output path not provided')
        if file_number is None:
            raise ValueError('file_number path not provided')

        mask = rand_split.pop('mask')
        data = self._cat_file.get_data()[mask]
        for i in rand_split.keys():
            rand_split_file = sc.FITSCatalog(output_path + i + file_number + '.fits', open_mode= sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=True)
            rand_split_file.save_as_fits(data=data[rand_split[i]], ext_name=ext_name, sex_cat_path=self._cat_filepath)

    def save_stat2(self, stat, output_path):
        """Save statistics

        Save the statistics in ASCII format.

        Parameters
        ----------
        stat : dict
            Dictionary containing the statistics
        output_path : str
            Path of the output dir

        """

        if stat is None:
            raise ValueError('stat not provided')
        if output_path is None:
            raise ValueError('output path not provided')

        f = open(output_path, 'w')
        f.write('# Statistics\n')

        for i in stat.keys():
            f.write(i + ' = ' + str(stat[i]) + '\n')

        f.close()


    def save_stat(self, mask_type, output_file):
        print('Writing stats to file {}'.format(output_file))
        f = open(output_file, 'w')
        print('# Statistics', file=f)

        for stat in self.stat[mask_type]:
            string = '{} = {}'.format(stat, self.stat[mask_type][stat])
            print(string, file=f)

        f.close()


    #####################
    # Bulding functions #
    #####################

    def _make_mask(self):
        """Make mask

        This function transforms the constraints from the config file to condition.

        """

        if len(self._mask) == 0:
            return None

        self.mask = {}
        for i in self._mask.keys():
            mask_tmp = np.ones(self._cat_size, dtype=bool)
            for j in self._mask[i]:
                if j == 'NO_SAVE':
                    continue
                mask_tmp &= sc.interpreter(j, self._cat_file.get_data(), make_compare= True, mask_dict= self.mask).result
            self.mask[i] = mask_tmp

    def _make_plot2(self):
        """Make plot

        This function interpret the different parameters for the plot.

        """

        if len(self._plot) == 0:
            return None

        self.plot = {}
        for i in self._plot.keys():
            self.plot[i] = {}
            for j in self._plot[i]:
                s = re.split('=', j)
                if len(s) != 2:
                    raise ValueError('Not a good format : {}'.format(j))
                ss = re.split('_', s[0])
                if len(ss) == 1:
                    self.plot[i][ss[0]] = {'0': s[1]}
                elif len(ss) == 2:
                    if ss[0] not in self.plot[i].keys():
                        self.plot[i][ss[0]] = {}
                    self.plot[i][ss[0]][ss[1]] = s[1]
                else:
                    raise ValueError('Not a good format : {}'.format(j))

    def _make_new_cat(self):
        """Make new catalog

        This function interpret the contents for each column of the new catalog.

        """

        if len(self._new_cat) == 0:
            return None

        self.new_cat = {}
        for i in self._new_cat.keys():
            self.new_cat[i] = {}
            for j in self._new_cat[i]:
                s = re.split('=', j)
                if len(s) == 2:
                    if s[0] == 'OUTPUT_FORMAT':
                        self.new_cat[i][s[0]] = s[1]
                    else:
                        self.new_cat[i][s[0]] = sc.interpreter(s[1], self._cat_file.get_data(), make_compare= False, mask_dict= self.mask).result
                else:
                    raise ValueError('Not a good format : {}'.format(j))

    def _make_rand_split(self):
        """Make random split

        This function create mask with random indices corresponding of the specfied ratio.

        """

        if len(self._rand_split) == 0:
            return None

        self.rand_split = {}
        mask = np.ones(self._cat_size, dtype=bool)
        for i in self._rand_split.keys():
            self.rand_split[i] = {}
            for j in self._rand_split[i]:
                s = re.split('=', j)
                if len(s) != 2:
                    raise ValueError('Not a good format : {}'.format(self._rand_split[i][0]))
                if s[0] == 'RATIO':
                    try:
                        ratio = float(s[1])
                    except:
                        raise ValueError('RATIO is not a number')
                    if ratio >= 1:
                        ratio /= 100.
                elif s[0] == 'MASK':
                    ss = re.split(',', s[1])
                    for k in ss:
                        try:
                            mask &= self.mask[k]
                        except:
                            raise ValueError('mask {0} does not exist'.format(k))

            cat_size = len(np.where(mask)[0])
            n_keep = int(ceil(cat_size*ratio))
            mask_ratio = []
            mask_left = range(0, cat_size)
            while(len(mask_ratio) != n_keep):
                j = np.random.randint(0, len(mask_left))
                mask_ratio.append(mask_left.pop(j))
            mask_ratio = np.array(mask_ratio)
            mask_left = np.array(mask_left)
            self.rand_split[i]['mask'] = mask
            self.rand_split[i]['ratio_{0}'.format(int(ratio*100))] = mask_ratio
            self.rand_split[i]['ratio_{0}'.format(100-int(ratio*100))] = mask_left

    def _make_stat2(self):
        """Make statistics

        This function interpret the different statistics required.

        """

        if len(self._stat) == 0:
            return None

        self.stat = {}
        for i in self._stat.keys():
            self.stat[i] = {}
            for j in self._stat[i]:
                s = re.split('=', j)
                if len(s) != 2:
                    raise ValueError('Not a good format : {}'.format(j))
                self.stat[i][s[0]] = sc.interpreter(s[1], self._cat_file.get_data(), make_compare= False, mask_dict= self.mask).result


    def _make_plot(self, file_number):
        """Produce plots, and writes them to disk as pdf files.
           TODO: Move saving of plots outside this method.

        Parameters
        ----------
        file_number: int
            running image number

        Returns
        -------
        None
        """

        self.plot = {}
        # Loop over mask types (different selections)
        for mask_type in self._plot:

            # Loop over lines in plot section, make plot for each line
            for plot_line in self._plot[mask_type]:

                fig, ax = self.create_new_plot()
                plot_type, plot_x, plot_y, cmds = self.get_plot_info_from_line(plot_line)

                if mask_type == 'ALL':

                    # Special mask type (plot all selections in one plot). Add to plot
                    # all objects (no mask), and all masks in config file
                    self.fill_plot_from_mask(ax, 'ALL', plot_type, plot_x, plot_y)
                    for mask_type_for_all in self._mask:
                        self.fill_plot_from_mask(ax, mask_type_for_all, plot_type, plot_x, plot_y)

                else:

                    self.fill_plot_from_mask(ax, mask_type, plot_type, plot_x, plot_y)

                self.finalize_plot(ax, plot_x, plot_y, cmds=cmds)
                self.save_plot(mask_type, '{}:{}:{}'.format(plot_type, plot_x, plot_y), file_number)

    def get_plot_info_from_line(self, plot_line):
        """Return plot information from line in PLOT section

        Parameters
        ----------
        plot_line: ':'-separated string
            plot information (line in config file)

        Returns
        -------
        plot_type: string
            plot type
        plot_x: string
            column name to be plotted as x axis
        plot_y: string
            column name to be plotted as y axis
        cmds: string
            extra plt commands to be run (can None)
        """

        # The following can be done with better parsing...
        line = re.split(':', plot_line)
        if len(line) < 3:
            raise Exception('Plot definition needs to have at least three entries (type, x, y), found {}'.format(len(line)))

        plot_type = line[0]
        plot_x    = line[1]
        plot_y    = line[2]
        if len(line) > 3:
            cmds = line[3]
        else:
            cmds = []

        return plot_type, plot_x, plot_y, cmds

    def create_new_plot(self):
        """Create a new plot.

        Parameters
        ----------
        None

        Returns
        -------
        fig: figure object
            current figure
        ax: axes object
            current axes
        """

        fig = plt.clf()
        ax  = plt.gca()

        return fig, ax

    def fill_plot_from_mask(self, ax, mask_type, plot_type, plot_x, plot_y):
        """Fill plot with data from mask/selection.

        Parameters
        ----------
        ax: axes object
            current axes
        mask_type: string
            mask (selection) type, define section in config file,
            'ALL' for no mask (select all objects)
        plot_type: string
            plot type
        plot_x: string
            column name to be plotted as x axis
        plot_y: string
            column name to be plotted as y axis

        Returns
        -------
        None
        """

        # Get objects in mask
        if mask_type == 'ALL':
            dat = self._cat_file.get_data()
        else:
            dat = self._cat_file.get_data()[self.mask[mask_type]]

        x   = dat[plot_x]
        y   = dat[plot_y]

        # Plot
        if plot_type == 'scatter':
            ax.scatter(x, y, 1, label=mask_type)
            #plt.plot(x, y, 'o', markersize=3, markeredgewidth=0.3, markerfacecolor='none', label=mask_type)
        else:
            raise Exception('Plot type \'{}\' not supported'.format(plot_type))

    def finalize_plot(self, ax, plot_x, plot_y, cmds=[]):
        """Performe cosmetics on plot at the end, should be called once, before saving plot.

        Parameters
        ----------
        plot_x: string
            column name to be plotted as x axis
        plot_y: string
            column name to be plotted as y axis
        cmds: string, optiona, default=[]
            extra plt commands

        Returns
        -------
        None
        """

        # TODO: How do we get the units?
        plt.xlabel(plot_x)
        plt.ylabel(plot_y)
        plt.legend()

        if len(cmds) > 0:
            cmd_list = re.split(';', cmds)
            for cmd in cmd_list:
                exec(cmd)

    def save_plot(self, mask_type, plot_name, file_number):
        """Save plot to file.

        Parameters
        ----------
        mask_type: string
            mask (selection) type
        plot_name: string
            plot information (line in config file)
        file_number: int
            running image number

        Returns
        -------
        None
        """

        # TODO: check if plot_name already exists, if yes create unique name
        out_name = '{}/{}-{}{}.pdf'.format(self._plot_output_dir, mask_type, plot_name, file_number)
        print('Saving plot to file {}'.format(out_name))
        plt.savefig('{}'.format(out_name))

    def _make_stat(self):
        """Compute statistics on output products.
           Fill self.stat with results of stats computations.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.stat = {}
        # Loop over mask types (different selections)
        for mask_type in self._stat:

            self.stat[mask_type] = {}

            # Loop over statistics (lines in stat section)
            for stat in self._stat[mask_type]:

                if stat == 'NUMBER_OBJ':
                    # Number of selected objects
                    res = len(self._cat_file.get_data()[self.mask[mask_type]])
                    self.stat[mask_type][stat] = res
                else:
                    print('Warning: statistic \'{}\' not supported, continuing...'.format(stat))


class SEPlot(object):
    """SEPlot class

    Tools to create plots.

    Parameters
    ----------
    plot_dict : dict
        Dictionary containing the parameters for the plot
    catalog : numpy.recarray or astropy.fits.fitsrec
        Array containing the full data
    output_path : str
        Path for the output
    mask_dict : dict
        Dictionary containing masks to apply

    Notes
    -----

    Types of plots available : plot, scatter, hist.

    """

    def __init__(self, plot_dict, catalog, output_path, mask_dict = None):

        if plot_dict is None:
            raise ValueError('plot_dict not provided')
        if catalog is None:
            raise ValueError('catalog not provided')
        if output_path is None:
            raise ValueError('output_path not provided')

        self._plot = plot_dict
        self._output_path = output_path
        self._cat = catalog
        self._mask_dict = mask_dict

        if 'TYPE' not in self._plot.keys():
            raise ValueError('Plot type not specified')

        if self._plot['TYPE']['0'] in ['plot', 'PLOT']:
            if ('X' not in self._plot.keys()) | ('Y' not in self._plot.keys()):
                raise ValueError('X and/or Y not provided')
            self._make_plot()
        elif self._plot['TYPE']['0'] in ['scatter', 'SCATTER']:
            if ('X' not in self._plot.keys()) | ('Y' not in self._plot.keys()):
                raise ValueError('X and/or Y not provided')
            self._make_scatter()
        elif self._plot['TYPE']['0'] in ['histogram', 'hist', 'HISTOGRAM', 'HIST']:
            if 'Y' not in self._plot.keys():
                raise ValueError('Y not provided')
            self._make_hist()
        else:
            ValueError('Type : {} not available'.format(self._plot['TYPE']['0']))


    def _make_plot(self):
        """Make plot

        This function call pyplot.plot.

        """

        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                ii = 1
                for i in s[1:-1]:
                    if ii%2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                    ii += 1
        else:
            title = ''

        self._fig.suptitle(title)

        for i in self._plot['Y'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][i]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj%2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                            jj += 1
                except:
                    label = None
            else:
                label = None
            if 'COLOR' in self._plot.keys():
                try:
                    color = self._plot['COLOR'][i]
                except:
                    color = None
            else:
                color = None
            if 'MARKER' in self._plot.keys():
                try:
                    marker = self._plot['MARKER'][i]
                except:
                    marker = '+'
            else:
                marker = '+'
            if 'LINE' in self._plot.keys():
                try:
                    line = self._plot['LINE'][i]
                except:
                    line = ''
            else:
                line = ''
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = self._plot['ALPHA'][i]
                except:
                    alpha = None
            else:
                alpha = None

            try:
                x = self._plot['X'][i]
            except:
                if len(self._plot['X']) == 1:
                    x = self._plot['X'][self._plot['X'].keys()[0]]
                else:
                    raise ValueError("You need to specified X for each Y provided if they dont have the same")

            plt.plot(sc.interpreter(x, self._cat, mask_dict= self._mask_dict).result,
                     sc.interpreter(self._plot['Y'][i], self._cat, mask_dict= self._mask_dict).result,
                     label= label, color= color, marker= marker, ls= line, alpha= alpha, figure= self._fig)

        if 'LABEL' in self._plot.keys():
            plt.legend()

        if 'XLABEL' in self._plot.keys():
            plt.xlabel(self._plot['XLABEL']['0'])
        if 'YLABEL' in self._plot.keys():
            plt.ylabel(self._plot['YLABEL']['0'])

        if 'FORMAT' in self._plot.keys():
            out_format = self._plot['FORMAT']['0']
        else:
            out_format = "PNG"

        self._fig.savefig(self._output_path + '.' + out_format.lower(), format= out_format)
        plt.close()


    def _make_scatter(self):
        """Make scatter

        This function call pyplot.scatter.

        """

        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                ii = 1
                for i in s[1:-1]:
                    if ii%2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                    ii += 1
        else:
            title = ''

        self._fig.suptitle(title)

        for i in self._plot['SCATTER'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][i]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj%2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                            jj += 1
                except:
                    label = None
            else:
                label = None
            if 'MARKER' in self._plot.keys():
                try:
                    marker = self._plot['MARKER'][i]
                except:
                    marker = '+'
            else:
                marker = '+'
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = self._plot['ALPHA'][i]
                except:
                    alpha = None
            else:
                alpha = None

            try:
                x = self._plot['X'][i]
            except:
                if len(self._plot['X']) == 1:
                    x = self._plot['X'][self._plot['X'].keys()[0]]
                else:
                    raise ValueError("You need to specified X for each SCATTER provided if they dont have the same")
            try:
                y = self._plot['Y'][i]
            except:
                if len(self._plot['Y']) == 1:
                    y = self._plot['Y'][self._plot['Y'].keys()[0]]
                else:
                    raise ValueError("You need to specified Y for each SCATTER provided if they dont have the same")

            plt.scatter(sc.interpreter(x, self._cat, mask_dict= self._mask_dict).result,
                        sc.interpreter(y, self._cat, mask_dict= self._mask_dict).result,
                        c = sc.interpreter(self._plot['SCATTER'][i], self._cat, mask_dict= self._mask_dict).result,
                        label= label, marker= marker, alpha= alpha, figure= self._fig)

        if 'LABEL' in self._plot.keys():
            plt.legend()

        if 'XLABEL' in self._plot.keys():
            plt.xlabel(self._plot['XLABEL']['0'])
        if 'YLABEL' in self._plot.keys():
            plt.ylabel(self._plot['YLABEL']['0'])

        plt.colorbar()

        if 'FORMAT' in self._plot.keys():
            out_format = self._plot['FORMAT']['0']
        else:
            out_format = "PNG"

        self._fig.savefig(self._output_path + '.' + out_format.lower(), format= out_format)
        plt.close()


    def _make_hist(self):
        """Make hist

        This function call pyplot.hist.

        """

        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                ii = 1
                for i in s[1:-1]:
                    if ii%2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                    ii += 1
        else:
            title = ''

        self._fig.suptitle(title)

        if 'HTYPE' in self._plot.keys():
            htype = self._plot['HTYPE']['0']
        else:
            htype = 'bar'
        if 'LOG' in self._plot.keys():
            if (self._plot['LOG']['0'] == 'True') | (self._plot['LOG']['0'] == 'true') | (self._plot['LOG']['0'] == '1'):
                log = True
            else:
                log = False
        else:
            log = False

        for i in self._plot['Y'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][i]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj%2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare= False, mask_dict= self._mask_dict).result)
                            jj += 1
                except:
                    label = None
            else:
                label = None
            if 'COLOR' in self._plot.keys():
                try:
                    color = self._plot['COLOR'][i]
                except:
                    color = None
            else:
                color = None
            if 'BIN' in self._plot.keys():
                try:
                    bins = int(self._plot['BIN'][i])
                except:
                    if len(self._plot['BIN']) == 1:
                        bins = int(self._plot['BIN'][self._plot['BIN'].keys()[0]])
            else:
                bins = 50
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = float(self._plot['ALPHA'][i])
                except:
                    alpha = None
            else:
                alpha = None

            plt.hist(sc.interpreter(self._plot['Y'][i], self._cat, mask_dict= self._mask_dict).result,
                     bins= bins, color= color, label= label, alpha= alpha, histtype= htype, log= log)

        if 'LABEL' in self._plot.keys():
            plt.legend()

        if 'XLABEL' in self._plot.keys():
            plt.xlabel(self._plot['XLABEL']['0'])
        if 'YLABEL' in self._plot.keys():
            plt.ylabel(self._plot['YLABEL']['0'])


        if 'FORMAT' in self._plot.keys():
            out_format = self._plot['FORMAT']['0']
        else:
            out_format = "PNG"

        self._fig.savefig(self._output_path + '.' + out_format.lower(), format= out_format)
        plt.close()
