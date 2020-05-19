# -*- coding: utf-8 -*-
"""SETOOLS SCRIPT

This script contain a class to handle operations on SExtractor output catalog.

:Authors: Axel Guinot

:Date: 16/01/2017

"""

import numpy as np
from math import ceil
import re
import operator
import string

import matplotlib.pylab as plt

import os

import shapepipe.pipeline.file_io as sc


def _mkdir(direc):
    """Create directory direc.
      This function should probably go somewhere further up.

    Parameters
    ----------
    direc: string
      Directory name

    Returns
    -------
    None
    """

    if not os.path.isdir(direc):
        try:
            os.mkdir(direc)
        except Exception as detail:
            raise Exception('SETools: Cannot create directory {}, error={}'.format(direc, detail))


class SETools(object):
    """SETools class

    Tools to analyse SExtractor catalogs.

    Parameters
    ----------
    cat: str, numpy.ndarray
        Path to SExtractor catalog (FITS_LDAC format) or numpy.ndarray (structure array)
    output_dir: str
        Path to pipeline result directory
    file_number_string: string
        input catalogue number/specifier
    config_filepath: str
        Path to config.setools file
    cat_file: boolean
        True if 'cat' is a path to a file. False otherwise

    """

    def __init__(self, cat, output_dir, file_number_string, config_filepath, cat_file=True):

        if cat_file:
            self._is_file = True
            self._cat_filepath = cat
            cat_file = sc.FITSCatalog(self._cat_filepath, SEx_catalog=True)
            cat_file.open()
            self._data = cat_file.get_data()
            cat_file.close()
        else:
            self._cat_filepath = None
            self._is_file = False
            self._data = cat
        self._cat_size = len(self._data)

        self._config_file = open(config_filepath)

        self._output_dir = output_dir

        self._file_number_string = file_number_string

    def process(self):
        """Process
        Main function called to process a SExtractor catalog.
        """

        if self._is_file:
            file_number = self._file_number_string
        else:
            file_number = ''

        self.read()

        # Processing: Create mask = filter input
        if len(self._mask) != 0:
            # MKDEBUG new, prevent to create directory, somehow this causes error sometimes
            # (in multi-process run?)
            direc = self._output_dir + '/mask'
            _mkdir(direc)
            # direc = self._output_dir
            self._make_mask()
            for i in self.mask.keys():
                if 'NO_SAVE' in self._mask[i]:
                    continue
                file_name = direc + '/' + i + file_number + '.fits'
                self.save_mask(self.mask[i], file_name)

        if len(self._plot) != 0:
            # MKDEBUG new, prevent to create directory, somehow this causes error sometimes
            # (in multi-process run?)
            # if self._plot_output_dir:
            #     direc = self._plot_output_dir
            # else:
            # Depreciated, should be removed, no mkdir here
            direc = self._output_dir + '/plot'
            _mkdir(direc)
            self._make_plot()
            for i in self.plot.keys():
                output_path = direc + '/' + i + file_number
                SEPlot(self.plot[i], self._data, output_path, self.mask)

        if len(self._new_cat) != 0:
            direc = self._output_dir + '/new_cat'
            _mkdir(direc)
            self._make_new_cat()
            for i in self.new_cat.keys():
                file_name = direc + '/' + i + file_number
                self.save_new_cat(self.new_cat[i], file_name)

        if len(self._rand_split) != 0:
            direc = self._output_dir + '/rand_split'
            _mkdir(direc)
            self._make_rand_split()
            for i in self.rand_split.keys():
                output_dir = direc + '/' + i + '_'
                self.save_rand_split(self.rand_split[i], output_dir, file_number)

        # if len(self._flag_split) != 0:
        #     direc = self._output_dir + '/flag_split'
        #     _mkdir(direc)
        #     self._make_flag_split()
        #     for i in self.flag_mask_dict.keys():
        #         output_path = direc + '/' + self._flag_split.keys()[0] + '_flag_' + i + file_number + '.fits'
        #         self.save_flag_split(self.flag_mask_dict[i], output_path)

        if len(self._stat) != 0:
            # MKDEBUG new, prevent to create directory, somehow this causes error sometimes
            # (in multi-process run?)
            # if self._stat_output_dir:
            #     direc = self._stat_output_dir
            # else:
            direc = self._output_dir + '/stat'
            _mkdir(direc)
            self._make_stat()
            for i in self.stat.keys():
                output_path = direc + '/' + i + file_number + '.txt'
                self.save_stat(self.stat[i], output_path)

    #####################
    # Reading functions #
    #####################

    def read(self):
        """Read the config file

        This function read the config file and create a dictionary for every task.

        """

        self._config_file.seek(0)

        self._mask = {}
        self._mask_key = []
        self._plot = {}
        self._stat = {}
        self._new_cat = {}
        self._rand_split = {}
        # self._flag_split={}
        in_section = 0
        while True:
            line_tmp = self._config_file.readline()

            if line_tmp == '':
                break

            line_tmp = self._clean_line(line_tmp)

            if line_tmp is None:
                continue

            # Loop over SETools file, look for section headers
            # [SECTION_TYPE:OBJECT_NAME], e.g.
            # [MASK:star_selection]

            if (in_section != 0) & (re.split(r'\[', line_tmp)[0] == ''):
                in_section = 0
            if not in_section:
                if (re.split(r'\[', line_tmp)[0] != ''):
                    raise Exception('No section found')

                sec = re.split(r'\[|\]', line_tmp)[1]
                if re.split(':', sec)[0] == 'MASK':
                    in_section = 1
                    try:
                        mask_name = re.split(':', sec)[1]
                    except:
                        mask_name = 'mask_{0}'.format(len(self._mask) + 1)
                    self._mask_key.append(mask_name)
                    self._mask[mask_name] = []
                elif re.split(':', sec)[0] == 'PLOT':
                    in_section = 2
                    try:
                        plot_name = re.split(':', sec)[1]
                    except:
                        plot_name = 'plot_{0}'.format(len(self._plot) + 1)
                    self._plot[plot_name] = []
                elif re.split(':', sec)[0] == 'STAT':
                    in_section = 3
                    try:
                        stat_name = re.split(':', sec)[1]
                    except:
                        stat_name = 'stat_{0}'.format(len(self._stat) + 1)
                    self._stat[stat_name] = []
                elif re.split(':', sec)[0] == 'NEW_CAT':
                    in_section = 4
                    try:
                        new_cat_name = re.split(':', sec)[1]
                    except:
                        new_cat_name = 'new_cat_{0}'.format(len(self._new_cat) + 1)
                    self._new_cat[new_cat_name] = []
                elif re.split(':', sec)[0] == 'RAND_SPLIT':
                    in_section = 5
                    try:
                        rand_split_name = re.split(':', sec)[1]
                    except:
                        rand_split_name = 'rand_split_{0}'.format(len(self._rand_split) + 1)
                    self._rand_split[rand_split_name] = []
                # elif re.split(':', sec)[0] == 'FLAG_SPLIT':
                #     in_section = 6
                #     try:
                #         flag_split_name = re.split(':', sec)[1]
                #     except:
                #         flag_split_name = 'flag_split_{0}'.format(len(self._flag_split)+1)
                #     self._flag_split[flag_split_name] = []
                else:
                    raise Exception("Section has to be in ['MASK','PLOT','STAT','NEW_CAT','RAND_SPLIT','FLAG_SPLIT']")
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
                # elif in_section == 6:
                #     self._flag_split[flag_split_name].append(line_tmp)

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

        if re.split('#', line_tmp)[0] == '':
            return None

        line_tmp = line_tmp.replace('\n', '')
        line_tmp = line_tmp.replace('\t', '')
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
        if len(mask) == 0:
            pass

        if output_path is None:
            raise ValueError('output path not provided')

        mask_file = sc.FITSCatalog(output_path, open_mode=sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=(self._cat_filepath is not None))
        mask_file.save_as_fits(data=self._data[mask], ext_name=ext_name, sex_cat_path=self._cat_filepath)

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
            new_file = sc.FITSCatalog(output_path + '.fits', open_mode=sc.BaseCatalog.OpenMode.ReadWrite)
            new_file.save_as_fits(data=new_cat, ext_name=ext_name)
        elif output_format == 'SEx_cat':
            new_file = sc.FITSCatalog(output_path + '.fits', open_mode=sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=(self._cat_filepath is not None))
            new_file.save_as_fits(data=new_cat, ext_name=ext_name, sex_cat_path=self._cat_filepath)
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
        data = self._data[mask]
        for i in rand_split.keys():
            rand_split_file = sc.FITSCatalog(output_path + i + file_number + '.fits', open_mode=sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog=(self._cat_filepath is not None))
            rand_split_file.save_as_fits(data=data[rand_split[i]], ext_name=ext_name, sex_cat_path=self._cat_filepath)

    # def save_flag_split(self, flag_split, output_path, ext_name='LDAC_OBJECTS'):
    #     """Save flag splitted catalogs

    #     Save catalogs following the flag split specified.

    #     Parameters
    #     ----------
    #     flag_split : dict
    #         Dictionary containing the indices for the split and mask to apply
    #     output_path : str
    #         Path of the output file
    #     ext_name : str
    #         Name of the extension where data are stored

    #     """

    #     if flag_split is None:
    #         raise ValueError('flag_split not provided')

    #     if output_path is None:
    #         raise ValueError('output path not provided')

    #     data = self._data[flag_split]
    #     flag_split_file = sc.FITSCatalog(output_path, open_mode= sc.BaseCatalog.OpenMode.ReadWrite, SEx_catalog= (self._cat_filepath is not None))
    #     flag_split_file.save_as_fits(data=data , ext_name=ext_name, sex_cat_path=self._cat_filepath)

    def save_stat(self, stat, output_path):
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

        for i in self._mask_key:
            global_mask = np.ones(self._cat_size, dtype=bool)
            global_ind = np.where(global_mask)[0]
            for j in self._mask[i]:
                s = re.split('{|}', j)
                if s[0] == '':
                    try:
                        global_mask &= self.mask[s[1]]
                        global_ind = np.where(global_mask)[0]
                        self._mask[i].pop(self._mask[i].index(j))
                    except:
                        raise ValueError("'{0}' not found in mask".format(s[1]))

            mask_tmp = None
            for j in self._mask[i]:
                if j == 'NO_SAVE':
                    continue
                tmp = sc.interpreter(j, self._data[global_mask], make_compare=True, mask_dict=self.mask).result
                if mask_tmp is None:
                    mask_tmp = np.ones(tmp.shape[0], dtype=bool)
                mask_tmp &= tmp

            new_ind = global_ind[mask_tmp]
            final_mask = np.zeros(self._cat_size, dtype=bool)
            final_mask[new_ind] = True

            self.mask[i] = final_mask

    def _make_plot(self):
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
                    raise ValueError('Plot option keyword/value not in correct format (key=val): {}'.format(j))
                ss = re.split('_', s[0])
                if len(ss) == 1:
                    self.plot[i][ss[0]] = {'0': s[1]}
                elif len(ss) == 2:
                    if ss[0] not in self.plot[i].keys():
                        self.plot[i][ss[0]] = {}
                    self.plot[i][ss[0]][ss[1]] = s[1]
                else:
                    raise ValueError('Plot keyword not in correct format (key or key_i): {}'.format(j))

    def _make_new_cat(self):
        """Make new catalog

        This function interprets the contents for each column of the new catalog.

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
                        self.new_cat[i][s[0]] = sc.interpreter(s[1], self._data, make_compare=False, mask_dict=self.mask).result
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
            n_keep = int(ceil(cat_size * ratio))
            mask_ratio = []
            mask_left = list(range(0, cat_size))
            while(len(mask_ratio) != n_keep):
                j = np.random.randint(0, len(mask_left))
                mask_ratio.append(mask_left.pop(j))
            mask_ratio = np.array(mask_ratio)
            mask_left = np.array(mask_left)
            self.rand_split[i]['mask'] = mask
            self.rand_split[i]['ratio_{0}'.format(int(ratio * 100))] = mask_ratio
            self.rand_split[i]['ratio_{0}'.format(100 - int(ratio * 100))] = mask_left

    # def _make_flag_split(self):
    #     """Make flag split

    #     This function create masks based on flags provided as a numpy extra file.

    #     """

    #     if len(self._flag_split) == 0:
    #         return None

    #     if self._extra_file is None:
    #         raise ValueError('An extra numpy file containing the flags has to provide for the flag split.')

    #     try:
    #         flags = np.load(self._extra_file, allow_pickle=True)
    #     except:
    #         raise ValueError("Can't open extra file : {}".format(self._extra_file))

    #     self.flag_mask_dict = {}
    #     for j in set(flags):
    #         self.flag_mask_dict[str(j)] = np.where(flags == j)

    def _make_stat(self):
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
                self.stat[i][s[0]] = sc.interpreter(s[1], self._data, make_compare=False, mask_dict=self.mask).result


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

    def __init__(self, plot_dict, catalog, output_path, mask_dict=None):

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
            self._check_key_for_plot(['X', 'Y'])
            # MKDEBUG: Following lines have been replaced by subroutine
            # if ('X' not in self._plot.keys()) | ('Y' not in self._plot.keys()):
            # raise ValueError('X and/or Y not provided')
            self._make_plot()
        elif self._plot['TYPE']['0'] in ['scatter', 'SCATTER']:
            self._check_key_for_plot(['X', 'Y'])
            # if ('X' not in self._plot.keys()) | ('Y' not in self._plot.keys()):
            # raise ValueError('X and/or Y not provided')
            self._make_scatter()
        elif self._plot['TYPE']['0'] in ['histogram', 'hist', 'HISTOGRAM', 'HIST']:
            self._check_key_for_plot(['Y'])
            # if 'Y' not in self._plot.keys():
            # raise ValueError('Y not provided')
            self._make_hist()
        else:
            ValueError('Type : {} not available'.format(self._plot['TYPE']['0']))

    def _check_key_for_plot(self, key_list):
        """Raise exception if keys not found in plot description.

        Parameters
        ----------
        key_list: list of strings
            list of keys

        Returns
        -------
        None
        """

        for key in key_list:
            if key not in self._plot.keys():
                raise ValueError('Key \'{}\' not provided for plot of type \'{}\''.format(key, self._plot['TYPE']['0']))

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
                    if ii % 2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
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
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
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
            if 'MARKERSIZE' in self._plot.keys():
                try:
                    markersize = self._plot['MARKERSIZE'][i]
                except:
                    markersize = 1
            else:
                markersize = 1
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

            plt.plot(sc.interpreter(x, self._cat, mask_dict=self._mask_dict).result,
                     sc.interpreter(self._plot['Y'][i], self._cat, mask_dict=self._mask_dict).result,
                     label=label, color=color, marker=marker, markersize=markersize, ls=line, alpha=alpha, figure=self._fig)

        # Set ploy limits for x and y
        for (lim, set_lim) in zip(['XLIM', 'YLIM'], [plt.xlim, plt.ylim]):
            if lim in self._plot.keys():
                try:
                    val = re.split(',', self._plot[lim]['0'])
                except:
                    raise ValueError('Plot {0} keyword/value not in correct format ({0}=lower,upper): {1}'.format(lim, self._plot[lim]['0']))

                set_lim(float(val[0]), float(val[1]))

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

        self._fig.savefig(self._output_path + '.' + out_format.lower(),
                          format=out_format)
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
                    if ii % 2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
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
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
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

            plt.scatter(sc.interpreter(x, self._cat, mask_dict=self._mask_dict).result,
                        sc.interpreter(y, self._cat, mask_dict=self._mask_dict).result,
                        c=sc.interpreter(self._plot['SCATTER'][i], self._cat, mask_dict=self._mask_dict).result,
                        label=label, marker=marker, alpha=alpha, figure=self._fig)

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

        self._fig.savefig(self._output_path + '.' + out_format.lower(), format=out_format)
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
                    if ii % 2 == 0:
                        title += i
                    else:
                        title += str(sc.interpreter(i, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
                    ii += 1
        else:
            title = ''

        self._fig.suptitle(title)

        if 'HTYPE' in self._plot.keys():
            htype = self._plot['HTYPE']['0']
        else:
            htype = 'bar'
        if 'LOG' in self._plot.keys():
            if self._plot['LOG']['0'] in ['True', 'true', '1']:
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
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(sc.interpreter(j, self._cat, make_compare=False, mask_dict=self._mask_dict).result)
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

            plt.hist(sc.interpreter(self._plot['Y'][i], self._cat, mask_dict=self._mask_dict).result,
                     bins=bins, color=color, label=label, alpha=alpha, histtype=htype, log=log)

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

        self._fig.savefig(self._output_path + '.' + out_format.lower(), format=out_format)
        plt.close()
