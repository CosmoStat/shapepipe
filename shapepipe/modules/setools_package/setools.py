"""SETOOLS.

This module contains a class to handle operations on SExtractor output
catalogues.

:Authors: Axel Guinot

"""

import operator
import os
import re
import string

import matplotlib.pylab as plt
import numpy as np

from shapepipe.pipeline import file_io
from shapepipe.pipeline.str_handler import StrInterpreter
from shapepipe.utilities.file_system import mkdir


class SETools(object):
    """The SETools Class.

    Tools to analyse SExtractor catalogues.

    Parameters
    ----------
    cat: str, numpy.ndarray
        Path to SExtractor catalogue (FITS_LDAC format) or numpy.ndarray
        (structured array)
    output_dir: str
        Path to pipeline result directory
    file_number_string: str
        input catalogue number/specifier
    config_filepath: str
        Path to ``config.setools`` file
    cat_file: bool, optional
        ``True`` if ``cat`` is a path to a file, ``False`` otherwise

    """

    def __init__(
        self,
        cat,
        output_dir,
        file_number_string,
        config_filepath,
        cat_file=True,
    ):

        if cat_file:
            self._is_file = True
            self._cat_filepath = cat
            cat_file = file_io.FITSCatalogue(
                self._cat_filepath,
                SEx_catalogue=True,
            )
            cat_file.open()
            try:
                self._data = cat_file.get_data()
            except:
                raise IOError(f"Could not load catalogue data from {cat}")
            cat_file.close()

        else:
            self._cat_filepath = None
            self._is_file = False
            self._data = cat

        self._cat_size = len(self._data)
        self._config_file = open(config_filepath)
        self._output_dir = output_dir
        self._file_number_string = file_number_string

    def process(self, w_log):
        """Process.

        Main function called to process a SExtractor catalogue.

        Parameters
        ----------
        w_log : logging.Logger
            Logging instance

        """
        if self._is_file:
            file_number = self._file_number_string
        else:
            file_number = ''

        self.read()

        # Processing: Create mask = filter input
        if len(self._mask) != 0:
            direc = f'{self._output_dir}/mask'
            mkdir(direc)
            self._make_mask()
            for key in self.mask.keys():
                if 'NO_SAVE' in self._mask[key]:
                    continue
                file_name = f'{direc}/{key}{file_number}.fits'
                self.save_mask(self.mask[key], file_name)

        if len(self._plot) != 0:
            direc = f'{self._output_dir}/plot'
            mkdir(direc)
            self._make_plot()
            for key in self.plot.keys():
                output_path = f'{direc}/{key}{file_number}'
                SEPlot(self.plot[key], self._data, output_path, self.mask)

        if len(self._new_cat) != 0:
            direc = f'{self._output_dir}/new_cat'
            mkdir(direc)
            self._make_new_cat()
            for key in self.new_cat.keys():
                file_name = f'{direc}/{key}{file_number}'
                self.save_new_cat(self.new_cat[key], file_name)

        if len(self._rand_split) != 0:
            direc = f'{self._output_dir}/rand_split'
            mkdir(direc)
            self._make_rand_split()

            for sample_type in self.rand_split.keys():
                empty_found = False
                for ratio in self.rand_split[sample_type].keys():
                    if (
                        not empty_found
                        and len(self.rand_split[sample_type][ratio]) == 0
                    ):
                        empty_found = True
                if empty_found:
                    w_log.info(
                        'At least one random-split catalogue is empty, no '
                        + 'random sub-samples written for sample_type='
                        + f'{sample_type}')
                    continue

                output_dir = f'{direc}/{sample_type}_'
                self.save_rand_split(
                    self.rand_split[sample_type],
                    output_dir,
                    file_number,
                )

        if len(self._stat) != 0:
            direc = f'{self._output_dir}/stat'
            mkdir(direc)
            self._make_stat()
            for key in self.stat.keys():
                output_path = f'{direc}/{key}{file_number}.txt'
                self.save_stat(self.stat[key], output_path)

    def read(self):
        """Read the Configuration File.

        This function read the config file and creates a dictionary for every
        task.

        Raises
        ------
        RuntimeError
            If not section found
        ValueError
            For invalid section

        """
        self._config_file.seek(0)

        self._mask = {}
        self._mask_key = []
        self._plot = {}
        self._stat = {}
        self._new_cat = {}
        self._rand_split = {}
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
                    raise RuntimeError('No section found')

                sec = re.split(r'\[|\]', line_tmp)[1]
                if re.split(':', sec)[0] == 'MASK':
                    in_section = 1
                    try:
                        mask_name = re.split(':', sec)[1]
                    except Exception:
                        mask_name = f'mask_{len(self._mask) + 1}'
                    self._mask_key.append(mask_name)
                    self._mask[mask_name] = []
                elif re.split(':', sec)[0] == 'PLOT':
                    in_section = 2
                    try:
                        plot_name = re.split(':', sec)[1]
                    except Exception:
                        plot_name = f'plot_{len(self._plot) + 1}'
                    self._plot[plot_name] = []
                elif re.split(':', sec)[0] == 'STAT':
                    in_section = 3
                    try:
                        stat_name = re.split(':', sec)[1]
                    except Exception:
                        stat_name = f'stat_{len(self._stat) + 1}'
                    self._stat[stat_name] = []
                elif re.split(':', sec)[0] == 'NEW_CAT':
                    in_section = 4
                    try:
                        new_cat_name = re.split(':', sec)[1]
                    except Exception:
                        new_cat_name = f'new_cat_{len(self._new_cat) + 1}'
                    self._new_cat[new_cat_name] = []
                elif re.split(':', sec)[0] == 'RAND_SPLIT':
                    in_section = 5
                    try:
                        rand_split_name = re.split(':', sec)[1]
                    except Exception:
                        rand_split_name = (
                            f'rand_split_{len(self._rand_split) + 1}'
                        )
                    self._rand_split[rand_split_name] = []
                else:
                    raise ValueError(
                        "Section has to be in ['MASK','PLOT','STAT',"
                        + "'NEW_CAT','RAND_SPLIT','FLAG_SPLIT']"
                    )
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
        """Clean Lines.

        This function is called during the reading process to clean lines
        from spaces, empty lines and ignore comments.

        Parameters
        ----------
        line : str
            Input line

        Returns
        -------
        str
            If the line is not empty or a comment return the contents and
            ``None`` otherwise

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

    def save_mask(self, mask, output_path, ext_name='LDAC_OBJECTS'):
        """Save Mask.

        This function will apply a mask to the data and save them into a new
        SExtractor-like FITS file.

        Parameters
        ----------
        mask : numpy.ndarray
            Numpy array of boolean containing
        output_path : str
            Path to the general output directory
        ext_name : str, optional
            Name of the HDU containing masked data; default is ``LDAC_OBJECTS``
            SExtractor name

        Raises
        ------
        ValueError
            If mask not provided
        ValueError
            If output path not provided

        """
        if mask is None:
            raise ValueError('mask not provided')
        if len(mask) == 0:
            pass

        if output_path is None:
            raise ValueError('output path not provided')

        mask_file = file_io.FITSCatalogue(
            output_path,
            open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            SEx_catalogue=(self._cat_filepath is not None),
        )
        mask_file.save_as_fits(
            data=self._data[mask],
            ext_name=ext_name,
            sex_cat_path=self._cat_filepath,
        )

    def save_new_cat(self, new_cat, output_path, ext_name='LDAC_OBJECTS'):
        """Save New Catalogue.

        This function creates a new catalogue with a specific format
        (FITS bin table, SExtractor-like FITS catalogue or ASCII).

        Parameters
        ----------
        new_cat : dict
            Dictionary containing the data and ``OUTPUT_FORMAT``
        output_path : str
            Path of the output
        ext_name : str
            Name of the extension for FITS bin table output

        Raises
        ------
        ValueError
            If ``OUTPUT_FORMAT`` not provided
        ValueError
            For invalid output format

        """
        try:
            output_format = new_cat.pop('OUTPUT_FORMAT')
        except Exception:
            raise ValueError('OUTPUT_FORMAT not provided')

        if output_format == 'fits':
            new_file = file_io.FITSCatalogue(
                output_path + '.fits',
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
            )
            new_file.save_as_fits(data=new_cat, ext_name=ext_name)
        elif output_format == 'SEx_cat':
            new_file = file_io.FITSCatalogue(
                output_path + '.fits',
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
                SEx_catalogue=(self._cat_filepath is not None),
            )
            new_file.save_as_fits(
                data=new_cat,
                ext_name=ext_name,
                sex_cat_path=self._cat_filepath,
            )
        elif (output_format == 'txt') | (output_format == 'ascii'):
            new_file = open(output_path + '.txt', 'w')
            new_file.write('# HEADER\n')
            new_file.write('# ')
            n_max = -1
            for key in new_cat.keys():
                if len(new_cat[key]) > n_max:
                    n_max = len(new_cat[key])
                new_file.write(f'{key}\t')
            new_file.write('\n')
            for idx in range(n_max):
                for key in new_cat.keys():
                    try:
                        new_file.write(f'{new_cat[key][idx]}\t')
                    except Exception:
                        new_file.write('\t')
                new_file.write('\n')
            new_file.close()
        else:
            raise ValueError("Format should be in ['fits', 'SEx_cat', 'txt']")

    def save_rand_split(
        self,
        rand_split,
        output_path,
        file_number,
        ext_name='LDAC_OBJECTS',
    ):
        """Save Random Split Catalogues.

        Save two catalogues following the random split specified.

        Parameters
        ----------
        rand_split : dict
            Dictionary containing the indices for the split and mask to apply
        output_path : str
            Path of the output dir
        file_number : str
            Numbering of the pipeline
        ext_name : str, optional
            Name of the extension where data are stored

        Raises
        ------
        ValueError
            If ``rand_split`` not provided
        ValueError
            If ``output_path`` not provided
        ValueError
            If ``file_number`` not provided

        """
        if rand_split is None:
            raise ValueError('rand_split not provided')
        if output_path is None:
            raise ValueError('output path not provided')
        if file_number is None:
            raise ValueError('file_number path not provided')

        mask = rand_split.pop('mask')
        data = self._data[mask]

        for idx in rand_split.keys():
            rand_split_file = file_io.FITSCatalogue(
                f'{output_path}{idx}{file_number}.fits',
                open_mode=file_io.BaseCatalogue.OpenMode.ReadWrite,
                SEx_catalogue=(self._cat_filepath is not None),
            )
            rand_split_file.save_as_fits(
                data=data[rand_split[idx]],
                ext_name=ext_name,
                sex_cat_path=self._cat_filepath,
            )

    def save_stat(self, stat, output_path):
        """Save Statistics.

        Save the statistics in ASCII format.

        Parameters
        ----------
        stat : dict
            Dictionary containing the statistics
        output_path : str
            Path of the output dir

        Raises
        ------
        ValueError
            If ``stat`` not provided
        ValueError
            If ``output`` path not provided

        """
        if stat is None:
            raise ValueError('stat not provided')
        if output_path is None:
            raise ValueError('output path not provided')

        file = open(output_path, 'w')
        file.write('# Statistics\n')

        for key in stat.keys():
            file.write(f'{key} = {str(stat[key])}\n')

        file.close()

    def _make_mask(self):
        """Make Mask.

        This function transforms the constraints from the config file to
        conditions.

        Raises
        ------
        RuntimeError
            If value not found in mask

        """
        if len(self._mask) == 0:
            return None

        self.mask = {}

        for key in self._mask_key:
            global_mask = np.ones(self._cat_size, dtype=bool)
            global_ind = np.where(global_mask)[0]
            for idx in self._mask[key]:
                s = re.split('{|}', idx)
                if s[0] == '':
                    try:
                        global_mask &= self.mask[s[1]]
                        global_ind = np.where(global_mask)[0]
                        self._mask[key].pop(self._mask[key].index(idx))
                    except Exception:
                        raise RuntimeError(f'"{s[1]}" not found in mask')

            mask_tmp = None
            for idx in self._mask[key]:
                if idx == 'NO_SAVE':
                    continue
                tmp = StrInterpreter(
                    idx,
                    self._data[global_mask],
                    make_compare=True,
                    mask_dict=self.mask,
                ).result
                if mask_tmp is None:
                    mask_tmp = np.ones(tmp.shape[0], dtype=bool)
                mask_tmp &= tmp

            new_ind = global_ind[mask_tmp]
            final_mask = np.zeros(self._cat_size, dtype=bool)
            final_mask[new_ind] = True

            self.mask[key] = final_mask

    def _make_plot(self):
        """Make Plot.

        This function interpret the different parameters for the plot.

        Raises
        ------
        ValueError
            If plot option keyword or value is not in the correct format
        ValueError
            If plot keyword not in the correct format

        """
        if len(self._plot) == 0:
            return None

        self.plot = {}
        for key in self._plot.keys():
            self.plot[key] = {}
            for idx in self._plot[key]:
                s = re.split('=', idx)
                if len(s) != 2:
                    raise ValueError(
                        'Plot option keyword/value not in correct format '
                        + f'(key=val): {idx}'
                    )
                ss = re.split('_', s[0])
                if len(ss) == 1:
                    self.plot[key][ss[0]] = {'0': s[1]}
                elif len(ss) == 2:
                    if ss[0] not in self.plot[key].keys():
                        self.plot[key][ss[0]] = {}
                    self.plot[key][ss[0]][ss[1]] = s[1]
                else:
                    raise ValueError(
                        'Plot keyword not in correct format (key or key_i)'
                        + f': {idx}'
                    )

    def _make_new_cat(self):
        """Make New Catalogue.

        This function interprets the contents for each column of the new
        catalogue.

        Raises
        ------
        ValueError
            For invalid output format

        """
        if len(self._new_cat) == 0:
            return None

        self.new_cat = {}
        for key in self._new_cat.keys():
            self.new_cat[key] = {}
            for idx in self._new_cat[key]:
                s = re.split('=', j)
                if len(s) == 2:
                    if s[0] == 'OUTPUT_FORMAT':
                        self.new_cat[key][s[0]] = s[1]
                    else:
                        self.new_cat[key][s[0]] = StrInterpreter(
                            s[1],
                            self._data,
                            make_compare=False,
                            mask_dict=self.mask,
                        ).result
                else:
                    raise ValueError(f'Not a valid format : {idx}')

    def _make_rand_split(self):
        """Make Random Split.

        This function creates mask with random indices corresponding to the
        specfied ratio.

        Raises
        ------
        ValueError
            For invalid random split format
        ValueError
            If ratio is not a number
        ValueError
            If mask does not exist

        """
        if len(self._rand_split) == 0:
            return None

        self.rand_split = {}
        mask = np.ones(self._cat_size, dtype=bool)
        for key in self._rand_split.keys():
            self.rand_split[key] = {}
            for idx in self._rand_split[key]:
                s = re.split('=', idx)
                if len(s) != 2:
                    raise ValueError(
                        f'Not a valid format : {self._rand_split[key][0]}'
                    )
                if s[0] == 'RATIO':
                    try:
                        ratio = float(s[1])
                    except Exception:
                        raise ValueError('RATIO is not a number')
                    if ratio >= 1:
                        ratio /= 100.
                elif s[0] == 'MASK':
                    ss = re.split(',', s[1])
                    for k in ss:
                        try:
                            mask &= self.mask[k]
                        except Exception:
                            raise ValueError(f'mask {k} does not exist')

            cat_size = len(np.where(mask)[0])
            n_keep = int(np.ceil(cat_size * ratio))
            mask_ratio = []
            mask_left = list(range(0, cat_size))
            while len(mask_ratio) != n_keep:
                idx = np.random.randint(0, len(mask_left))
                mask_ratio.append(mask_left.pop(idx))
            mask_ratio = np.array(mask_ratio)
            mask_left = np.array(mask_left)
            self.rand_split[key]['mask'] = mask
            self.rand_split[key][f'ratio_{int(ratio * 100)}'] = mask_ratio
            self.rand_split[key][f'ratio_{100 - int(ratio * 100)}'] = mask_left

    def _make_stat(self):
        """Make Statistics.

        This function interprets the different statistics required.

        Raises
        ------
        ValueError
            For invalid statistics format

        """
        if len(self._stat) == 0:
            return None

        self.stat = {}
        for key in self._stat.keys():
            self.stat[key] = {}
            for idx in self._stat[key]:
                s = re.split('=', idx)
                if len(s) != 2:
                    raise ValueError(f'Not a valid format : {idx}')
                self.stat[key][s[0]] = StrInterpreter(
                    s[1],
                    self._data,
                    make_compare=False,
                    mask_dict=self.mask,
                ).result


class SEPlot(object):
    """The SEPlot Class.

    Tools to create plots.

    Parameters
    ----------
    plot_dict : dict
        Dictionary containing the parameters for the plot
    catalogue : numpy.recarray or astropy.fits.fitsrec
        Array containing the full data
    output_path : str
        Path for the output
    mask_dict : dict, optional
        Dictionary containing masks to apply

    Raises
    ------
    ValueError
        If ``plot_dict`` not provided
    ValueError
        If ``catalogue`` not provided
    ValueError
        If ``output_path`` not provided
    ValueError
        If plot type not specified
    ValueError
        For invalid plot type

    Notes
    -----
    Types of plots available : plot, scatter, hist.

    """

    def __init__(self, plot_dict, catalogue, output_path, mask_dict=None):

        if plot_dict is None:
            raise ValueError('plot_dict not provided')
        if catalogue is None:
            raise ValueError('catalogue not provided')
        if output_path is None:
            raise ValueError('output_path not provided')

        self._plot = plot_dict
        self._output_path = output_path
        self._cat = catalogue
        self._mask_dict = mask_dict

        if 'TYPE' not in self._plot.keys():
            raise ValueError('Plot type not specified')

        if self._plot['TYPE']['0'] in ['plot', 'PLOT']:
            self._check_key_for_plot(['X', 'Y'])
            self._make_plot()
        elif self._plot['TYPE']['0'] in ['scatter', 'SCATTER']:
            self._check_key_for_plot(['X', 'Y'])
            self._make_scatter()
        elif self._plot['TYPE']['0'] in [
            'histogram',
            'hist',
            'HISTOGRAM',
            'HIST'
        ]:
            self._check_key_for_plot(['Y'])
            self._make_hist()
        else:
            raise ValueError(f'Type : {self._plot["TYPE"]["0"]} not available')

    def _check_key_for_plot(self, key_list):
        """Check Key for Plot.

        Raise exception if keys not found in plot description.

        Parameters
        ----------
        key_list: list
            List of keys

        Raises
        ------
        ValueError
            If key not provided for plot type

        """
        for key in key_list:
            if key not in self._plot.keys():
                raise ValueError(
                    f'Key \'{key}\' not provided for plot of type '
                    + f'\'{self._plot["TYPE"]["0"]}\''
                )

    def _make_plot(self):
        """Make Plot.

        This function calls ``matplotlib.pyplot.plot``.

        Raises
        ------
        ValueError
            If X not provided for each value of Y
        ValueError
            If plot keyword not in the correct format

        See Also
        --------
        matplotlib.pyplot.plot

        """
        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                ii = 1
                for idx in s[1:-1]:
                    if ii % 2 == 0:
                        title += idx
                    else:
                        title += str(StrInterpreter(
                            idx,
                            self._cat,
                            make_compare=False,
                            mask_dict=self._mask_dict,
                        ).result)
                    ii += 1
        else:
            title = ''

        self._fig.suptitle(title)

        for key in self._plot['Y'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][key]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(StrInterpreter(
                                    j,
                                    self._cat,
                                    make_compare=False,
                                    mask_dict=self._mask_dict,
                                ).result)
                            jj += 1
                except Exception:
                    label = None
            else:
                label = None
            if 'COLOR' in self._plot.keys():
                try:
                    color = self._plot['COLOR'][key]
                except Exception:
                    color = None
            else:
                color = None
            if 'MARKER' in self._plot.keys():
                try:
                    marker = self._plot['MARKER'][key]
                except Exception:
                    marker = '+'
            else:
                marker = '+'
            if 'MARKERSIZE' in self._plot.keys():
                try:
                    markersize = self._plot['MARKERSIZE'][key]
                except Exception:
                    markersize = 1
            else:
                markersize = 1
            if 'LINE' in self._plot.keys():
                try:
                    line = self._plot['LINE'][key]
                except Exception:
                    line = ''
            else:
                line = ''
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = self._plot['ALPHA'][key]
                except Exception:
                    alpha = None
            else:
                alpha = None

            try:
                x = self._plot['X'][key]
            except Exception:
                if len(self._plot['X']) == 1:
                    x = self._plot['X'][self._plot['X'].keys()[0]]
                else:
                    raise ValueError(
                        'You need to specify X for each Y provided if they '
                        + 'dont have the same'
                    )

            plt.plot(
                StrInterpreter(x, self._cat, mask_dict=self._mask_dict).result,
                StrInterpreter(
                    self._plot['Y'][key],
                    self._cat,
                    mask_dict=self._mask_dict,
                ).result,
                label=label,
                color=color,
                marker=marker,
                markersize=markersize,
                ls=line,
                alpha=alpha,
                figure=self._fig,
            )

        # Set ploy limits for x and y
        for (lim, set_lim) in zip(['XLIM', 'YLIM'], [plt.xlim, plt.ylim]):
            if lim in self._plot.keys():
                try:
                    val = re.split(',', self._plot[lim]['0'])
                except Exception:
                    raise ValueError(
                        f'Plot {lim} keyword/value not in correct format '
                        + f'({lim}=lower,upper): {self._plot[lim]["0"]}'
                    )

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

        self._fig.savefig(
            f'{self._output_path}.{out_format.lower()}',
            format=out_format,
        )
        plt.close()

    def _make_scatter(self):
        """Make Scatter.

        This function calls ``matplotlib.pyplot.scatter``.

        Raises
        ------
        ValueError
            If X not provided for every ``SCATTER``
        ValueError
            If Y not provided for every ``SCATTER``

        See Also
        --------
        matplotlib.pyplot.scatter

        """
        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                counter = 1
                for idx in s[1:-1]:
                    if counter % 2 == 0:
                        title += idx
                    else:
                        title += str(StrInterpreter(
                            idx,
                            self._cat,
                            make_compare=False,
                            mask_dict=self._mask_dict,
                        ).result)
                    counter += 1
        else:
            title = ''

        self._fig.suptitle(title)

        for key in self._plot['SCATTER'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][key]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(StrInterpreter(
                                    j,
                                    self._cat,
                                    make_compare=False,
                                    mask_dict=self._mask_dict,
                                ).result)
                            jj += 1
                except Exception:
                    label = None
            else:
                label = None
            if 'MARKER' in self._plot.keys():
                try:
                    marker = self._plot['MARKER'][key]
                except Exception:
                    marker = '+'
            else:
                marker = '+'
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = self._plot['ALPHA'][key]
                except Exception:
                    alpha = None
            else:
                alpha = None

            try:
                x = self._plot['X'][key]
            except Exception:
                if len(self._plot['X']) == 1:
                    x = self._plot['X'][self._plot['X'].keys()[0]]
                else:
                    raise ValueError(
                        'You need to specify X for each SCATTER provided if '
                        + 'they dont have the same'
                    )
            try:
                y = self._plot['Y'][key]
            except Exception:
                if len(self._plot['Y']) == 1:
                    y = self._plot['Y'][self._plot['Y'].keys()[0]]
                else:
                    raise ValueError(
                        'You need to specify Y for each SCATTER provided if '
                        + 'they dont have the same'
                    )

            plt.scatter(
                StrInterpreter(x, self._cat, mask_dict=self._mask_dict).result,
                StrInterpreter(y, self._cat, mask_dict=self._mask_dict).result,
                c=StrInterpreter(
                    self._plot['SCATTER'][key],
                    self._cat,
                    mask_dict=self._mask_dict,
                ).result,
                label=label,
                marker=marker,
                alpha=alpha,
                figure=self._fig,
            )

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
            out_format = 'PNG'

        self._fig.savefig(
            f'{self._output_path}.{out_format.lower()}',
            format=out_format,
        )
        plt.close()

    def _make_hist(self):
        """Make Hist.

        This function calls ``matplotlib.pyplot.hist``.

        See Also
        --------
        matplotlib.pyplot.hist

        """
        self._fig = plt.figure()

        if 'TITLE' in self._plot.keys():
            title = self._plot['TITLE']['0']
            s = re.split('@', title)
            if len(s) >= 3:
                title = s[0]
                counter = 1
                for idx in s[1:-1]:
                    if counter % 2 == 0:
                        title += idx
                    else:
                        title += str(StrInterpreter(
                            idx,
                            self._cat,
                            make_compare=False,
                            mask_dict=self._mask_dict,
                        ).result)
                    counter += 1
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

        for key in self._plot['Y'].keys():
            if 'LABEL' in self._plot.keys():
                try:
                    label = self._plot['LABEL'][key]
                    s = re.split('@', label)
                    if len(s) >= 3:
                        label = s[0]
                        jj = 1
                        for j in s[1:-1]:
                            if jj % 2 == 0:
                                label += j
                            else:
                                label += str(StrInterpreter(
                                    j,
                                    self._cat,
                                    make_compare=False,
                                    mask_dict=self._mask_dict,
                                ).result)
                            jj += 1
                except Exception:
                    label = None
            else:
                label = None
            if 'COLOR' in self._plot.keys():
                try:
                    color = self._plot['COLOR'][key]
                except Exception:
                    color = None
            else:
                color = None
            if 'BIN' in self._plot.keys():
                try:
                    bins = int(self._plot['BIN'][key])
                except Exception:
                    if len(self._plot['BIN']) == 1:
                        bins = int(
                            self._plot['BIN'][self._plot['BIN'].keys()[0]]
                        )
            else:
                bins = 50
            if 'ALPHA' in self._plot.keys():
                try:
                    alpha = float(self._plot['ALPHA'][key])
                except Exception:
                    alpha = None
            else:
                alpha = None

            plt.hist(
                StrInterpreter(
                    self._plot['Y'][key],
                    self._cat,
                    mask_dict=self._mask_dict,
                ).result,
                bins=bins,
                color=color,
                label=label,
                alpha=alpha,
                histtype=htype,
                log=log,
            )

        if 'LABEL' in self._plot.keys():
            plt.legend()
        if 'XLABEL' in self._plot.keys():
            plt.xlabel(self._plot['XLABEL']['0'])
        if 'YLABEL' in self._plot.keys():
            plt.ylabel(self._plot['YLABEL']['0'])
        if 'FORMAT' in self._plot.keys():
            out_format = self._plot['FORMAT']['0']
        else:
            out_format = 'PNG'

        self._fig.savefig(
            f'{self._output_path}.{out_format.lower()}',
            format=out_format,
        )
        plt.close()
