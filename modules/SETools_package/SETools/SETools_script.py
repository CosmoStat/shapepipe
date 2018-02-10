# -*- coding: utf-8 -*-
"""SETOOLS SCRIPT

This script contain a class to handle operation on SExtractor output catalog.

:Authors: Axel Guinot

:Date: 16/01/2017

"""

# Compability with python2.x for x>6
from __future__ import print_function


import numpy as np
import re
import operator

import os

import scatalog as sc


class SETools(object):
    """SETools class

    Tools to analyse SExtractor catalogs.

    Parameters
    ----------
    cat_filepath : str
        Path to SExtractor catalog (FITS_LDAC format)
    config_filepath : str
        Path to config.setools file
    output_dir : str
        Path to pipeline result directory
    stat_output_dir : str, optional, default=None
        Path to pipeline statistics output directory

    """

    def __init__(self, cat_filepath, config_filepath, output_dir, stat_output_dir=None):

        self._cat_filepath = cat_filepath
        self._config_filepath = config_filepath
        self._output_dir = output_dir
        if stat_output_dir:
            self._stat_output_dir = stat_output_dir

        self._cat_file = sc.FITSCatalog(self._cat_filepath, SEx_catalog=True)
        self._cat_file.open()
        self._cat_size = len(self._cat_file.get_data())

        self._config_file = open(self._config_filepath)

        self._init_stat_function()
        self._comp_dict = {'<': operator.lt,
                           '>': operator.gt,
                           '<=': operator.le,
                           '>=': operator.ge,
                           '==': operator.eq,
                           '!=': operator.ne}


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
                file_name = self._output_dir + '/mask/' + i + file_number + '.fits'
                self.save_mask(self.mask[i], file_name)

        ### Post-processing of output products
        # Statistics
        self._make_stat()
        for mask_type in self.stat.keys():
            file_name = self._stat_output_dir + '/' + mask_type + file_number + '.txt'
            self.save_stat(mask_type, file_name)
            


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
        self._hist={}
        self._stat={}
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
                elif re.split(':', sec)[0] =='HIST':
                    in_section = 3
                    try:
                        hist_name = re.split(':', sec)[1]
                    except:
                        hist_name = 'hist_{0}'.format(len(self._hist)+1)
                    self._hist[hist_name] = []
                elif re.split(':', sec)[0] == 'STAT':
                    in_section = 4
                    try:
                        stat_name = re.split(':', sec)[1]
                    except:
                        stat_name = 'stat_{0}'.format(len(self._stat)+1)
                    self._stat[stat_name] = []
                else:
                    raise Exception("Section has to be in ['MASK','PLOT','HIST','STAT']")
            else:
                if in_section == 1:
                    self._mask[mask_name].append(line_tmp)
                elif in_section == 2:
                    self._plot[plot_name].append(line_tmp)
                elif in_section == 3:
                    self._hist[hist_name].append(line_tmp)
                elif in_section == 4:
                    self._stat[stat_name].append(line_tmp)



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

        line_tmp = line.replace(' ','')
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
                mask_tmp &= self._compare(j)
            self.mask[i] = mask_tmp


    def _make_plot(self):
        pass

    def _make_hist(self):
        pass

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

            # Loop over statistics
            for stat in self._stat[mask_type]:

                if stat == 'NUMBER_OBJ':
                    # Number of selected objects
                    res = len(self._cat_file.get_data()[self.mask[mask_type]])
                    print('make_stat {} {}'.format(mask_type, res))
                    self.stat[mask_type][stat] = res
                else:
                    print('Warning: statistic \'{}\' not supported, continuing...'.format(stat))



    ##################
    # Stat functions #
    ##################

    def _apply_func(self, string):
        """Apply function

        Look for a function in a string and apply it.

        Parameters
        ----------
        string : str
            String containing the function.

        Returns
        -------
        float
            Result of the function.

        """

        s = re.split('\(|\)', string)

        if len(s) == 1:
            return self._param_value(s[0])
        elif len(s) == 3:
            try:
                return self._stat_func[s[0]](self._param_value(s[1]))
            except:
                raise Exception('Unknown function : {0}'.format(s[0]))
        else:
            raise Exception('Only one function can be applied.'
                            'Problem with the term : {0}'.format(string))


    def _init_stat_function(self):
        """Initialise available stat functions

        Create a dictionary containing the functions.

        """

        self._stat_func = {}

        self._stat_func['mean'] = np.mean
        self._stat_func['median'] = np.median
        self._stat_func['mode'] = self._mode


    def _mode(self, input, eps=0.001):
        """Compute Mode

        Compute the most frequent value of a continuous distribution.

        Parameters
        ----------
        input : numpy.ndarray
            Numpy array containing the data.
        eps : float, optional
            Accuracy to achieve (default is 0.001)

        Note
        ----

        The input array must have 10 or more elements.

        """

        if self._cat_size > 100:
            bins = int(float(self._cat_size)/10.)
        elif self._cat_size >= 10:
            bins = int(float(self._cat_size)/3.)
        else:
            raise ValueError("Can't compute with less than 10 elements in the input")

        data = input
        diff = eps+1.

        while diff>eps:
            hist = np.histogram(data, bins)

            b_min = hist[1][hist[0].argmax()]
            b_max = hist[1][hist[0].argmax()+1]

            diff = b_max - b_min

            data = data[(data > b_min) & (data < b_max)]

        return (b_min + b_max) / 2.


    #############################
    # String handling functions #
    #############################

    def _compare(self, string):
        """Handle comparison in a string

        This function transform condition in a string to real condition.

        Parameters
        ----------
        string : str
            strind containing the comparison.

        """

        comp='<|>|<=|>=|==|!='

        if len(re.split(comp,string))!=2:
            raise Exception('Only one comparison in [<, >, <=, >=, ==, !=] per lines')

        for i in ['<=', '>=', '<','>', '==', '!=']:
            terms = re.split(i, string)
            if len(terms) == 2:
                first = self._apply_func(terms[0])
                second = self._apply_func(terms[1])

                return self._comp_dict[i](first, second)


    def _param_value(self, param):
        """Param Value

        Return the value of a parameter from the header.

        Parameters
        ----------
        param : str
            Parameter of the catalog. (see Note)

        Returns
        -------
        float
            Value of the parameter. (see Note)

        Note
        ----
        You can give a number as parameter it will return this number as float.
        Also, you can enter a linear combinaition of parameters with numbers in it.

        Example
        -------
            param = 'GAIN+NEXP*TEXP+10'

        """

        if param is None:
            raise ValueError("Parameter not specified")

        param=re.sub(' ','',param)
        param_split=re.split('\*|\/|\-|\+',param)
        if len(param_split)==1:
            return self._get_value(param)
        else:
            return self._operate(param,param_split)


    def _get_value(self, param):
        """Get Value

        Return the value of the corresponding parameter. Or return a float with a number as parameter.

        Parameters
        ----------
        param : str
            Parameter of the catalog.

        Returns
        -------
        float
            Vvalue of the parameter. Or float

        Note
        ----
        You can't make operation here !

        """

        if param is None:
            raise ValueError("Parameter not specified")

        try:
            param_value=float(param)
            return param_value
        except:
            try:
                return self._cat_file.get_data()[param]
            except:
                raise ValueError('param has to be a float or a catalog parameter. {0} not found'.format(param))


    def _operate(self, param, param_split):
        """Handle operation in a string

        Make operation between catalog's parameters and/or numbers

        Parameters
        ----------
        param : str
            Parameter or linear combination of parameters.
        param_split : str
            The different parameter splitted using '\*|\/|\-|\+' as delimiter.

        Returns
        -------
        float
            Result of the operation

        Note
        ----
        It's used as a recursive function

        """

        op='\*|\/|\-|\+'
        if param is None:
            raise ValueError("Parameter not specified")
        if param_split is None:
            raise ValueError("Parameters splited not specified")

        if len(re.split(op,param))==1:
            return self._get_value(param)

        tmp = self._param_op_func(re.split('\+',param), param_split, operator.add, 0)
        if tmp != 'pass':
            return tmp
        else:
            tmp = self._param_op_func(re.split('\-',param), param_split, operator.sub, 'init')
            if tmp != 'pass':
                return tmp
            else:
                tmp = self._param_op_func(re.split('\*',param), param_split, operator.mul, 1)
                if tmp != 'pass':
                    return tmp
                else:
                    return self._param_op_func(re.split('\/',param), param_split, operator.div, 'init')


    def _param_op_func(self, param_op, param_split, op, tmp):
        """Make a specified operation

        This function handle the posible operation between parameters.

        Parameters
        ----------
        param_op : list
            List of parameters to operate.
        param_split : list
            The different parameter splitted using '\*|\/|\-|\+' as delimiter.
        op : func
            The kind of operation provide as an operator function
            (Example : operator.sub).
        tmp : str or float
            Temporary result of the global operation or value to
            initiate operation.

        Returns
        -------
        float or 'pass'
            Result of the operation or 'pass' if their is remaining operations.

        """

        if len(param_op) > 2:
            for i in param_op:
                if tmp == 'init':
                    tmp = self._operate(i, param_split)
                else:
                    tmp = op(tmp, self._operate(i, param_split))
            return tmp
        elif len(param_op) == 2:
            if param_op[0] in param_split:
                first = self._get_value(param_op[0])
            else:
                first = self._operate(param_op[0], param_split)
            if param_op[1] in param_split:
                second = self._get_value(param_op[1])
            else:
                second = self._operate(param_op[1], param_split)
            return op(first, second)
        else:
            return 'pass'
