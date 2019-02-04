# -*- coding: utf-8 -*-

"""FILE HANDLER

This module defines a class for handling pipeline files.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
import re
from glob import glob
from shapepipe.pipeline.run_log import RunLog
from shapepipe.modules.module_runners import get_module_runners


def find_files(path, pattern='*', ext='*', empty_error=False):
        """ Find Files

        This method recursively retrieves file names from a given path that
        match a given pattern and/or have a given extension.

        Parameters
        ----------
        path : str
            Full path to files
        pattern : str, optional
            File pattern, default is '*'
        ext : str, optional
            File extension, default is '*'
        empty_error : bool, optional
            Raise error if file list is empty, default is False

        Returns
        -------
        list
            List of file names

        Raises
        ------
        ValueError
            For '*' in pattern
        ValueError
            For '*' in extension
        ValueError
            For invalid extension format
        RuntimeError
            For empty file list

        """

        dot = '.'
        star = '*'

        if pattern != star and star in pattern:
            raise ValueError('Do not include "*" in pattern.')

        if ext != star and star in ext:
            raise ValueError('Do not include "*" in extension.')

        if (not ext.startswith(dot) and dot in ext) or (ext.count(dot) > 1):
            raise ValueError('Invalid extension format: "{}".'.format(ext))

        if ext != star and not ext.startswith(dot):
            ext = dot + ext

        search_string = '{}/**/*{}*{}'.format(path, pattern, ext)

        file_list = sorted(glob(search_string, recursive=True))

        if empty_error and not file_list:
            raise RuntimeError('No files found matching the conditions in {}'
                               '.'.format(path))

        return [file for file in file_list if not os.path.isdir(file)]


class FileHandler(object):
    """ File Handler

    This class manages the files used and produced during a pipeline run.

    Parameters
    ----------
    run_name : str
        Run name
    module_list : list
        List of modules to be run
    config : CustomParser
        Configuaration parser instance

    """

    def __init__(self, run_name, modules, config):

        self._run_name = run_name
        self._module_list = modules
        self._config = config
        self.module_runners = get_module_runners(self._module_list)
        self._input_list = config.getlist('FILE', 'INPUT_DIR')
        self._output_dir = config.getexpanded('FILE', 'OUTPUT_DIR')
        self._log_name = config.get('FILE', 'LOG_NAME')
        self._run_log_file = self.format(self._output_dir,
                                         config.get('FILE', 'RUN_LOG_NAME'),
                                         '.txt')
        self._module_dict = {}

        if config.has_option('FILE', 'FILE_PATTERN'):
            self._file_pattern = config.getlist('FILE', 'FILE_PATTERN')
        else:
            self._file_pattern = None
        if config.has_option('FILE', 'FILE_EXT'):
            self._file_ext = config.getlist('FILE', 'FILE_EXT')
        else:
            self._file_ext = None
        if config.has_option('FILE', 'NUMBERING_SCHEME'):
            self._numbering_scheme = config.get('FILE', 'NUMBERING_SCHEME')
        else:
            self._numbering_scheme = r'RE:\_\d+'

    @property
    def run_dir(self):
        """ Run Directory

        This method defines the run directory.

        """

        return self._run_dir

    @run_dir.setter
    def run_dir(self, value):

        self.check_dir(value)

        self._run_dir = value

    @property
    def input_dir(self):
        """ Input Directory

        This method defines the input directory.

        """

        return self._input_dir

    @input_dir.setter
    def input_dir(self, value):

        self.check_dir(value)

        self._input_dir = value

    @staticmethod
    def check_dir(dir_name):
        """ Check Directory

        Raise error if directory exists.

        Parameters
        ----------
        dir_name : str
            Directory name

        Raises
        ------
        OSError
            If directory already exists

        """

        if os.path.isdir(dir_name):
            raise OSError('Directory {} already exists.'.format(dir_name))

    @classmethod
    def mkdir(cls, dir_name):
        """ Make Directory

        This method creates a directory at the specified path.

        Parameters
        ----------
        dir_name : str
            Directory name with full path

        """

        cls.check_dir(dir_name)
        os.mkdir(dir_name)

    @staticmethod
    def format(path, name, ext=''):
        """ Format Path Name

        This method appends the file/directory name to the input path.

        Parameters
        ----------
        path : str
            Full path
        name : str
            File or directory name
        ext : str, optional
            File extension, default is ''

        Returns
        -------
        str
            Formated path

        """

        return '{}/{}{}'.format(path, name, ext)

    def _get_input_dir(self):
        """ Get Input Directory

        This method sets the module input directory

        """

        input_dir = []

        for dir in self._input_list:

            if os.path.isdir(dir):
                input_dir.append(dir)

            elif 'last' in dir.lower():
                module = dir.lower().split(':')[1]
                input_dir.append(self.format(self.format(
                                 self._run_log.get_last(module),
                                 module), 'output'))

            elif ':' in dir.lower():
                string, module = dir.lower().split(':')
                input_dir.append(self.format(self.format(
                                 self._run_log.get_run(string), module),
                                 'output'))

            else:
                raise ValueError('Invalid INPUT_DIR. Make sure the paths '
                                 'provided are valid directories or use the '
                                 'allowed special keys.')

        self._input_dir = input_dir

    def create_global_run_dirs(self):
        """ Create Global Run Directories

        This method creates the pipeline output directories for a given run.

        """

        self.run_dir = self.format(self._output_dir, self._run_name)
        self._log_dir = self.format(self.run_dir, 'logs')
        self.log_name = self.format(self._log_dir, self._log_name)
        self._run_log = RunLog(self._run_log_file, self._module_list,
                               self.run_dir)

        self.mkdir(self.run_dir)
        self.mkdir(self._log_dir)

        self._get_input_dir()

    def get_add_module_property(self, module, property):
        """ Get Additional Module Properties

        Get a list of additional module property values.

        Parameters
        ----------
        module : str
            Module name
        property : str
            Property name

        Returns
        -------
        list
            Additional module property values

        """

        if (self._config.has_option(module.upper(), 'ADD_{}'.format(
                property.upper()))):

            return self._config.getlist(module.upper(), 'ADD_{}'.format(
                                        property.upper()))

    def _set_module_property(self, module, property):
        """ Get Module Property

        Get a module property from either the configuration file or the module
        runner.

        Parameters
        ----------
        module : str
            Module name
        property : str
            Property name

        """

        if self._config.has_option(module.upper(), property.upper()):
            prop_val = self._config.get(module.upper(), property.upper())

        else:
            prop_val = getattr(self.module_runners[module], property)

        if len(self._module_dict) == 1 or isinstance(prop_val, type(None)):
            prop_val = getattr(self, '_{}'.format(property))

        print('MKDEBUG set property {} to {}'.format(property, prop_val))

        self._module_dict[module][property] = prop_val

    def _set_module_list_property(self, module, property):
        """ Get Module List Property

        Get a module list property from either the configuration file or the
        module runner.

        Parameters
        ----------
        module : str
            Module name
        property : str
            Property name

        """

        if self._config.has_option(module.upper(), property.upper()):
            prop_list = self._config.getlist(module.upper(), property.upper())

        elif (property in ('file_pattern', 'file_ext') and not
                isinstance(getattr(self, '_{}'.format(property)), type(None))
                and len(self._module_dict) == 1):
            prop_list = getattr(self, '_{}'.format(property))

        else:
            prop_list = getattr(self.module_runners[module], property)

        if self.get_add_module_property(module, property):
            prop_list += self.get_add_module_property(module, property)

        self._module_dict[module][property] = prop_list

    def _set_module_properties(self, module):
        """ Get Module Properties

        Get module properties defined in module runner wrapper.

        Parameters
        ----------
        module : str
            Module name

        """

        module_props = ('numbering_scheme',)
        module_list_props = ('input_module', 'file_pattern', 'file_ext',
                             'depends', 'executes')

        [self._set_module_property(module, property) for property in
         module_props]
        [self._set_module_list_property(module, property) for property in
         module_list_props]

        # Make sure the number of patterns and extensions match
        if ((len(self._module_dict[module]['file_ext']) == 1) and
                (len(self._module_dict[module]['file_pattern']) > 1)):
            self._module_dict[module]['file_ext'] = \
                [self._module_dict[module]['file_ext'][0] for i in
                 self._module_dict[module]['file_pattern']]

        if (len(self._module_dict[module]['file_ext']) !=
                len(self._module_dict[module]['file_pattern'])):
            raise ValueError('The number of file_ext values does not match '
                             'the number of file_pattern values.')

    def _create_module_run_dirs(self, module):
        """ Create Module Run Directories

        This method creates the module output directories for a given run.

        Parameters
        ----------
        module : str
            Module name

        """

        self._module_dict[module]['run_dir'] = \
            (self.format(self._run_dir, module))
        self._module_dict[module]['log_dir'] = \
            (self.format(self._module_dict[module]['run_dir'], 'logs'))
        self._module_dict[module]['output_dir'] = \
            (self.format(self._module_dict[module]['run_dir'], 'output'))

        self.mkdir(self._module_dict[module]['run_dir'])
        self.mkdir(self._module_dict[module]['log_dir'])
        self.mkdir(self._module_dict[module]['output_dir'])

        # Set current output directory to module output directory
        self.output_dir = self._module_dict[module]['output_dir']

    def _set_module_input_dir(self, module):
        """ Set Module Input Directory

        Set the module input directory. If the module specified is the
        first module in the pipeline or does not have any input modules then
        only the `INPUT_DIR` from `[FILE]` is used, otherwise the output
        directory from the preceding module is used.

        Additional input directories can be specified with `INPUT_DIR` from
        `[MODULE]`.

        Parameters
        ----------
        module : str
            Module name

        """

        if (isinstance(self._module_dict[module]['input_module'], type(None))
                or len(self._module_dict) == 1):
            self._module_dict[module]['input_dir'] = self._input_dir

        else:
            self._module_dict[module]['input_dir'] = \
                ([self._module_dict[input_mod]['output_dir'] for input_mod in
                  self._module_dict[module]['input_module']])

            if self._config.has_option(module.upper(), 'INPUT_DIR'):
                self._module_dict[module]['input_dir'] = \
                    self._config.getlist(module.upper(), 'INPUT_DIR')

        if self.get_add_module_property(module, 'input_dir'):
            self._module_dict[module]['input_dir'] += \
                self.get_add_module_property(module, 'input_dir')

    @classmethod
    def _all_pattern_per_dir(cls, dir_list, pattern_list, ext_list):
        """ All Patterns Per Directory

        Find module files assuming all file patterns are applied to all input
        directories.

        Parameters
        ----------
        dir_list : list
            List of input directories
        pattern_list : list
            List of file patterns
        ext_list : list
            List of file extensions

        Returns
        -------
        list
            List of input files

        """

        file_list = []

        for pattern, ext in zip(pattern_list, ext_list):

            sub_file_list = [find_files(dir, pattern, ext) for dir in dir_list]
            sub_file_list = cls._flatten_list(sub_file_list)

            if not sub_file_list:
                raise RuntimeError('No files found matching patterns ({})'
                                   ' and or extensions ({}) in directories '
                                   ' provided ({})'.format(pattern, ext,
                                                           dir_list))
            else:
                file_list.append(sub_file_list)

        return file_list

    @staticmethod
    def _generate_re_pattern(match_pattern):
        """ Generate Regular Expression Pattern

        Generate a regular expression pattern from an input string.

        Parameters
        ----------
        match_pattern : str
            Pattern string

        Returns
        -------
        _sre.SRE_Pattern
            Regular expression pattern

        Raises
        ------
        TypeError
            For invalid input type

        """

        if not isinstance(match_pattern, str):
            TypeError('Match pattern must be a string.')

        chars = [char for char in match_pattern if not char.isalnum()]
        split_pattern = '|'.join(chars).replace('.', r'\.')
        chars = ['\\{}'.format(char) for char in chars] + ['']
        num_length = ['\\d{{{}}}'.format(len(digits)) for digits in
                      re.split(split_pattern, match_pattern)]
        re_pattern = r''.join([a for b in zip(num_length, chars)
                               for a in b]).replace('{1}', '+')

        return re.compile(re_pattern)

    @staticmethod
    def _strip_dir_from_file(file_name, dir_list):
        """ Strip Directory from File Name

        Remove the directory string from the file name.

        Parameters
        ----------
        file_name : str
            File name
        dir_list : list
            Input directory list

        Returns
        -------
        str
            File name

        """

        return [file_name.replace(_dir + '/', '') for _dir in dir_list
                if _dir in file_name][0]

    @classmethod
    def _get_file_pattern(cls, file_name, dir_list, re_pattern):
        """ Get File Pattern

        Get the string component of the input file name matching the specified
        regular expression.

        Parameters
        ----------
        file_name : str
            File name
        dir_list : list
            Input directory list
        re_pattern : str
            Regular rexpression pattern

        Returns
        -------
        str
            File pattern

        """

        file_name = cls._strip_dir_from_file(file_name, dir_list)

        if not re.search(re_pattern, file_name):
            raise RuntimeError('No files found matching {}'.format(re_pattern))

        return re.search(re_pattern, file_name).group()

    @staticmethod
    def _flatten_list(input_list):
        """ Flatten List

        Flatten a list of lists.

        Parameters
        ----------
        input_list : list
            A list of lists

        Returns
        -------
        list
            Flattened list

        """

        return [item for sublist in input_list for item in sublist]

    @classmethod
    def _check_pattern(cls, file_list, dir_list, pattern):
        """ Check Pattern

        Find files in file list that match the input pattern.

        Parameters
        ----------
        file_list : list
            List of file names
        pattern : str
            File pattern

        Returns
        -------
        list
            List of files matching the pattern

        """

        new_list = []

        for file in cls._flatten_list(file_list):

            item = cls._strip_dir_from_file(file, dir_list)

            if (len(item.split(pattern)) > 1 and
                (not item.split(pattern)[1] or not
                 item.split(pattern)[1][0].isdigit())):
                new_list.append(file)

        return new_list

    @classmethod
    def _match_list_items(cls, file_list, dir_list, re_pattern):
        """ Match List Items

        Match files names in a list of lists.

        Parameters
        ----------
        file_list : list
            List of file names
        dir_list : list
            Input directory list
        re_pattern : str
            Regular rexpression pattern

        Returns
        -------
        tuple
            List of matched file names, list of unmatched file names

        Raises
        ------
        TypeError
            For invalid input type

        """

        if not isinstance(file_list, list):
            TypeError('File list must be a list.')

        all_patterns = [cls._get_file_pattern(file, dir_list, re_pattern)
                        for file in max(file_list, key=len)]
        new_list = [cls._check_pattern(file_list, dir_list, pattern)
                    for pattern in all_patterns]
        new_list = [item for item in new_list if item]

        file_dict = dict(zip(all_patterns, new_list))

        max_n_cols = max([len(flist) for flist in new_list])
        matched = dict([(key, val) for key, val in file_dict.items() if
                        len(val) == max_n_cols])
        missing = dict([(key, val) for key, val in file_dict.items() if
                        len(val) < max_n_cols])

        return matched, missing

    @classmethod
    def _check_file_list(cls, file_list, dir_list, num_scheme):
        """ Check File List

        Check the file list for missing files.

        Parameters
        ----------
        file_list : list
            List of file names
        dir_list : list
            Input directory list
        num_scheme : str
            Numbering scheme

        Returns
        -------
        tuple
            List of matched file names, list of unmatched file names

        """

        if num_scheme.startswith('RE:'):

            re_pattern = num_scheme.replace('RE:', '')

        else:

            re_pattern = cls._generate_re_pattern(num_scheme)

        return cls._match_list_items(file_list, dir_list, re_pattern)

    def _get_module_input_files(self, module):
        """ Get Module Input Files

        Retrieve the module input files names from the input directory.

        Parameters
        ----------
        module : str
            Module name

        """

        dir_list = self._module_dict[module]['input_dir']
        pattern_list = self._module_dict[module]['file_pattern']
        ext_list = self._module_dict[module]['file_ext']
        num_scheme = self._module_dict[module]['numbering_scheme']

        file_list = self._all_pattern_per_dir(dir_list, pattern_list,
                                              ext_list)

        self.process_list, self.missed = self._check_file_list(file_list,
                                                               dir_list,
                                                               num_scheme)

    def set_up_module(self, module):
        """ Set Up Module

        Set up module parameters for file handler.

        Parameters
        ----------
        module : str
            Module name

        """

        self._module_dict[module] = {}
        self._set_module_properties(module)
        self._create_module_run_dirs(module)
        self._set_module_input_dir(module)
        self._get_module_input_files(module)

    def get_worker_log_name(self, module, job_name, file_number_string):
        """ Get Worker Log Name

        This method generates a worker log name.

        Parameters
        ----------
        module : str
            Module name
        job_name : str
            Job name
        file_number_string : str
            File numbering in output

        Returns
        -------
        str
            Worker log file name

        """

        return '{}/{}_file{}'.format(self._module_dict[module]['log_dir'],
                                     job_name, file_number_string)
