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


def find_files(path, pattern='*', ext='*'):
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

        if not file_list:
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
            self._num_scheme = config.get('FILE', 'NUMBERING_SCHEME')
        else:
            self._num_scheme = r'RE:\_\d+'

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

    def _get_module_properties(self, module):
        """ Get Module Properties

        Get module properties defined in module runner wrapper.

        Parameters
        ----------
        module : str
            Module name

        """

        # Get the name of the input module from module runner
        self._module_dict[module]['input_module'] = \
            self.module_runners[module].input_module

        # Get the input file pattern from module runner (or config file)
        if (not isinstance(self._file_pattern, type(None))
                and len(self._module_dict) == 1):
            self._module_dict[module]['file_pattern'] = self._file_pattern
        else:
            self._module_dict[module]['file_pattern'] = \
                self.module_runners[module].file_pattern

        # Get the input file extesion from module runner (or config file)
        if (not isinstance(self._file_ext, type(None))
                and len(self._module_dict) == 1):
            self._module_dict[module]['file_ext'] = self._file_ext
        else:
            self._module_dict[module]['file_ext'] = \
                self.module_runners[module].file_ext

        # Make sure the number of patterns and extensions match
        if ((len(self._module_dict[module]['file_ext']) == 1) and
                (len(self._module_dict[module]['file_pattern']) > 1)):
            self._module_dict[module]['file_ext'] = \
                [self._module_dict[module]['file_ext'][0] for i in
                 self._module_dict[module]['file_pattern']]

        elif ((len(self._module_dict[module]['file_pattern']) == 1) and
                (len(self._module_dict[module]['file_ext']) > 1)):
            self._module_dict[module]['file_pattern'] = \
                [self._module_dict[module]['file_pattern'][0] for i in
                 self._module_dict[module]['file_ext']]

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
                self._module_dict[module]['input_dir'] += \
                    self._config.getlist(module.upper(), 'INPUT_DIR')

    @staticmethod
    def _one_pattern_per_dir(dir_list, pattern_list, ext_list):
        """ One Pattern Per Directory

        Find module files assuming one file pattern per input directory.

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

        return [find_files(dir, pattern, ext) for dir, pattern, ext in
                zip(dir_list, pattern_list, ext_list)]

    @staticmethod
    def _all_pattern_per_dir(dir_list, pattern_list, ext_list):
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

        return [find_files(dir, pattern, ext) for pattern, ext in
                zip(pattern_list, ext_list) for dir in dir_list]

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

    def _strip_dir_from_file(self, file_name, dir_list):
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

        return [file_name.replace(_dir, '') for _dir in dir_list
                if _dir in file_name][0]

    def _get_file_pattern(self, file_name, dir_list):
        """ Get File Pattern

        Get the string component of the input file name matching the specified
        regular expression.

        Parameters
        ----------
        file_name : str
            File name
        dir_list : list
            Input directory list

        Returns
        -------
        str
            File pattern

        """

        file_name = self._strip_dir_from_file(file_name, dir_list)

        return re.search(self._re_pattern, file_name).group()

    def _match_list_items(self, file_list, dir_list):
        """ Match List Items

        Match files names in a list of lists.

        Parameters
        ----------
        file_list : list
            List of file names
        dir_list : list
            Input directory list

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

        all_patterns = [self._get_file_pattern(file, dir_list) for file in
                        max(file_list, key=len)]
        new_list = [[item for sublist in file_list for item in sublist if
                     pattern in self._strip_dir_from_file(item, dir_list)]
                    for pattern in all_patterns]

        file_dict = dict(zip(all_patterns, new_list))

        max_n_cols = max([len(flist) for flist in new_list])
        matched = dict([(key, val) for key, val in file_dict.items() if
                        len(val) == max_n_cols])
        missing = dict([(key, val) for key, val in file_dict.items() if
                        len(val) < max_n_cols])

        return matched, missing

    def _check_file_list(self, file_list, dir_list):
        """ Check File List

        Check the file list for missing files.

        Parameters
        ----------
        file_list : list
            List of file names
        dir_list : list
            Input directory list

        Returns
        -------
        tuple
            List of matched file names, list of unmatched file names

        """

        if self._num_scheme.startswith('RE:'):

            self._re_pattern = self._num_scheme.replace('RE:', '')

        else:

            self._re_pattern = self._generate_re_pattern(self._num_scheme)

        return self._match_list_items(file_list, dir_list)

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

        if len(dir_list) == len(pattern_list):
            file_list = self._one_pattern_per_dir(dir_list, pattern_list,
                                                  ext_list)
        else:
            file_list = self._all_pattern_per_dir(dir_list, pattern_list,
                                                  ext_list)

        self.process_list, self.missed = self._check_file_list(file_list,
                                                               dir_list)

    def set_up_module(self, module):
        """ Set Up Module

        Set up module parameters for file handler.

        Parameters
        ----------
        module : str
            Module name

        """

        self._module_dict[module] = {}
        self._get_module_properties(module)
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
