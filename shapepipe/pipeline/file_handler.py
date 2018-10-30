# -*- coding: utf-8 -*-

"""FILE HANDLER

This module defines a class for handling pipeline files.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
from shapepipe.pipeline.run_log import RunLog
from shapepipe.modules import module_runners


class FileHandler(object):
    """ File Handler

    This class manages the files used and produced during a pipeline run.

    Parameters
    ----------
    run_name : str
        Run name
    config : CustomParser
        Configuaration parser instance

    """

    def __init__(self, run_name, config):

        self._run_name = run_name
        self._input_str = config.getexpanded('FILE', 'INPUT_DIR')
        self._output_dir = config.getexpanded('FILE', 'OUTPUT_DIR')
        self._log_name = config.get('FILE', 'LOG_NAME')
        self._run_log_file = self.format(self._output_dir,
                                         config.get('FILE', 'RUN_LOG_NAME'),
                                         '.txt')
        self._module_dict = {}

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

        if os.path.isdir(self._input_str):
            input_dir = self._input_str

        elif 'last' in self._input_str.lower():
            module = self._input_str.lower().split(',')[1]
            input_dir = self.format(self.format(self._run_log.get_last(),
                                    module), 'output')

        else:
            string, module = self._input_str.lower().split(',')
            input_dir = self.format(self.format(self._run_log.get_run(string),
                                    module), 'output')

        self._input_dir = input_dir

    def create_global_run_dirs(self):
        """ Create Global Run Directories

        This method creates the pipeline output directories for a given run.

        """

        self.run_dir = self.format(self._output_dir, self._run_name)
        self._log_dir = self.format(self.run_dir, 'logs')
        self.log_name = self.format(self._log_dir, self._log_name)
        self._run_log = RunLog(self._run_log_file, self.run_dir)

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

        self._module_dict[module]['n_inputs'] = \
            getattr(module_runners, module).n_inputs
        self._module_dict[module]['input_module'] = \
            getattr(module_runners, module).input_module
        self._module_dict[module]['ext'] = getattr(module_runners, module).ext

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

        Specify the module input directory.

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
                (self._module_dict[self._module_dict[module]
                 ['input_module']]['output_dir'])

    def _get_module_input_files(self, module):
        """ Get Module Input Files

        Retrieve the module input files names from the input directory.

        Parameters
        ----------
        module : str
            Module name

        """

        self._module_dict[module]['files'] = \
            self._get_files_by_ext(self._module_dict[module]['input_dir'],
                                   self._module_dict[module]['ext'])

        self.process_list = self._module_dict[module]['files']

    def _get_files_by_ext(self, path, ext):
        """ Get Files by Extension

        This method retrieves file names from a given path with a given file
        extension.

        Parameters
        ----------
        path : str
            Full path to files
        ext : str
            File extension

        Returns
        -------
        list
            List of file names

        Raises
        ------
        RuntimeError
            For empty file list

        """

        file_list = [self.format(path, file) for file in os.listdir(path)
                     if file.endswith(ext)]

        if not file_list:
            raise RuntimeError('No files found matching the conditions in {}'
                               '.'.format(path))

        return file_list

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

    def get_worker_log_name(self, module, job_name):
        """ Get Worker Log Name

        This method generates a worker log name.

        Parameters
        ----------
        job_name : str
            Job name

        Returns
        -------
        str
            Worker log file name

        """

        return self.format(self._module_dict[module]['log_dir'], job_name)
