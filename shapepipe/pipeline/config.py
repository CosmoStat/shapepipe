# -*- coding: utf-8 -*-

"""CONFIGURATION FILE HANDLING

This module defines methods for handling the pipeline configuration file.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
from configparser import ConfigParser


class CustomParser(ConfigParser):
    """ Custom Parser

    This class adds functionality to the ``ConfigParser`` class.

    """

    def getexpanded(self, section, option):
        """ Get Expanded

        This method expands enviroment varibles obtaiened using the get method.

        Parameters
        ----------
        section : str
            Configuration file section
        option : str
            Configuration file option

        Returns
        -------
        str
            Expanded enviroment variables

        """

        return self._get(section, os.path.expandvars, option)

    def getlist(self, section, option, delimiter=','):
        """ Get List

        This method retrieves a list of strings separated by a given
        delimiter.

        Parameters
        ----------
        section : str
            Configuration file section
        option : str
            Configuration file option
        delimiter : str, optional
            Delimiter between list entries, default is ','

        Returns
        -------
        list
            List of strings

        """

        return [opt.strip() for opt in self.getexpanded(section,
                option).split(delimiter)]


class SetUpParser(object):
    """ Set Up Parser

    This class sets up an instance of ``CustomParser`` and checks the
    pipeline related parameters.

    Parameters
    ----------
    file_name : str
        Configuration file name

    """

    def __init__(self, file_name):

        self.file_name = file_name
        self.config = CustomParser()

    @property
    def file_name(self):
        """ File name

        This sets the configuration file name.

        Raises
        ------
        IOError
            For non existent configuration file

        """

        return self._file_name

    @file_name.setter
    def file_name(self, value):

        if not os.path.exists(value):
            raise IOError('Configuration file {} does not exist.'.format(
                          value))

        self._file_name = value

    def _set_defaults(self):
        """ Set Defaults

        Set default configuration options.

        """

        if not self.config.has_option('DEFAULT', 'RUN_NAME'):
            self.config.set('DEFAULT', 'RUN_NAME', 'shapepipe_run')

        if not self.config.has_option('DEFAULT', 'RUN_DATETIME'):
            self.config.set('DEFAULT', 'RUN_DATETIME', 'True')

        if not self.config.has_option('DEFAULT', 'VERBOSE'):
            self.config.set('DEFAULT', 'VERBOSE', 'True')

    def _set_execution_options(self):
        """ Set Execution Options

        This method checks the execution options in the configuration file.

        Raises
        ------
        RuntimeError
            For no module runner specified
        ImportError
            For non-existent module runner

        """

        if not self.config.has_option('EXECUTION', 'MODULE'):
            raise RuntimeError('No module(s) specified')

        if not self.config.has_option('EXECUTION', 'MODE'):
            self.config.set('EXECUTION', 'MODE', 'smp')

    def _set_file_options(self):
        """ Set File Options

        This module checks the file options in the configuration file.

        Raises
        ------
        RuntimeError
            For no input directory specified
        OSError
            For non-existent input directory
        RuntimeError
            For no output directory specified
        OSError
            For non-existent output directory

        """

        if not self.config.has_option('FILE', 'LOG_NAME'):
            self.config.set('FILE', 'LOG_NAME', 'shapepipe')

        if not self.config.has_option('FILE', 'RUN_LOG_NAME'):
            self.config.set('FILE', 'RUN_LOG_NAME', 'shapepipe_runs')

        if not self.config.has_option('FILE', 'INPUT_DIR'):
            raise RuntimeError('Not input directory specified')

        if not self.config.has_option('FILE', 'OUTPUT_DIR'):
            raise RuntimeError('Not output directory specified')

        elif not os.path.isdir(self.config.getexpanded('FILE', 'OUTPUT_DIR')):
            raise OSError('Directory {} not found.'.format(
                          self.config.getexpanded('FILE', 'OUTPUT_DIR')))

        if not self.config.has_option('FILE', 'CORRECT_FILE_PATTERN'):
            self.config.set('FILE', 'CORRECT_FILE_PATTERN', 'True')

    def _set_worker_options(self):
        """ Set Worker Options

        This module checks the worker options in the configuration file.

        """

        if not self.config.has_section('WORKER'):
            self.config.add_section('WORKER')

        if not self.config.has_option('WORKER', 'PROCESS_PRINT_LIMIT'):
            self.config.set('WORKER', 'PROCESS_PRINT_LIMIT', '200')

    def get_parser(self):
        """ Get Parser

        Return a configuration file parser instance.

        Returns
        -------
        CustomParser
            Custom configuration file parser

        """

        self.config.read(self.file_name)
        self._set_defaults()
        self._set_execution_options()
        self._set_file_options()
        self._set_worker_options()

        return self.config


def create_config_parser(file_name):
    """ Create Configuration Parser

    This method creates a configuration file parser instance.

    Parameters
    ----------
    file_name : str
        Configuration file name

    Returns
    -------
    CustomParser
        Custom configuration file parser

    """

    parser = SetUpParser(file_name).get_parser()
    parser.file_name = file_name

    return parser
