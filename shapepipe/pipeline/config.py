# -*- coding: utf-8 -*-

"""CONFIGURATION FILE HANDLING

This module defines methods for handling the pipeline configuration file.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os
from configparser import ConfigParser
from ..modules import module_runners


class CustomParser(ConfigParser):
    """ Custom Parser

    This class adds functionality to the ConfigParser class.

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
        str with expanded enviroment variables

        """

        return self._get(section, os.path.expandvars, option)

    def getlist(self, section, option, delimiter=','):

        return [opt for opt in self.get(section, option).split(delimiter)]


def create_config_parser(file_name):
    """ Create Configuration Parser

    This method creates a configuration file parser instance.

    Parameters
    ----------
    file_name : str
        Configuration file name

    Returns
    -------
    ExpandingParser instance

    Raises
    ------
    IOError
        For non existent configuration file
    RuntimeError
        For no module runner specified
    RuntimeError
        For non-existent module runner

    """

    if not os.path.exists(file_name):
        raise IOError('Configuration file {} does not exist.'.format(
                      file_name))

    config = CustomParser()
    config.read(file_name)

    # Set default options
    if not config.has_option('DEFAULT', 'RUN_NAME'):
        config.set('DEFAULT', 'RUN_NAME', 'shapepipe_run')

    if not config.has_option('DEFAULT', 'RUN_DATETIME'):
        config.set('DEFAULT', 'RUN_DATETIME', 'True')

    if not config.has_option('DEFAULT', 'VERBOSE'):
        config.set('DEFAULT', 'VERBOSE', 'True')

    # Check execution options
    if not config.has_option('EXECUTION', 'MODULE'):
        raise RuntimeError('Not module(s) specified')

    for module in config.getlist('EXECUTION', 'MODULE'):
        if not hasattr(module_runners, module):
            raise ImportError('Module runner {} not found in {}.'.format(
                              module, module_runners.__name__))

    # Check file options
    if not config.has_option('FILE', 'LOG_NAME'):
        config.set('FILE', 'LOG_NAME', 'shapepipe')

    if not config.has_option('FILE', 'INPUT_DIR'):
        raise RuntimeError('Not input directory specified')

    elif not os.path.isdir(config.getexpanded('FILE', 'INPUT_DIR')):
        raise OSError('Directory {} not found.'.format(
                      config.getexpanded('FILE', 'INPUT_DIR')))

    if not config.has_option('FILE', 'OUTPUT_DIR'):
        raise RuntimeError('Not output directory specified')

    elif not os.path.isdir(config.getexpanded('FILE', 'OUTPUT_DIR')):
        raise OSError('Directory {} not found.'.format(
                      config.getexpanded('FILE', 'OUTPUT_DIR')))

    return config
