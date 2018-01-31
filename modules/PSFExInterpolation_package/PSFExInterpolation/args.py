# -*- coding: utf-8 -*-

"""ARGUMENTS MODULE

This module contain a class for defining the package command line arguments.

:Authors: Samuel Farrens and Marc Gentile

:Date: 06/11/2017

"""

# -- External Modules
from mpfg.mp_args import Args


class PackageArgs(Args):

    """Package Arguments

    This class defines the package command line arguments.

    Notes
    -----
    All of the class members simply override the defaults from Args in
    mp_args.py.

    """

    def __init__(self):

        Args.__init__(self)

    @property
    def default_config_dir(self):
        """Default Configuration File Directory

        Set the default configuration file directory.

        Notes
        -----
        This property overrides the property of the same name from Args in
        mp_args.py.

        """

        return './config'

    @property
    def default_config_filename(self):
        """Default Configuration File Name

        Set the default configuration file name.

        Notes
        -----
        This property overrides the property of the same name from Args in
        mp_args.py.

        """

        return 'config.cfg'

    def print_usage(self, prog_name):
        """Print Usage

        This method prints the usage message.

        Notes
        -----
        This method overrides the method of the same name from Args in
        mp_args.py.

        """

        # Print usage information
        print('\nUsage: {0} [options]'.format(prog_name))
        print('\nHelp:')
        print('-h,  --help\tprint this help')

        # Optinal execution arguments
        print('\nOptions:')
        print('-c,  --config-file\tconfiguration file name '
              '(default: {})'.format(self.default_config_filename))
        print('-d,  --config-dir\tconfiguration directory '
              '(default: {})'.format(self.default_config_dir))
        print('\n')
