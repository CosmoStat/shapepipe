# -*- coding: utf-8 -*-

"""ARGUMENT HANDLING

This module defines methods for handling the pipeline arguments.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import argparse as ap
from inspect import getmembers, isfunction
from shapepipe.modules import module_runners
from shapepipe.info import shapepipe_logo, __version__


class cutomFormatter(ap.ArgumentDefaultsHelpFormatter,
                     ap.RawDescriptionHelpFormatter):
    pass


def print_message(additional_arg):

    class customAction(ap.Action):

        def __init__(self, option_strings, version=None, dest=ap.SUPPRESS,
                     default=ap.SUPPRESS, help=help):
            super(customAction, self).__init__(
                option_strings=option_strings,
                dest=dest,
                default=default,
                nargs=0,
                help=help)

        def __call__(self, parser, args, values, option_string=None):
            print(additional_arg)
            exit()

    return customAction


def create_arg_parser():
    """ Create Argument Parser

    This method returns an argument parser.

    Returns
    -------
    ArgumentParser

    """

    # Create parser
    parser = ap.ArgumentParser(add_help=False, description=shapepipe_logo(),
                               formatter_class=cutomFormatter)
    optional = parser.add_argument_group('Optional Arguments')

    # Add arguments
    optional.add_argument('-h', '--help', action='help',
                          help='show this help message and exit')

    optional.add_argument('-v', '--version', action='version',
                          version='%(prog)s v{}'.format(__version__))

    optional.add_argument('-l', '--list_modules',
                          action=print_message('ShapePipe modules currently '
                                               'available: '
                                               '{}'.format(get_module_list())),
                          help='list modules currently available and exit')

    optional.add_argument('-c', '--config', default='config.ini',
                          help='configuration file name')

    # Return parser
    return parser.parse_args()


def get_module_list():
    """ Get Module List

    This method returns a list of the modules current available in
    module_runners.

    Returns
    -------
    list of module names

    """

    return list(zip(*getmembers(module_runners, isfunction)))[0]
