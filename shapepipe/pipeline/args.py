"""ARGUMENT HANDLING.

This module defines methods for handling the pipeline arguments.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import argparse as ap

from shapepipe.info import __version__, shapepipe_logo
from shapepipe.modules import __module_list__


class cutomFormatter(
    ap.ArgumentDefaultsHelpFormatter,
    ap.RawDescriptionHelpFormatter,
):
    """Custom Formatter.

    This class combines the argparse ``ArgumentDefaultsHelpFormatter`` and
    ``RawDescriptionHelpFormatter`` formatters.

    """

    pass


def print_message(message):
    """Print Message.

    This method returns a custom argparse action for printing a message.

    Parameters
    ----------
    message : str
        Message to be displayed

    Returns
    -------
    customAction
        Custom action class object

    """
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
            print(message)
            exit()

    return customAction


def module_str():
    """Format Module String.

    Format the list of modules as a single string.

    Returns
    -------
    str
        Formatted string of module names

    """
    string = ''

    for module in __module_list__:
        string += f' - {module}\n'

    return string


def create_arg_parser():
    """Create Argument Parser.

    This method returns an argument parser.

    Returns
    -------
    argparse.Namespace
        Argument parser

    """
    # Create parser
    parser = ap.ArgumentParser(
        add_help=False, description=shapepipe_logo(),
        formatter_class=cutomFormatter,
    )
    optional = parser.add_argument_group('Optional Arguments')

    # Add arguments
    optional.add_argument(
        '-h',
        '--help',
        action='help',
        help='show this help message and exit',
    )

    optional.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'%(prog)s v{__version__}'
    )

    optional.add_argument(
        '-l',
        '--list_modules',
        action=print_message(
            f'ShapePipe modules currently available:\n{module_str()}'
        ),
        help='list modules currently available and exit',
    )

    optional.add_argument(
        '-c',
        '--config',
        default='config.ini',
        help='configuration file name',
    )

    optional.add_argument(
        '-e',
        '--exclusive',
        help='exclusive input file number string',
    )
    # Return parser
    return parser.parse_args()
