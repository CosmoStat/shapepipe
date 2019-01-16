# -*- coding: utf-8 -*-

"""MODULE TOOLS

This module defines methods for handling pipeline modules.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from inspect import getmembers, isfunction
from shapepipe.modules import module_runners


def get_module_version(module):
    """ Get Module Version

    This method returns the version for a given module.

    Parameters
    ----------
    module : str
        Module name

    Returns
    -------
    str
        Module version

    """

    if not isinstance(module, str):
        raise TypeError('Module name must be a string.')

    return getattr(module_runners, module).version


def get_module_list():
    """ Get Module List

    This method returns a list of the modules current available in
    module_runners.

    Returns
    -------
    list
        List of module names

    """

    return list(zip(*getmembers(module_runners, isfunction)))[0]


def get_module_depends():
    """ Get Module Dependencies

    This method returns a list of the module dependencies.

    Returns
    -------
    list
        List of module dependencies

    """

    depends = []
    executes = []

    for module in get_module_list():
        depends += getattr(module_runners, module).depends
        executes += getattr(module_runners, module).executes

    return depends, executes
