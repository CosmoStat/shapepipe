# -*- coding: utf-8 -*-

"""SHAPEPIPE INFO

This module provides some basic information about the ShapePipe package.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

try:
    from termcolor import colored
except ImportError:
    import_fail = True
else:
    import_fail = False


# Package Info
version_info = (0, 0, 2)
__version__ = '.'.join(str(c) for c in version_info)
__name__ = 'shapepipe'
__author__ = 'Samuel Farrens'
__email__ = 'samuel.farrens@cea.fr'
__about__ = ('ShapePipe is a shape measurement pipeline developed with the'
             'CosmoStat lab at CEA Paris-Saclay.')
__setups__ = ['pytest-runner']
__installs__ = ['joblib>=0.13',
                'modopt>=1.2',
                'numpy>=1.14']
__tests__ = ['pytest',
             'pytest-cov',
             'pytest-pycodestyle']


def shapepipe_logo(colour=False):
    """ ShapePipe Logo

    Returns
    -------
    str logo string

    """

    shape = r'''
 _______  __   __  _______  _______  _______  _______  ___   _______  _______
|       ||  | |  ||   _   ||       ||       ||       ||   | |       ||       |
|  _____||  |_|  ||  |_|  ||    _  ||    ___||    _  ||   | |    _  ||    ___|
| |_____ |       ||       ||   |_| ||   |___ |   |_| ||   | |   |_| ||   |___
|_____  ||       ||       ||    ___||    ___||    ___||   | |    ___||    ___|
 _____| ||   _   ||   _   ||   |    |   |___ |   |    |   | |   |    |   |___
|_______||__| |__||__| |__||___|    |_______||___|    |___| |___|    |_______|
    '''

    if not import_fail and colour:
        shape = colored(shape, 'cyan', attrs=['bold'])

    logo = r'''
-------------------------------------------------------------------------------
{}


    Shape measurement pipeline developed at CosmoStat.

    Authors: Samuel Farrens   <samuel.farrens@cea.fr>
             Axel Guinot      <axel.guinot@cea.fr>

    Main Contributors:
             Martin Kilbinger
             Tobias Liaudat
             Morgan Schmitz

    Version: {}

-------------------------------------------------------------------------------
    '''.format(shape, __version__)

    return logo


def line():
    """ Line

    Returns
    -------
    str a horizontal line

    """

    line = r'''
-------------------------------------------------------------------------------
    '''

    return line
