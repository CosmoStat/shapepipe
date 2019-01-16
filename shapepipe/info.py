# -*- coding: utf-8 -*-

"""SHAPEPIPE INFO

This module provides some basic information about the ShapePipe package.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

# Package Info
version_info = (0, 0, 1)
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
             'pytest-cov']


def shapepipe_logo():
    """ ShapePipe Logo

    Returns
    -------
    str logo string

    """

    logo = r'''
-------------------------------------------------------------------------------
 _______  __   __  _______  _______  _______  _______  ___   _______  _______
|       ||  | |  ||   _   ||       ||       ||       ||   | |       ||       |
|  _____||  |_|  ||  |_|  ||    _  ||    ___||    _  ||   | |    _  ||    ___|
| |_____ |       ||       ||   |_| ||   |___ |   |_| ||   | |   |_| ||   |___
|_____  ||       ||       ||    ___||    ___||    ___||   | |    ___||    ___|
 _____| ||   _   ||   _   ||   |    |   |___ |   |    |   | |   |    |   |___
|_______||__| |__||__| |__||___|    |_______||___|    |___| |___|    |_______|


    Shape measurement pipeline developed at CosmoStat.

    Author: Samuel Farrens <samuel.farrens@cea.fr>

    Version: {}

-------------------------------------------------------------------------------
    '''.format(__version__)

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
