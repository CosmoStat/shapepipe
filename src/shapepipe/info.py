"""SHAPEPIPE INFO.

This module provides some basic information about the ShapePipe package.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

try:
    from termcolor import colored
except ImportError:
    import_fail = True
else:
    import_fail = False
    
from shapepipe import __version__

def shapepipe_logo(colour=False):
    """Get ShapePipe Logo.

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
             Martin Kilbinger <martin.kilbinger@cea.fr>

    Main Contributors:
             Tobias Liaudat
             Morgan Schmitz
             Andre Zamorano Vitorelli
             Francois Lanusse
             Xavier Jimenez

    Version: {}

-------------------------------------------------------------------------------
    '''.format(shape, __version__)

    return logo


def line():
    """Get Horizontal Line.

    Returns
    -------
    str a horizontal line

    """
    line = r'''
-------------------------------------------------------------------------------
    '''

    return line
