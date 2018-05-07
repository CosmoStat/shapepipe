# -*- coding: utf-8 -*-

"""INFO

Set the package information

:Authors: Samuel Farrens

:Date: 10/04/2018
"""

__version__ = "1.0.0"
__whoami__ = "MER"

# Python dependencies
pipe_depend = ['mpfg', 'mpfx', 'scatalog', 'sconfig', 'shapepipe_base',
               'slogger']
external_depend = ['numpy']

__python_depend__ = pipe_depend + external_depend
