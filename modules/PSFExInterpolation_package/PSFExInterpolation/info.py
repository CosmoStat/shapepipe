# -*- coding: utf-8 -*-

"""INFO

Set the package information

:Authors: Samuel Farrens

:Date: 16/04/2018
"""

__version__ = "1.1.1"
__whoami__ = "PSFExInterpolation"

# Python dependencies
pipe_depend = ['mpfg', 'mpfx', 'scatalog', 'sconfig', 'shapepipe_base',
               'slogger']
external_depend = ['numpy', 'galsim', 'astropy', 'psfex']

__python_depend__ = pipe_depend + external_depend

# System dependencies

__system_depend__ = []
