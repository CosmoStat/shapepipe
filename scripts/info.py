# -*- coding: utf-8 -*-
  
"""INFO

Set the package information

:Authors: Martin Kilbinger

:Date: 11/10/18

"""

__version__ = "1.0.0"
__author__ = 'Martin Kilbinger, Sam Farrens'
__email__  = 'martin.kilbinger@cea.fr'


# Python dependencies
internal_depend = ['cfis', 'stuff', 'scatalog']

external_depend = ['numpy', 'astropy', 'optparse', 'platform', 'pylab', 'sip_tpv']

__python_depend__ = internal_depend + external_depend

# System dependencies

__system_depend__ = []


