"""SHAPEPIPE PACKAGE.

ShapePipe is a galaxy shape measurement pipeline.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

__all__ = ['modules', 'pipeline', 'utilities']

from importlib.metadata import metadata, version
from . import *

__version__ = version("shapepipe")
__about__ = metadata("shapepipe").get("Summary")