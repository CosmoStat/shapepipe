"""SHAPEPIPE MODULES.

This module contains sub-modules that can be run with ShapePipe.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

import os

# Get a list of all files and directories in the modules directory
modules_dir = os.listdir(os.path.dirname(os.path.abspath(__file__)))
# List all ShapePipe module packages
__all__ = [dir for dir in modules_dir if dir.endswith("_package")]
# List all ShapePipe module runners
__module_list__ = [file for file in modules_dir if file.endswith("_runner.py")]
