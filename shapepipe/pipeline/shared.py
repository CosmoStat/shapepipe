# -*- coding: utf-8 -*-

"""SHARED

The module defines functions that can be shared between pipeline modules.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def split_module_run(module_str):
    """Split Module Run

    Extract module name and run from input string.

    Parameters
    ----------
    module_str : str
        Module name or run string

    Returns
    -------
    tuple
        Module name and module run string

    Raises
    ------
    TypeError
        If input is not a string

    """
    if not isinstance(module_str, str):
        raise TypeError(
            f'Input module_str must be a string not {type(module_str)}.'
        )

    run_split = '_run_'
    module_run = module_str

    if run_split in module_str:
        module_name = module_str.split(run_split)[0]
    else:
        module_name = module_str

    return module_name, module_run
