# -*- coding: utf-8 -*-

"""MODULE DECORATOR

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(n_inputs=1, input_module=None, ext=None):
    """ Module Runner Wrapper

    This method adds properties to module runners.

    Parameters
    ----------
    n_inputs : int, optional
        Number of input files, default is 1
    input_module : str, optional
        Input module name, default is None
    ext : str, optional
        File extension, default is None

    """

    def decorator(func):

        func.n_inputs = n_inputs
        func.input_module = input_module
        func.ext = ext

        return func

    return decorator
