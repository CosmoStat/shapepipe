# -*- coding: utf-8 -*-

"""MODULE DECORATOR

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(n_inputs=1, input_module=None, file_pattern='', file_ext=''):
    """ Module Runner Wrapper

    This method adds properties to module runners.

    Parameters
    ----------
    n_inputs : int, optional
        Number of input files, default is 1
    input_module : str, optional
        Input module name, default is None
    file_pattern : str, optional
        File pattern, default is ''
    file_ext : str, optional
        File extension, default is ''

    """

    def decorator(func):

        func.n_inputs = n_inputs
        func.input_module = input_module
        func.file_pattern = file_pattern
        func.file_ext = file_ext

        return func

    return decorator
