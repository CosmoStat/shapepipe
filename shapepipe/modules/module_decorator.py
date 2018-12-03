# -*- coding: utf-8 -*-

"""MODULE DECORATOR

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(input_module=None, file_pattern='', file_ext=''):
    """ Module Runner Wrapper

    This method adds properties to module runners.

    Parameters
    ----------
    input_module : str, optional
        Input module name, default is None
    file_pattern : str or list, optional
        File pattern, default is ''
    file_ext : str or list, optional
        File extension, default is ''

    """

    if isinstance(file_pattern, str):
        file_pattern = [file_pattern]

    if isinstance(file_ext, str):
        file_ext = [file_ext]

    if (len(file_ext) == 1) and (len(file_pattern) > 1):
        file_ext = [file_ext[0] for i in file_pattern]
    elif (len(file_pattern) == 1) and (len(file_ext) > 1):
        file_pattern = [file_pattern[0] for i in file_ext]

    if len(file_ext) != len(file_pattern):
        raise ValueError('The number of file_ext values does not match the '
                         'number of file_pattern values in the module '
                         'decorator.')

    def decorator(func):

        func.input_module = input_module
        func.file_pattern = file_pattern
        func.file_ext = file_ext

        return func

    return decorator
