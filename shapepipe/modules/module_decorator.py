# -*- coding: utf-8 -*-

"""MODULE DECORATOR

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(input_module=None, version='0.0', file_pattern='',
                  file_ext='', depends=[], executes=[]):
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

    if not isinstance(input_module, (str, type(None))):
        raise TypeError('Module name must be a string.')

    if not isinstance(version, str):
        raise TypeError('Module version must be a string.')

    if isinstance(file_pattern, str):
        file_pattern = [file_pattern]
    elif not isinstance(file_pattern, list):
        raise TypeError('File pattern must be a string or a list of strings')

    if isinstance(file_ext, str):
        file_ext = [file_ext]
    elif not isinstance(file_ext, list):
        raise TypeError('File extension must be a string or a list of strings')

    if isinstance(depends, str):
        depends = [depends]
    elif not isinstance(depends, list):
        raise TypeError('Dependencies must be a string or a list of strings')

    if isinstance(executes, str):
        executes = [executes]
    elif not isinstance(depends, list):
        raise TypeError('Executables must be a string or a list of strings')

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
        func.version = version
        func.file_pattern = file_pattern
        func.file_ext = file_ext
        func.depends = depends
        func.executes = executes

        return func

    return decorator
