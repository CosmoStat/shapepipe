# -*- coding: utf-8 -*-

"""MODULE DECORATOR

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(input_module=None, version='0.0', file_pattern='',
                  file_ext='', depends=[], executes=[], numbering_scheme=None,
                  run_method='parallel'):
    """ Module Runner Wrapper

    This method adds properties to module runners.

    Parameters
    ----------
    input_module : str or list, optional
        Input module name, default is None
    version : str, optional
        Module version string
    file_pattern : str or list, optional
        File pattern, default is ''
    file_ext : str or list, optional
        File extension, default is ''
    depends : str or list, optional
        Module dependencies, default is []
    executes : str or list, optional
        Module executables, default is []
    numbering_scheme : str, optional
        Module numbering scheme, default is None

    """

    if isinstance(input_module, str):
        input_module = [input_module]
    elif not isinstance(input_module, (list, type(None))):
        raise TypeError('Input module must be a list or a string.')

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
    elif not isinstance(executes, list):
        raise TypeError('Executables must be a string or a list of strings')

    if not isinstance(numbering_scheme, (str, type(None))):
        raise TypeError('Numbering scheme must be a string, found \'{}\'.'
                        ''.format(numbering_scheme))

    if (len(file_ext) == 1) and (len(file_pattern) > 1):
        file_ext = [file_ext[0] for i in file_pattern]

    if len(file_ext) != len(file_pattern):
        raise ValueError('The number of file_ext values ({}) does not match '
                         'the number of file_pattern values ({}) in the '
                         'module decorator.'
                         ''.format(len(file_ext), len(file_pattern)))

    def decorator(func):

        func.input_module = input_module
        func.version = version
        func.file_pattern = file_pattern
        func.file_ext = file_ext
        func.depends = depends
        func.executes = executes
        func.numbering_scheme = numbering_scheme
        func.run_method = run_method

        return func

    return decorator
