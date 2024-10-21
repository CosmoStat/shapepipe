"""MODULE DECORATOR.

This module defines the module runner dectorator.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""


def module_runner(
    version="0.0",
    input_module=None,
    file_pattern="",
    file_ext="",
    depends=[],
    executes=[],
    numbering_scheme=None,
    run_method="parallel",
):
    """Wrap Module Runners.

    This method adds properties to module runners.

    Parameters
    ----------
    version : str, optional
        Module version string, default is ``0.0``
    input_module : str or list, optional
        Input module name, default is ``None``
    file_pattern : str or list, optional
        File pattern, default is ``''``
    file_ext : str or list, optional
        File extension, default is ``''``
    depends : str or list, optional
        Module dependencies, default is ``[]``
    executes : str or list, optional
        Module executables, default is ``[]``
    numbering_scheme : str, optional
        Module numbering scheme, default is ``None``
    run_method : {'parallel', 'serial'}, optional
        Module run method, default is ``'parallel'``

    Raises
    ------
    TypeError
        If ``version`` is not a string
    TypeError
        If ``input_module`` is not a list or a string
    TypeError
        If ``file_pattern`` is not a list or a string
    TypeError
        If ``file_ext`` is not a list or a string
    TypeError
        If ``depends`` is not a list or a string
    TypeError
        If ``executes`` is not a list or a string
    TypeError
        If ``numbering_scheme`` is not a string
    ValueError
        If ``run_method`` is not valid
    ValueError
        If length of ``file_pattern`` and ``file_ext`` lists do not match

    """
    if not isinstance(version, str):
        raise TypeError("Module version must be a string.")

    if isinstance(input_module, str):
        input_module = [input_module]
    elif not isinstance(input_module, (list, type(None))):
        raise TypeError("Input module must be a list or a string.")

    if isinstance(file_pattern, str):
        file_pattern = [file_pattern]
    elif not isinstance(file_pattern, list):
        raise TypeError("File pattern must be a string or a list of strings.")

    if isinstance(file_ext, str):
        file_ext = [file_ext]
    elif not isinstance(file_ext, list):
        raise TypeError("File extension must be a string or a list of strings.")

    if isinstance(depends, str):
        depends = [depends]
    elif not isinstance(depends, list):
        raise TypeError("Dependencies must be a string or a list of strings.")

    if isinstance(executes, str):
        executes = [executes]
    elif not isinstance(executes, list):
        raise TypeError("Executables must be a string or a list of strings.")

    if not isinstance(numbering_scheme, (str, type(None))):
        raise TypeError(
            f"Numbering scheme must be a string, found '{numbering_scheme}'."
        )

    if run_method not in {"parallel", "serial"}:
        raise ValueError('Run method must be a "parallel" or "serial".')

    if (len(file_ext) == 1) and (len(file_pattern) > 1):
        file_ext = [file_ext[0] for i in file_pattern]

    if len(file_ext) != len(file_pattern):
        raise ValueError(
            f"The number of file_ext values ({len(file_ext)}) does not match "
            + f"the number of file_pattern values ({len(file_pattern)}) in "
            + "the module decorator."
        )

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
