"""SHARED.

The module defines functions that can be shared between pipeline modules.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from glob import glob


def check_duplicate(input_list):
    """Check Duplicate.

    Check whether input list contains at least one duplicate.

    Parameters
    ----------
    input_list : list
        input list

    Returns
    -------
    str
        Duplicate element, empty string if none found

    """
    input_set = set()

    for elem in input_list:
        if elem in input_set:
            return elem
        else:
            input_set.add(elem)

    return ""


def find_files(path, pattern="*", ext="*"):
    """Find Files.

    This method recursively retrieves file names from a given path that
    match a given pattern and/or have a given extension.

    Parameters
    ----------
    path : str
        Full path to files
    pattern : str, optional
        File pattern, default is '*'
    ext : str, optional
        File extension, default is '*'

    Returns
    -------
    list
        List of file names

    Raises
    ------
    ValueError
        For '*' in pattern
    ValueError
        For '*' in extension
    ValueError
        For invalid extension format

    """
    dot = "."
    star = "*"

    if pattern != star and star in pattern:
        raise ValueError('Do not include "*" in pattern.')

    if ext != star and star in ext:
        raise ValueError('Do not include "*" in extension.')

    if (not ext.startswith(dot) and dot in ext) or (ext.count(dot) > 1):
        raise ValueError(f'Invalid extension format: "{ext}".')

    if ext != star and not ext.startswith(dot):
        ext = dot + ext

    search_string = f"{path}/**/*{pattern}*{ext}"

    return glob(search_string, recursive=True)


def split_module_run(module_str):
    """Split Module Run.

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
            f"Input module_str must be a string not {type(module_str)}."
        )

    run_split = "_run_"
    module_run = module_str

    if run_split in module_str:
        module_name = module_str.split(run_split)[0]
    else:
        module_name = module_str

    return module_name, module_run
