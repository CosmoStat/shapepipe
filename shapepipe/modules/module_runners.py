"""MODULE RUNNERS.

This module defines methods for running the pipeline modules.

:Author: Samuel Farrens <samuel.farrens@cea.fr>

"""

from importlib import import_module


def get_module_runners(modules):
    """Get Module Runners.

    Import the specified module runners.

    Parameters
    ----------
    modules : list
        List of module names

    Returns
    -------
    dict
        Dictionary of module runners

    """
    package = "shapepipe.modules"

    module_runners = dict(
        [
            (
                module,
                getattr(
                    import_module(f".{module}", package=package),
                    module,
                ),
            )
            for module in modules
        ]
    )

    return module_runners
