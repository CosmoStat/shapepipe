"""Import every package submodule."""

import importlib
import pkgutil

import pytest
import shapepipe


def _submodules():

    return sorted(
        module_info.name
        for module_info in pkgutil.walk_packages(
            shapepipe.__path__,
            prefix=f"{shapepipe.__name__}.",
        )
    )


@pytest.mark.parametrize("module_name", _submodules())
def test_submodule_imports(module_name):

    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError as err:
        if err.name == "fitsio":
            pytest.skip(
                "fitsio is available in the dev image, not the test extra"
            )
        raise
