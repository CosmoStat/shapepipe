"""Every public ``shapepipe.*`` submodule must import cleanly.

Catches module-level syntax errors, typos in import paths, circular
imports, and (crucially) the pattern where a module unconditionally
imports an optional dependency — a ``try/except ImportError`` guard
placed *after* a failing top-level import is unreachable.

Assumes the shapepipe package and its declared dependencies are
installed (the official container, or ``pip install -e .``).
"""

import importlib
import pkgutil

import pytest

import shapepipe

# Modules with known-broken transitive imports — tracked so the suite
# lands green but the failure is discoverable and auto-notifies (via
# ``strict=True``) once the upstream fix lands.
KNOWN_XFAIL = {
    # IndentationError at canfar_monitor.py:55 (docstring/body mismatch).
    "shapepipe.canfar.canfar_monitor":
        "IndentationError in canfar_monitor.py:55",
    "shapepipe.canfar_run":
        "IndentationError in canfar_monitor.py:55 (transitive)",
    # stile v0.1 hard-imports treecorr.corr2, which newer treecorr removed.
    "shapepipe.modules.mccd_package.mccd_plot_utilities":
        "stile v0.1 imports removed treecorr.corr2",
    "shapepipe.modules.mccd_plots_runner":
        "stile v0.1 imports removed treecorr.corr2",
    "shapepipe.modules.random_cat_package.random_cat":
        "stile v0.1 imports removed treecorr.corr2",
    "shapepipe.modules.random_cat_runner":
        "stile v0.1 imports removed treecorr.corr2",
    # Modules added in v1.x Dockerfile rewrites that import deps not in
    # the develop image (astroquery, numba, fitsio) or removed from
    # stdlib/setuptools (pkg_resources).
    "shapepipe.modules.mask_package.mask":
        "astroquery not in develop Docker image",
    "shapepipe.modules.mask_runner":
        "astroquery not in develop Docker image (transitive)",
    "shapepipe.modules.ngmix_package.ngmix":
        "numba not in develop Docker image",
    "shapepipe.modules.ngmix_runner":
        "numba not in develop Docker image (transitive)",
    "shapepipe.modules.uncompress_fits_package.uncompress_fits":
        "fitsio not installed (declared as optional extra)",
    "shapepipe.modules.uncompress_fits_runner":
        "fitsio not installed (transitive)",
}


def _iter_shapepipe_modules() -> list[str]:
    return sorted(
        m.name
        for m in pkgutil.walk_packages(
            shapepipe.__path__, prefix="shapepipe."
        )
        if "tests" not in m.name.split(".")
    )


def _params():
    for name in _iter_shapepipe_modules():
        if name in KNOWN_XFAIL:
            yield pytest.param(
                name,
                marks=pytest.mark.xfail(
                    reason=KNOWN_XFAIL[name], strict=True, raises=Exception
                ),
            )
        else:
            yield name


def pytest_generate_tests(metafunc):
    if "module_name" in metafunc.fixturenames:
        metafunc.parametrize("module_name", list(_params()))


def test_module_imports(module_name: str) -> None:
    importlib.import_module(module_name)
