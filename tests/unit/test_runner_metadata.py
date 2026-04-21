"""Every ``*_runner.py`` module must export a function decorated with
``@module_runner``.

The decorator attaches introspection metadata the pipeline needs
(``version``, ``input_module``, ``file_pattern``, ``file_ext``,
``depends``, ``executes``, ``numbering_scheme``, ``run_method``). A new
runner that forgets the decorator loads silently but fails at dispatch
time. This test catches the oversight at import time.

We locate the runner by looking for any top-level function that carries
the metadata, since the function name doesn't always match the file
name (e.g. ``pastecat_runner.py`` exports ``paste_cat_runner``).
"""

import importlib
import inspect
import pkgutil

import pytest

import shapepipe.modules

REQUIRED_ATTRS = (
    "version",
    "input_module",
    "file_pattern",
    "file_ext",
    "depends",
    "executes",
    "numbering_scheme",
    "run_method",
)

# Modules imported here but broken upstream — tracked separately.
# Mapping: module name → xfail reason.
KNOWN_XFAIL = {
    # Both reach stile, which imports treecorr.corr2 (removed in treecorr 5.x).
    "shapepipe.modules.mccd_plots_runner":
        "stile v0.1 imports removed treecorr.corr2",
    "shapepipe.modules.random_cat_runner":
        "stile v0.1 imports removed treecorr.corr2",
}


def _runner_modules() -> list[str]:
    return sorted(
        m.name
        for m in pkgutil.iter_modules(
            shapepipe.modules.__path__, prefix="shapepipe.modules."
        )
        if m.name.endswith("_runner")
    )


def _params():
    for name in _runner_modules():
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
    if "runner_module" in metafunc.fixturenames:
        metafunc.parametrize("runner_module", list(_params()))


def test_runner_has_metadata(runner_module: str) -> None:
    mod = importlib.import_module(runner_module)
    decorated = [
        obj
        for _, obj in inspect.getmembers(mod, inspect.isfunction)
        if obj.__module__ == runner_module
        and all(hasattr(obj, a) for a in REQUIRED_ATTRS)
    ]
    assert decorated, (
        f"{runner_module} exports no function decorated with @module_runner "
        f"(required attrs: {REQUIRED_ATTRS})"
    )
