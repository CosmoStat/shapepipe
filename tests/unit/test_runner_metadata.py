"""Validate module-runner decorator metadata."""

import importlib
from pathlib import Path

import pytest


RUNNER_PATHS = sorted(Path("src/shapepipe/modules").glob("*_runner.py"))
RUNNER_ATTRIBUTES = [
    "version",
    "input_module",
    "file_pattern",
    "file_ext",
    "depends",
    "executes",
    "numbering_scheme",
    "run_method",
]


def _decorated_functions(module):

    return [
        obj
        for obj in vars(module).values()
        if callable(obj)
        and all(hasattr(obj, attribute) for attribute in RUNNER_ATTRIBUTES)
    ]


@pytest.mark.parametrize(
    "runner_path",
    RUNNER_PATHS,
    ids=lambda path: path.name,
)
def test_runner_module_exposes_decorated_run_function(runner_path):

    try:
        module = importlib.import_module(
            "shapepipe.modules." + runner_path.with_suffix("").name
        )
    except ModuleNotFoundError as err:
        if err.name == "fitsio":
            pytest.skip(
                "fitsio is available in the dev image, not the test extra"
            )
        raise

    if runner_path.name == "module_runners.py":
        return

    runners = _decorated_functions(module)

    assert len(runners) == 1
    runner = runners[0]

    assert runner.run_method in {"parallel", "serial"}
    assert len(runner.file_pattern) == len(runner.file_ext)
