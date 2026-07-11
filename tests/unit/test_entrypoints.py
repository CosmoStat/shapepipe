"""Smoke-test public console scripts."""

import subprocess
import sys
import tomllib
from pathlib import Path

import pytest


with Path("pyproject.toml").open("rb") as pyproject:
    ENTRYPOINTS = sorted(
        tomllib.load(pyproject)["project"]["scripts"].items(),
    )


@pytest.mark.parametrize(
    "entrypoint",
    ENTRYPOINTS,
    ids=lambda entrypoint: entrypoint[0],
)
def test_console_entrypoint_help_runs(entrypoint):

    _, target = entrypoint
    module_name, _, function_name = target.partition(":")
    command = (
        [sys.executable, "-m", module_name, "-h"]
        if function_name == "main"
        else [
            sys.executable,
            "-c",
            (
                f"from {module_name} import {function_name}; "
                f"{function_name}()"
            ),
            "-h",
        ]
    )

    result = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout.lower()


@pytest.mark.parametrize("exit_code", [1, None])
def test_main_propagates_run_exit_code(monkeypatch, exit_code):
    """``main`` must forward ``run``'s return value.

    ``run`` returns 1 when it catches an error; if ``main`` drops that,
    ``exit(main())`` becomes ``exit(0)`` and every handled failure looks like
    success to the batch system.
    """
    import shapepipe.shapepipe_run as entry

    monkeypatch.setattr(entry, "run", lambda args=None: exit_code)

    assert entry.main() == exit_code
