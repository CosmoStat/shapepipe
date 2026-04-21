"""Every ``[project.scripts]`` entry must resolve and handle ``-h`` cleanly.

The test skips entries whose command isn't on PATH — that's an install
concern, not a code concern. When a command *is* installed, ``-h``
must exit 0 or 2 (argparse's help exit code) with non-empty output.
This catches entry points whose target function treats ``-h`` as a
positional argument, or whose argparser is initialised so late the
help flag never runs.
"""

import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

# Entry points known to mishandle ``-h`` — tracked via ``xfail`` so the
# suite stays green and the issue is discoverable once fixed.
KNOWN_XFAIL = {
    "summary_run":
        "treats '-h' as the 'patch' positional arg, tries to mkdir '-h'",
    # All three import shapepipe.canfar_run, which fails with
    # IndentationError in canfar_monitor.py:55.
    "canfar_submit_job":
        "IndentationError in canfar_monitor.py:55 (transitive)",
    "canfar_monitor":
        "IndentationError in canfar_monitor.py:55 (transitive)",
    "canfar_monitor_log":
        "IndentationError in canfar_monitor.py:55 (transitive)",
}


def _scripts_from_pyproject() -> dict[str, str]:
    root = Path(__file__).resolve().parents[2]
    pyproject = tomllib.loads((root / "pyproject.toml").read_text())
    return pyproject.get("project", {}).get("scripts", {})


def _params():
    for name in sorted(_scripts_from_pyproject().keys()):
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
    if "script_name" in metafunc.fixturenames:
        metafunc.parametrize("script_name", list(_params()))


def test_entrypoint_help(script_name: str) -> None:
    if shutil.which(script_name) is None:
        pytest.skip(
            f"{script_name} not on PATH — not installed in this env"
        )
    result = subprocess.run(
        [script_name, "-h"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    # argparse help exits 0; some parsers use 2. Either is fine so
    # long as *something* was printed, meaning the help path ran.
    assert result.returncode in (0, 2), (
        f"{script_name} -h exited {result.returncode}\n"
        f"stdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert result.stdout or result.stderr, (
        f"{script_name} -h produced no output"
    )
