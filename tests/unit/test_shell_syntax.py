"""Every shell script under scripts/ must parse under ``bash -n``.

Catches unclosed quotes, unbalanced brackets, stray parentheses, and the
whole family of errors that only surface the moment someone runs the
script. Parametrized one-test-per-file so CI output names the offender.
"""

import subprocess
from pathlib import Path

import pytest


def _shell_scripts(root: Path) -> list[Path]:
    """Return every ``.sh``/``.bash`` under ``scripts/``."""
    scripts_dir = root / "scripts"
    return sorted(
        p for p in scripts_dir.rglob("*")
        if p.is_file() and p.suffix in {".sh", ".bash"}
    )


def pytest_generate_tests(metafunc):
    """Parametrize per-script at collection time so we can name each file."""
    if "script_path" in metafunc.fixturenames:
        root = Path(__file__).resolve().parents[2]
        scripts = _shell_scripts(root)
        metafunc.parametrize(
            "script_path",
            scripts,
            ids=[str(p.relative_to(root)) for p in scripts],
        )


def test_shell_script_parses(script_path: Path) -> None:
    """bash -n must succeed."""
    result = subprocess.run(
        ["bash", "-n", str(script_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"bash -n failed on {script_path}:\n{result.stderr}"
    )
