"""Smoke-test bundled shell scripts for parse errors."""

import subprocess
from pathlib import Path

import pytest


SHELL_SCRIPTS = sorted(Path("scripts/sh").glob("*"))


@pytest.mark.parametrize("script_path", SHELL_SCRIPTS, ids=str)
def test_shell_script_parses(script_path):

    subprocess.run(["bash", "-n", script_path], check=True)
