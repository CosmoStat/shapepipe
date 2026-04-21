"""Shared pytest fixtures for the repository-level test suite."""

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root.

    Tests live at ``<root>/tests/``, so the root is two parents up
    from this conftest file.
    """
    return Path(__file__).resolve().parent.parent
