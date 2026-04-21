"""Every example INI config must parse as a valid configparser file.

Uses ``RawConfigParser`` so ``$SP_RUN``-style variable references don't
trigger interpolation errors. This catches malformed section headers,
duplicate keys within a section, and missing ``=``/``:`` separators.
"""

import configparser
from pathlib import Path

import pytest


def _config_files(root: Path) -> list[Path]:
    return sorted((root / "example").rglob("*.ini"))


def pytest_generate_tests(metafunc):
    if "config_path" in metafunc.fixturenames:
        root = Path(__file__).resolve().parents[2]
        configs = _config_files(root)
        metafunc.parametrize(
            "config_path",
            configs,
            ids=[str(p.relative_to(root)) for p in configs],
        )


def test_config_parses(config_path: Path) -> None:
    parser = configparser.RawConfigParser(strict=True)
    parser.read(config_path, encoding="utf-8")
    assert parser.sections(), (
        f"{config_path} parsed but has no sections — likely malformed"
    )
