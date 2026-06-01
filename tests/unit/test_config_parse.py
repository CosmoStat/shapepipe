"""Smoke-test example configuration files."""

import configparser
from pathlib import Path

import pytest


CONFIG_FILES = sorted(Path("example").glob("**/*.ini"))


@pytest.mark.parametrize("config_path", CONFIG_FILES, ids=str)
def test_example_config_is_parseable(config_path):

    parser = configparser.ConfigParser()

    assert parser.read(config_path) == [str(config_path)]
    assert parser.sections()
