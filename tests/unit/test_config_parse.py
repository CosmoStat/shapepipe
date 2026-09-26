"""Smoke-test example and workflow configuration files."""

import configparser
import re
from pathlib import Path

import pytest


CONFIG_FILES = sorted(Path("example").glob("**/*.ini"))

WORKFLOW_CONFIG_FILES = sorted(Path("workflow/config/cfis").glob("*.ini"))

# Every module runner ShapePipe ships, by its bare name (e.g. "mask_query_runner"),
# as named by its file under src/shapepipe/modules/.
RUNNER_NAMES = {
    path.stem for path in Path("src/shapepipe/modules").glob("*_runner.py")
}

# A bare "<name>_runner" token anywhere in a config value: MODULE lists,
# INPUT_MODULE (with or without a "last:" prefix), a runner name embedded in an
# INPUT_DIR path (".../<runner>/output"), or a *_RUNNER option pointing at
# another runner (e.g. ME_DOT_PSF_RUNNER).
RUNNER_TOKEN_RE = re.compile(r"\b[a-z][a-z0-9]*(?:_[a-z0-9]+)*_runner\b")


@pytest.mark.parametrize("config_path", CONFIG_FILES, ids=str)
def test_example_config_is_parseable(config_path):

    parser = configparser.ConfigParser()

    assert parser.read(config_path) == [str(config_path)]
    assert parser.sections()


@pytest.mark.parametrize("config_path", WORKFLOW_CONFIG_FILES, ids=str)
def test_workflow_config_is_parseable(config_path):

    parser = configparser.ConfigParser()

    assert parser.read(config_path) == [str(config_path)]
    assert parser.sections()


@pytest.mark.parametrize("config_path", WORKFLOW_CONFIG_FILES, ids=str)
def test_workflow_config_runner_names_exist(config_path):
    """Every ``*_runner`` token in a workflow config names a real runner.

    A config that names a deleted or misspelled runner in MODULE, an
    INPUT_MODULE, an INPUT_DIR path, or a ``*_RUNNER`` option fails to run at
    pipeline start-up; catching it here does not require a pipeline run.
    """
    text = config_path.read_text()
    tokens = set(RUNNER_TOKEN_RE.findall(text))

    unknown = tokens - RUNNER_NAMES

    assert not unknown, (
        f"{config_path} references runner(s) not found under "
        f"src/shapepipe/modules/: {sorted(unknown)}"
    )
