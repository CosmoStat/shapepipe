"""Smoke-test example and workflow configuration files."""

import ast
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


def _required_config_keys(module_name):
    """Config keys a runner reads unconditionally from its own section.

    Static approximation of "unconditional": a ``config.get*(module_config_sec,
    "KEY")`` call that is a direct statement of the runner function's body, not
    nested inside an ``if``/``for``/``while``/``try``/``with`` — a key read only
    inside one of those (an ``if config.has_option(...)`` guard, an
    ``if some_flag:`` branch with a fallback default, ...) is, by construction,
    not required by every config. Returns ``None`` if the runner file or its
    decorated function cannot be found (callers then skip the check for it).
    """
    path = Path("src/shapepipe/modules") / f"{module_name}.py"
    if not path.is_file():
        return None

    tree = ast.parse(path.read_text())
    func = next(
        (
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.FunctionDef) and node.name == module_name
        ),
        None,
    )
    if func is None:
        return None

    required = set()
    for stmt in func.body:
        if isinstance(stmt, (ast.If, ast.For, ast.While, ast.Try, ast.With)):
            continue
        for node in ast.walk(stmt):
            if not isinstance(node, ast.Call):
                continue
            callee = node.func
            if not (
                isinstance(callee, ast.Attribute)
                and isinstance(callee.value, ast.Name)
                and callee.value.id == "config"
                and callee.attr.startswith("get")
            ):
                continue
            has_fallback = any(kw.arg == "fallback" for kw in node.keywords)
            if (
                not has_fallback
                and len(node.args) >= 2
                and isinstance(node.args[1], ast.Constant)
                and isinstance(node.args[1].value, str)
            ):
                required.add(node.args[1].value)

    return required


def _module_config_sections(config_path):
    """Module name -> config section name, for every module in MODULE.

    Mirrors FileHandler.set_up_module / get_module_config_sec: a module named
    once uses its own upper-cased name as its section; a module repeated in
    MODULE gets ``<module>_run_<n>`` (1-indexed) upper-cased instead. A module
    whose name is not resolvable statically (e.g. ``${SP_PSF}_interp_runner``)
    is skipped.
    """
    parser = configparser.ConfigParser()
    parser.read(config_path)

    if not parser.has_option("EXECUTION", "MODULE"):
        # Not a ShapePipe pipeline config (e.g. config_MCCD.ini is the MCCD
        # library's own params file, read by mccd_preprocessing_runner/
        # mccd_fit_val_runner via CONFIG_PATH, not by ShapePipe's [EXECUTION]).
        return {}

    raw = parser.get("EXECUTION", "MODULE")
    modules = [m.strip() for m in raw.replace("\n", ",").split(",") if m.strip()]

    counts = {}
    sections = {}
    for module in modules:
        if "$" in module:
            continue
        counts[module] = counts.get(module, 0) + 1
        run_name = (
            module if modules.count(module) == 1 else f"{module}_run_{counts[module]}"
        )
        sections[module] = run_name.upper()

    return sections


@pytest.mark.parametrize("config_path", WORKFLOW_CONFIG_FILES, ids=str)
def test_workflow_config_module_sections_have_required_keys(config_path):
    """Every module in MODULE has its config section, with the keys it needs.

    Catches the same defect class as the runner-name check by a different
    route: a config whose MODULE chain names a real runner, but whose section
    is missing (wrong name) or missing a key that runner reads unconditionally,
    also fails at pipeline start-up.
    """
    parser = configparser.ConfigParser()
    parser.read(config_path)

    problems = []
    for module, section in _module_config_sections(config_path).items():
        required = _required_config_keys(module)
        if required is None:
            continue

        if not parser.has_section(section):
            problems.append(f"{module}: missing section [{section}]")
            continue

        missing = sorted(
            key for key in required if not parser.has_option(section, key)
        )
        if missing:
            problems.append(f"[{section}] ({module}) missing key(s): {missing}")

    assert not problems, f"{config_path}: " + "; ".join(problems)
