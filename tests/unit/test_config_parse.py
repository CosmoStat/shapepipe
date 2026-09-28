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


class _ConfigReadVisitor(ast.NodeVisitor):
    """Find a runner function's reads of its own config section.

    Heuristic, not exhaustive: it flags a ``config.get*(module_config_sec,
    "KEY")`` call as unconditionally required only if it sits outside every
    conditional construct (``if``/``for``/``while``/``try``/``with``, and a
    ternary ``... if ... else ...``) and takes no ``fallback=`` kwarg — and
    then drops it anyway if a ``config.has_option(module_config_sec, "KEY")``
    check on the same key appears anywhere in the function, at any nesting
    (this also covers an early-return guard, where the check precedes the read
    as a sibling statement rather than wrapping it). It only ever looks at
    calls whose section argument is literally the ``module_config_sec``
    parameter; a read of some other section is ignored. A false negative here
    (a genuinely required key missed) is expected and safe — it just narrows
    what the test can catch; a false positive would break unrelated configs,
    which the guards above are meant to prevent.
    """

    def __init__(self):
        self.unconditional = set()
        self.guarded_anywhere = set()
        self._cond_depth = 0

    def _visit_conditional(self, node):
        self._cond_depth += 1
        self.generic_visit(node)
        self._cond_depth -= 1

    visit_If = _visit_conditional
    visit_For = _visit_conditional
    visit_AsyncFor = _visit_conditional
    visit_While = _visit_conditional
    visit_Try = _visit_conditional
    visit_With = _visit_conditional
    visit_AsyncWith = _visit_conditional
    visit_IfExp = _visit_conditional

    def visit_Call(self, node):
        callee = node.func
        is_config_call = (
            isinstance(callee, ast.Attribute)
            and isinstance(callee.value, ast.Name)
            and callee.value.id == "config"
        )
        if is_config_call and len(node.args) >= 2:
            sec_arg, key_arg = node.args[0], node.args[1]
            on_module_sec = (
                isinstance(sec_arg, ast.Name) and sec_arg.id == "module_config_sec"
            )
            is_key_literal = isinstance(key_arg, ast.Constant) and isinstance(
                key_arg.value, str
            )
            if on_module_sec and is_key_literal:
                key = key_arg.value
                if callee.attr == "has_option":
                    self.guarded_anywhere.add(key)
                elif callee.attr.startswith("get"):
                    has_fallback = any(kw.arg == "fallback" for kw in node.keywords)
                    if not has_fallback and self._cond_depth == 0:
                        self.unconditional.add(key)

        self.generic_visit(node)


def _required_config_keys(module_name):
    """Config keys a runner's source reads unconditionally (see the visitor).

    Returns ``None`` if the runner file or its decorated function cannot be
    found (callers then skip the check for it).
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

    visitor = _ConfigReadVisitor()
    visitor.visit(func)

    return visitor.unconditional - visitor.guarded_anywhere


def _module_config_sections(raw_module_list):
    """(module name, config section name) for every invocation in MODULE.

    Mirrors FileHandler.set_up_module / get_module_config_sec: a module named
    once uses its own upper-cased name as its section; a module repeated in
    MODULE gets ``<module>_run_<n>`` (1-indexed) upper-cased instead — every
    invocation is returned, so a module used twice (e.g. ``vignetmaker_runner``
    for ``_RUN_1`` and ``_RUN_2``) yields two pairs, not one. A module whose
    name is not resolvable statically (e.g. ``${SP_PSF}_interp_runner``) is
    skipped.
    """
    modules = [
        m.strip() for m in raw_module_list.replace("\n", ",").split(",") if m.strip()
    ]

    counts = {}
    pairs = []
    for module in modules:
        if "$" in module:
            continue
        counts[module] = counts.get(module, 0) + 1
        run_name = (
            module if modules.count(module) == 1 else f"{module}_run_{counts[module]}"
        )
        pairs.append((module, run_name.upper()))

    return pairs


def _missing_required_keys(parser):
    """Every problem `test_workflow_config_module_sections_have_required_keys`
    checks for, given an already-populated ``ConfigParser``. Returns a list of
    human-readable problem strings (empty if none).
    """
    if not parser.has_option("EXECUTION", "MODULE"):
        # Not a ShapePipe pipeline config (e.g. config_MCCD.ini is the MCCD
        # library's own params file, read by mccd_preprocessing_runner/
        # mccd_fit_val_runner via CONFIG_PATH, not by ShapePipe's [EXECUTION]).
        return []

    problems = []
    for module, section in _module_config_sections(parser.get("EXECUTION", "MODULE")):
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

    return problems


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

    problems = _missing_required_keys(parser)

    assert not problems, f"{config_path}: " + "; ".join(problems)


def test_module_sections_check_catches_a_missing_repeated_section():
    """A module used twice must have both of its sections checked.

    Regression guard for _module_config_sections: it must return one
    (module, section) pair per invocation of a repeated module, not
    collapse them, or deleting one of the two sections a module needs
    (here vignetmaker_runner's _RUN_1) would go unnoticed.
    """
    config_path = next(
        p for p in WORKFLOW_CONFIG_FILES if p.name == "config_tile_PiViVi_mccd.ini"
    )

    parser = configparser.ConfigParser()
    parser.read(config_path)
    assert parser.has_section("VIGNETMAKER_RUNNER_RUN_1")

    parser.remove_section("VIGNETMAKER_RUNNER_RUN_1")

    problems = _missing_required_keys(parser)

    assert any("VIGNETMAKER_RUNNER_RUN_1" in problem for problem in problems)
