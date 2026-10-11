"""The run config's ngmix defect options reach the ngmix runner.

``defect_fill`` and ``defect_weighting`` each travel the same path: the
workflow reads the key (Snakefile) and validates it against the options
tuple in ngmix.py, tile_ngmix exports it as ``SP_<KEY>``, and
config_tile_Ng_template.ini expands that into ``[NGMIX_RUNNER] <KEY>``;
empty takes the module default. These files never see each other at run
time, so their agreement is asserted here, statically, plus the expansion
through ShapePipe's own config parser.
"""

import ast
import re
from pathlib import Path

import pytest

from shapepipe.modules.ngmix_package import ngmix
from shapepipe.pipeline.config import CustomParser

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / "workflow"
TEMPLATE = WORKFLOW / "config" / "cfis" / "config_tile_Ng_template.ini"

# (run-config key, the options tuple, a non-default option)
OPTIONS = [
    ("defect_fill", {"interpolate", "noise"}, "noise"),
    ("defect_weighting", {"des_y6", "fourfold_zero", "hole", "full"},
     "des_y6"),
]


def _names(key):
    upper = key.upper()
    return upper, upper + "S", "SP_" + upper


@pytest.mark.parametrize("key, options, _", OPTIONS)
def test_the_default_is_an_option(key, options, _):
    default, tuple_name, _env = _names(key)
    assert getattr(ngmix, default) in getattr(ngmix, tuple_name)
    assert set(getattr(ngmix, tuple_name)) == options


@pytest.mark.parametrize("key, _, __", OPTIONS)
def test_the_snakefile_reads_the_options_from_the_module_it_runs(key, _, __):
    """The Snakefile validates the key against the tuple it parses out of
    ngmix.py; parse it the same way here."""
    default, tuple_name, _env = _names(key)
    snakefile = (WORKFLOW / "Snakefile").read_text()
    assert f'config.get("{key}")' in snakefile
    assert f'_ngmix_options("{tuple_name}")' in snakefile
    assert re.search(rf'\("{key}", {default}, {tuple_name}\)', snakefile)
    tree = ast.parse(
        (REPO_ROOT / "src" / "shapepipe" / "modules" / "ngmix_package"
         / "ngmix.py").read_text()
    )
    parsed = next(
        ast.literal_eval(node.value) for node in tree.body
        if isinstance(node, ast.Assign)
        and [getattr(t, "id", None) for t in node.targets] == [tuple_name]
    )
    assert parsed == getattr(ngmix, tuple_name)


@pytest.mark.parametrize("key, _, __", OPTIONS)
def test_tile_ngmix_exports_it_and_tile_vignets_carries_it(key, _, __):
    """tile_ngmix exports SP_<KEY>; tile_vignets carries the value as a
    param, so a change reruns the whole group, never the chunks alone."""
    default, _tuple, env = _names(key)
    rules = (WORKFLOW / "rules" / "tile.smk").read_text()
    ngmix_rule = rules[rules.index("rule tile_ngmix:"):]
    vignets = rules[rules.index("rule tile_vignets:"):
                    rules.index("rule tile_ngmix:")]
    assert re.search(rf'"{env}":\s*{default}\b', ngmix_rule)
    assert re.search(rf"{key}\s*=\s*{default}\b", vignets)


@pytest.mark.parametrize("key, _, exported", OPTIONS)
def test_the_template_expands_the_export(monkeypatch, key, _, exported):
    option, _tuple, env = _names(key)
    parser = CustomParser()
    parser.read(TEMPLATE)
    for value in ("", exported):
        monkeypatch.setenv(env, value)
        assert parser.getexpanded("NGMIX_RUNNER", option) == value
    monkeypatch.delenv(env)
    assert parser.getexpanded("NGMIX_RUNNER", option) == ""
