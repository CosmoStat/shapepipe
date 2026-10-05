"""The run config's ``defect_weighting`` reaches ngmix's DEFECT_WEIGHTING.

The workflow reads ``defect_weighting`` (Snakefile), tile_ngmix exports it as
``SP_DEFECT_WEIGHTING``, and config_tile_Ng_template.ini expands it into
``[NGMIX_RUNNER] DEFECT_WEIGHTING``; empty takes the module default. These
files never see each other at run time, so their agreement is asserted here,
statically, plus the expansion through ShapePipe's own config parser.
"""

import ast
import re
from pathlib import Path

import pytest

from shapepipe.modules.ngmix_package.ngmix import (
    DEFECT_WEIGHTING,
    DEFECT_WEIGHTINGS,
)
from shapepipe.pipeline.config import CustomParser

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / "workflow"
TEMPLATE = WORKFLOW / "config" / "cfis" / "config_tile_Ng_template.ini"


def test_the_default_is_an_option():
    assert DEFECT_WEIGHTING in DEFECT_WEIGHTINGS
    assert set(DEFECT_WEIGHTINGS) == {"des_y6", "fourfold_zero", "hole",
                                      "full"}


def test_the_snakefile_reads_the_options_from_the_module_it_runs():
    """The Snakefile validates ``defect_weighting`` against the tuple it
    parses out of ngmix.py; parse it the same way here."""
    snakefile = (WORKFLOW / "Snakefile").read_text()
    assert 'config.get("defect_weighting")' in snakefile
    assert "DEFECT_WEIGHTINGS" in snakefile
    tree = ast.parse(
        (REPO_ROOT / "src" / "shapepipe" / "modules" / "ngmix_package"
         / "ngmix.py").read_text()
    )
    parsed = next(
        ast.literal_eval(node.value) for node in tree.body
        if isinstance(node, ast.Assign)
        and [getattr(t, "id", None) for t in node.targets]
        == ["DEFECT_WEIGHTINGS"]
    )
    assert parsed == DEFECT_WEIGHTINGS


def test_tile_ngmix_exports_it_and_tile_vignets_carries_it():
    """tile_ngmix exports SP_DEFECT_WEIGHTING; tile_vignets carries the
    value as a param, so a change reruns the whole group, never the chunks
    alone."""
    rules = (WORKFLOW / "rules" / "tile.smk").read_text()
    ngmix = rules[rules.index("rule tile_ngmix:"):]
    vignets = rules[rules.index("rule tile_vignets:"):
                    rules.index("rule tile_ngmix:")]
    assert re.search(r'"SP_DEFECT_WEIGHTING":\s*DEFECT_WEIGHTING', ngmix)
    assert re.search(r"defect_weighting\s*=\s*DEFECT_WEIGHTING", vignets)


@pytest.mark.parametrize("exported", ["", "des_y6"])
def test_the_template_expands_the_export(monkeypatch, exported):
    parser = CustomParser()
    parser.read(TEMPLATE)
    monkeypatch.setenv("SP_DEFECT_WEIGHTING", exported)
    assert parser.getexpanded("NGMIX_RUNNER", "DEFECT_WEIGHTING") == exported
    monkeypatch.delenv("SP_DEFECT_WEIGHTING")
    assert parser.getexpanded("NGMIX_RUNNER", "DEFECT_WEIGHTING") == ""
