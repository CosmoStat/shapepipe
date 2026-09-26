"""Fitted-PSF gates on footprints, reclamation and campaign coverage.

Failure modes: fake PSFs request nonexistent persistence products; a fitted
PSF loses either durable edge before cleanup; coverage with fake PSFs reaches
submission instead of failing at parse time. Lift the Snakefile's functions
and the rule's input lambdas into sentinel roots, as test_campaign_lineage
does, without parsing a second workflow DAG.
"""

import re
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest
from snakemake.exceptions import WorkflowError

WORKFLOW = Path(__file__).resolve().parents[2] / "workflow"
SNAKEFILE = WORKFLOW / "Snakefile"
PRODUCTS = "/sentinel/products"
EXPOSURES = ["2605805", "2700001"]


def _namespace(psf_model, tmp_path, *, main=True):
    """Bind the real PSF gate and target helpers to a small indexed campaign."""
    tile_exp = {
        "210.282": [EXPOSURES[1], EXPOSURES[0]],
        "211.282": [EXPOSURES[0]],
        "999.999": ["9900001"],
    }
    ns = {
        "PSF_MODEL": psf_model,
        "PRODUCTS_DIR": Path(PRODUCTS),
        "RUN_DIR": tmp_path / "scratch",
        "CAMPAIGN": "campaign-sentinel",
        "TILES_READY": ["210.282", "211.282"],
        "READY_SET": {"210.282", "211.282"},
        "tile_exposures": tile_exp.__getitem__,
        "clean_consumers": lambda exp: ["210.282", "999.999"],
        "workflow": SimpleNamespace(is_main_process=main),
        "WorkflowError": WorkflowError,
        "Path": Path,
        "script_hash": lambda name: "hash-sentinel",
    }
    text = SNAKEFILE.read_text()
    assignment = re.search(r"^PERSISTS_PSF\s*=.*$", text, re.M)
    assert assignment, "Snakefile must bind PERSISTS_PSF"
    exec(assignment.group(0), ns)
    for name in ("prod_exp_dir", "prod_exp_manifest", "exp_dir", "tombstone",
                 "tile_dir", "tile_manifest", "psf_exposures",
                 "footprint_targets", "flag"):
        definition = re.search(
            rf"^def {name}\(.*?(?=^\S|\Z)", text, re.M | re.S)
        assert definition, f"Snakefile must define {name}()"
        exec(definition.group(0), ns)
    return ns


def _clean_inputs(ns):
    """Evaluate clean_exposure's actual input functions, not a copy of them."""
    text = (WORKFLOW / "rules" / "exposure.smk").read_text()
    rule = re.search(r"^rule clean_exposure:\n(.*?)(?=^rule |\Z)",
                     text, re.M | re.S)
    assert rule, "exposure.smk must define clean_exposure"
    inputs = re.search(r"^    input:\n(.*?)(?=^    \w+:)",
                       rule.group(1), re.M | re.S)
    assert inputs, "clean_exposure must declare its ordering edges"
    return eval("[\n" + textwrap.dedent(inputs.group(1)) + "\n]", ns)


def _coverage_preamble(ns, enabled):
    """Execute the coverage constants, guard and helper before its rule."""
    text = (WORKFLOW / "rules" / "coverage.smk").read_text()
    preamble, separator, _ = text.partition("\nrule coverage_map:")
    assert separator, "coverage.smk must define coverage_map"
    ns["config"] = {"coverage": {"enabled": enabled}}
    exec(preamble, ns)


@pytest.mark.parametrize("psf_model", ["fake", "psfex"])
@pytest.mark.parametrize("main", [True, False])
def test_footprint_targets_require_fitted_psf(psf_model, main, tmp_path):
    """Only a head-process fitted-PSF run requests the ready exposures."""
    ns = _namespace(psf_model, tmp_path, main=main)
    expected = []
    if psf_model == "psfex" and main:
        expected = [
            f"{PRODUCTS}/exp/{exp[:2]}/{exp}/manifests/exp_footprint.json"
            for exp in EXPOSURES
        ]
    assert ns["footprint_targets"]() == expected


@pytest.mark.parametrize("psf_model", ["fake", "psfex"])
def test_clean_exposure_waits_on_persist_iff_psf(psf_model, tmp_path):
    """Reclamation waits for both durable reads iff a PSF is fitted."""
    ns = _namespace(psf_model, tmp_path)
    wc = SimpleNamespace(exp=EXPOSURES[0])
    edges = [path for function in _clean_inputs(ns) for path in function(wc)]
    expected = [str(tmp_path / "scratch" / "tiles" / "21" / "210.282"
                    / "manifests" / "tile_vignets.json")]
    if psf_model == "psfex":
        expected += [
            f"{PRODUCTS}/exp/26/2605805/manifests/{stage}.json"
            for stage in ("exp_persist", "exp_footprint")
        ]
    assert edges == expected


@pytest.mark.parametrize("psf_model", ["fake", "psfex"])
@pytest.mark.parametrize("enabled", [True, False, "true", "false"])
def test_coverage_refuses_enabled_fake_at_parse_time(
        psf_model, enabled, tmp_path):
    """Fake PSFs are refused only when coverage is enabled, before any rule."""
    ns = _namespace(psf_model, tmp_path)
    is_enabled = enabled in (True, "true")
    if psf_model == "fake" and is_enabled:
        with pytest.raises(WorkflowError,
                           match=r"coverage.enabled.*psf_model=fake"):
            _coverage_preamble(ns, enabled)
    else:
        _coverage_preamble(ns, enabled)
        expected = ([ns["COVERAGE_HSP"], ns["COVERAGE_MANIFEST"]]
                    if is_enabled else [])
        assert ns["coverage_targets"]() == expected
