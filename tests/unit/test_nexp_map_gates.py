"""Fitted-PSF gates on footprints, reclamation and the exposure-count map.

Failure modes: fake PSFs request nonexistent persistence products; a fitted
PSF loses a durable edge before cleanup; a fake-PSF (sim) campaign requests
an exposure-count map it cannot build, or a fitted one silently loses it. Lift the
Snakefile's functions and the rule's input lambdas into sentinel roots, as
test_campaign_lineage does, without parsing a second workflow DAG.
"""

import re
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

WORKFLOW = Path(__file__).resolve().parents[2] / "workflow"
SNAKEFILE = WORKFLOW / "Snakefile"
PRODUCTS = "/sentinel/products"
EXPOSURES = ["2605805", "2700001"]


class WorkflowError(Exception):
    """Stand-in for snakemake.exceptions.WorkflowError, which the image lacks."""


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
        "glob": __import__("glob"),
        "logger": SimpleNamespace(warning=lambda message: None),
    }
    text = SNAKEFILE.read_text()
    assignment = re.search(r"^PERSISTS_PSF\s*=.*$", text, re.M)
    assert assignment, "Snakefile must bind PERSISTS_PSF"
    exec(assignment.group(0), ns)
    # Fake PSFs are only legal with image sims, which rasterize no defects.
    ns["INPUT_TYPE"] = "image_sims" if psf_model == "fake" else "data"
    defects = re.search(r"^MAPS_DEFECTS\s*=.*$", text, re.M)
    assert defects, "Snakefile must bind MAPS_DEFECTS"
    exec(defects.group(0), ns)
    for name in ("prod_exp_dir", "prod_exp_manifest", "exp_dir", "tombstone",
                 "exp_manifest", "exp_store_reclaimed", "footprint_edge",
                 "tile_dir", "tile_manifest", "psf_exposures",
                 "footprint_targets", "flag", "prod_exp_tar",
                 "prod_exp_fragment", "nexp_map", "nexp_map_manifest",
                 "nexp_map_exposures", "nexp_map_targets"):
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


def _nexp_gate(ns, nexp):
    """Execute the Snakefile's exposure-count-map gate against a `nexp:` block."""
    text = SNAKEFILE.read_text()
    gate = re.search(r"^_NEXP\s*=.*\n^MAPS_NEXP\s*=.*\n", text, re.M)
    assert gate, "Snakefile must gate exposure_maps.nexp on a fitted PSF"
    ns["_MAPS"] = {} if nexp is None else {"nexp": nexp}
    exec(gate.group(0), ns)


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
        stages = ["exp_persist"]
        if ns["MAPS_DEFECTS"]:
            stages.append("exp_defect_map")
        stages.append("exp_footprint")
        expected += [
            f"{PRODUCTS}/exp/26/2605805/manifests/{stage}.json"
            for stage in stages
        ]
    assert edges == expected


@pytest.mark.parametrize("psf_model", ["fake", "psfex"])
@pytest.mark.parametrize("nexp", [None, {}, {"enabled": True},
                                  {"enabled": False}, {"enabled": "true"},
                                  {"enabled": "false"}])
def test_nexp_map_follows_the_psf_model(psf_model, nexp, tmp_path):
    """On by default for a fitted PSF, off only by explicit opt-out; a fake-PSF
    (sim) campaign skips it silently, whatever the block says."""
    ns = _namespace(psf_model, tmp_path)
    _nexp_gate(ns, nexp)
    opted_out = (nexp or {}).get("enabled") in (False, "false")
    built = psf_model == "psfex" and not opted_out
    expected = [ns["nexp_map"](), ns["nexp_map_manifest"]()] if built else []
    assert ns["nexp_map_targets"]() == expected
