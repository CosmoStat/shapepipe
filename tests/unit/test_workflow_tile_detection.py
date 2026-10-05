"""The two tile_detection modes run one tile_detect, and agree on its inputs.

Both modes run SExtractor (config_tile_Sx.ini) in the ``run_sp_tile_Sx`` stage
dir; ``tile_detection: unions_catalogue`` adds ``tile_get_catalogue`` and
points the ini's MATCH_CATALOGUE at its output through the workflow's
SP_MATCH_CATALOGUE. The agreements that make that work are between files that
never see each other at run time: the inis' RUN_NAMEs and fetch patterns,
``completeness``'s tables, ``run_report``'s stage list and tile.smk's path.
They are asserted here, statically. Container-free: the scripts are
stdlib-only.
"""

import configparser
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
CONFIG_DIR = REPO_ROOT / "workflow" / "config" / "cfis"


def _load(name):
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}", SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


def _ini(name):
    parser = configparser.ConfigParser()
    assert parser.read(CONFIG_DIR / name) == [str(CONFIG_DIR / name)]
    return parser


completeness = _load("completeness")


def test_modes_are_the_two_the_rules_know():
    assert completeness.TILE_DETECTIONS == ("sextractor", "unions_catalogue")


def test_detection_ini_writes_the_stage_dir():
    level, subdir = completeness.STAGE_DIR["tile_detect"]
    assert level == "tile"
    assert _ini("config_tile_Sx.ini")["DEFAULT"]["RUN_NAME"].strip() == subdir


def test_detection_joins_only_when_the_workflow_says_so():
    """MATCH_CATALOGUE is empty unless SP_MATCH_CATALOGUE is exported."""
    sx = _ini("config_tile_Sx.ini")["SEXTRACTOR_RUNNER"]
    assert sx["MATCH_CATALOGUE"].strip() == "${SP_MATCH_CATALOGUE:-}"
    assert float(sx["MATCH_RADIUS"]) == 1.0
    assert 0.9 < float(sx["MATCH_MIN_FRACTION"]) <= 1
    # The join renumbers before the post-processing keys epochs on NUMBER.
    assert sx["MAKE_POST_PROCESS"].strip() == "True"


def test_fetch_ini_writes_the_catalogue_tile_smk_joins():
    """Gic's output is the path tile.smk exports as SP_MATCH_CATALOGUE."""
    level, gic = completeness.STAGE_DIR["tile_get_catalogue"]
    assert level == "tile"
    assert _ini("config_tile_Gic.ini")["DEFAULT"]["RUN_NAME"].strip() == gic
    gi = _ini("config_tile_Gic.ini")["GET_IMAGES_RUNNER"]
    assert gi["OUTPUT_FILE_PATTERN"].strip() == "CFIS_cat-"
    assert gi["INPUT_FILE_EXT"].strip() == ".cat"
    smk = (REPO_ROOT / "workflow" / "rules" / "tile.smk").read_text()
    assert f'/output/{gic}/get_images_runner/output"' in smk
    assert 'f"{gic}/CFIS_cat{unit_num(tile)}.cat"' in smk


def _stage_dir(tmp_path, runner, n):
    out = tmp_path / runner / "output"
    out.mkdir(parents=True)
    for k in range(n):
        (out / f"f{k}").write_text("")
    return tmp_path


@pytest.mark.parametrize("stage, runner, expect", [
    ("tile_detect", "sextractor_runner", 2),
    ("tile_get_catalogue", "get_images_runner", 1),
])
def test_stage_counts(tmp_path, stage, runner, expect):
    ok, details = completeness.check_counts(
        stage, _stage_dir(tmp_path, runner, expect))
    assert ok and details == [(runner, expect, expect, False)]


@pytest.mark.parametrize("mode, present", [
    ("sextractor", False), ("unions_catalogue", True), (None, False),
])
def test_report_lists_the_fetch_stage_only_for_catalogue_runs(monkeypatch, mode, present):
    if mode is None:
        monkeypatch.delenv("SP_TILE_DETECTION", raising=False)
    else:
        monkeypatch.setenv("SP_TILE_DETECTION", mode)
    stages = _load("run_report").TILE_STAGES
    assert ("tile_get_catalogue" in stages) is present
    if present:
        assert stages.index("tile_get_catalogue") == stages.index("tile_detect") - 1


@pytest.mark.parametrize("machine", ["nibi", "candide"])
@pytest.mark.parametrize("input_type", ["data", "image_sims"])
def test_defaults_pair_the_catalogue_with_its_source(
        machine, input_type, monkeypatch, tmp_path):
    """Real data defaults to the catalogue, and every machine says where it is.

    tile_detection follows from the input type (the input_types: table), but
    `tile_detection: unions_catalogue` is refused by the Snakefile without
    `inputs.catalogues`, which is a per-machine path: every machine with a
    data entry has to declare one. Image sims have no UNIONS catalogue and
    use SExtractor.
    """
    pytest.importorskip("yaml")
    run_config = _load("run_config")
    monkeypatch.setenv("SP_PROFILE", machine)
    over = tmp_path / "run.yaml"
    over.write_text(f"input_type: {input_type}\n")
    config = run_config.load(
        str(REPO_ROOT / "workflow" / "config.yaml"), str(over))
    if input_type == "data":
        assert config["tile_detection"] == "unions_catalogue"
        assert run_config.catalogue_source(config)[0]
    else:
        assert config["tile_detection"] == "sextractor"


@pytest.mark.parametrize("catalogues", [None, "", "TBD"])
def test_catalogue_run_without_a_source_fails_at_parse(catalogues):
    """An unset or placeholder `inputs.catalogues` is refused before any job runs."""
    run_config = _load("run_config")
    config = {"tile_detection": "unions_catalogue",
              "inputs": {} if catalogues is None else {"catalogues": catalogues}}
    with pytest.raises(ValueError, match="inputs.catalogues"):
        run_config.catalogue_source(config)
    config["tile_detection"] = "sextractor"
    assert run_config.catalogue_source(config) == ("", "symlink")


@pytest.mark.parametrize("source, retrieve", [
    ("vos:cfis/tiles_DR6", "vos"), ("/data/tiles_DR6", "symlink"),
])
def test_catalogue_retrieve_follows_the_prefix(source, retrieve):
    run_config = _load("run_config")
    config = {"tile_detection": "unions_catalogue",
              "inputs": {"catalogues": source}}
    assert run_config.catalogue_source(config) == (source, retrieve)
