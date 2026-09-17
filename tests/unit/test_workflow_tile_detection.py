"""The two tile_detect modes agree on where the sexcat is.

``tile_detection: unions_catalogue`` swaps the SExtractor rule for a fetch plus
a conversion, and everything downstream of ``tile_detect`` is unchanged only
because both modes write ``run_sp_tile_Sx`` and both are checked as the
``tile_detect`` stage. The agreements that make that true are between files
that never see each other at run time: the two inis' RUN_NAMEs,
``completeness.STAGE_DIR``, the flavoured ``COMPLETENESS['tile_detect']``
table, and ``run_report``'s stage list. They are asserted here, statically.
Container-free: the scripts are stdlib-only.
"""

import configparser
import importlib.util
import re
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
    assert set(completeness.COMPLETENESS["tile_detect"]) == set(
        completeness.TILE_DETECTIONS)


def test_both_detection_inis_write_the_stage_dir():
    """config_tile_Sx and config_tile_Uc share the run dir STAGE_DIR names."""
    level, subdir = completeness.STAGE_DIR["tile_detect"]
    assert level == "tile"
    for name in ("config_tile_Sx.ini", "config_tile_Uc.ini"):
        assert _ini(name)["DEFAULT"]["RUN_NAME"].strip() == subdir, (
            f"{name} writes a run dir other than {subdir}: unit_pre would clear "
            "the wrong directory and the chain downstream would read nothing.")
    assert _ini("config_tile_Uc.ini")["DEFAULT"]["RUN_DATETIME"].strip() == "False"


def test_fetch_ini_matches_its_stage_and_feeds_the_converter():
    level, gic = completeness.STAGE_DIR["tile_get_catalogue"]
    assert level == "tile"
    assert _ini("config_tile_Gic.ini")["DEFAULT"]["RUN_NAME"].strip() == gic
    uc = _ini("config_tile_Uc.ini")["READ_EXT_SEXCAT_RUNNER"]
    assert f"$SP_RUN/output/{gic}/get_images_runner/output" in uc["INPUT_DIR"]
    assert uc["FILE_PATTERN"].split(",")[0].strip() == "CFIS_cat"
    # The multi-epoch post-processing is what gives the sexcat its EPOCH_k
    # extensions; ngmix_range.py refuses a sexcat without them.
    assert uc["MAKE_POST_PROCESS"].strip() == "True"


def _stage_dir(tmp_path, runner, n):
    out = tmp_path / runner / "output"
    out.mkdir(parents=True)
    for k in range(n):
        (out / f"f{k}").write_text("")
    return tmp_path


@pytest.mark.parametrize("mode, runner, expect", [
    ("sextractor", "sextractor_runner", 2),
    ("unions_catalogue", "read_ext_sexcat_runner", 1),
])
def test_tile_detect_is_checked_per_mode(tmp_path, monkeypatch, mode, runner, expect):
    monkeypatch.setenv("SP_TILE_DETECTION", mode)
    ok, details = completeness.check_counts(
        "tile_detect", _stage_dir(tmp_path, runner, expect))
    assert ok and details == [(runner, expect, expect, False)]


def test_unset_mode_is_sextractor(tmp_path, monkeypatch):
    """A prologue without the export checks as before: data runs are unchanged."""
    monkeypatch.delenv("SP_TILE_DETECTION", raising=False)
    ok, details = completeness.check_counts(
        "tile_detect", _stage_dir(tmp_path, "sextractor_runner", 2))
    assert ok and details[0][0] == "sextractor_runner"


def test_invalid_mode_is_fatal(tmp_path, monkeypatch):
    monkeypatch.setenv("SP_TILE_DETECTION", "steven")
    with pytest.raises(ValueError, match="SP_TILE_DETECTION"):
        completeness.check_counts("tile_detect", tmp_path)


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


# --- the segmentation half: unions_catalogue + uberseg ----------------------


def test_blend_handlings_mirror_the_ngmix_module():
    """completeness.BLEND_HANDLINGS is a copy; the copy must stay true.

    The Snakefile validates `blend_handling:` against it in the launcher venv,
    outside the container where shapepipe is importable, so the tuple is
    mirrored rather than imported. Read out of the source text here for the
    same reason: this file is container-free.
    """
    src = (REPO_ROOT / "src" / "shapepipe" / "modules" / "ngmix_package"
           / "ngmix.py").read_text()
    match = re.search(r"^BLEND_HANDLINGS = \(([^)]*)\)", src, re.M)
    assert match, "ngmix.py no longer defines BLEND_HANDLINGS"
    assert completeness.BLEND_HANDLINGS == tuple(
        part.strip().strip('"\'') for part in match.group(1).split(",")
        if part.strip())


def test_segmentation_ini_matches_its_stage_and_yields_the_map():
    level, sg = completeness.STAGE_DIR["tile_segmentation"]
    assert level == "tile"
    ini = _ini("config_tile_Sg.ini")
    assert ini["DEFAULT"]["RUN_NAME"].strip() == sg
    assert ini["DEFAULT"]["RUN_DATETIME"].strip() == "False"
    sx = ini["SEXTRACTOR_RUNNER"]
    # The check image is the whole point, and it is the only one asked for.
    assert [c.strip() for c in sx["CHECKIMAGE"].split(",")] == ["SEGMENTATION"]
    # No multi-epoch post-processing: those extensions ride on the sexcat
    # config_tile_Uc.ini writes, and without it the WCS sqlite is not an input.
    assert sx["MAKE_POST_PROCESS"].strip() == "False"
    assert "log_exp_headers" not in sx["FILE_PATTERN"]
    # Same detection settings as the SExtractor mode, so the footprints are
    # the ones that mode would have drawn.
    detect = _ini("config_tile_Sx.ini")["SEXTRACTOR_RUNNER"]
    for key in ("DOT_SEX_FILE", "DOT_PARAM_FILE", "DOT_CONV_FILE",
                "WEIGHT_IMAGE", "FLAG_IMAGE"):
        assert sx[key].strip() == detect[key].strip(), key


def test_segmentation_stage_is_checked(tmp_path):
    """Two products: the SEGMENTATION check image and SExtractor's own sexcat."""
    ok, details = completeness.check_counts(
        "tile_segmentation", _stage_dir(tmp_path, "sextractor_runner", 2))
    assert ok and details == [("sextractor_runner", 2, 2, False)]


@pytest.mark.parametrize("detection, blend, present", [
    ("unions_catalogue", "uberseg", True),
    ("unions_catalogue", "noisefill", False),
    ("sextractor", "uberseg", False),
])
def test_report_lists_the_segmentation_stage_only_for_the_pair(
        monkeypatch, detection, blend, present):
    monkeypatch.setenv("SP_TILE_DETECTION", detection)
    monkeypatch.setenv("SP_BLEND_HANDLING", blend)
    stages = _load("run_report").TILE_STAGES
    assert ("tile_segmentation" in stages) is present
    if present:
        assert stages.index("tile_segmentation") == stages.index(
            "tile_detect") + 1


def test_run_config_defaults_to_the_catalogue_and_declares_its_source():
    """The committed run config must parse under its own defaults.

    `tile_detection: unions_catalogue` is refused by the Snakefile without
    `inputs.catalogues`, so the two settings travel together.
    """
    text = (REPO_ROOT / "workflow" / "config.yaml").read_text()
    config = {}
    for line in text.splitlines():
        if line.startswith("tile_detection:") or line.startswith(
                "blend_handling:"):
            key, _, value = line.partition(":")
            config[key] = value.strip()
    assert config["tile_detection"] == "unions_catalogue"
    assert config["blend_handling"] in completeness.BLEND_HANDLINGS
    assert re.search(r"^  catalogues: \S+", text, re.M)
