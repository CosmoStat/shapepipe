"""Check recorded values, not just locations, without importing the pipeline.

Failure modes: plausible numeric drift, wrong key/section/case, commented or
ambiguous settings, changed cut operators, lost list elements, and executing
Python while trying to inspect it. Fixtures exercise each through the same
reader as the record; location resolution must remain independent of equality.
"""

from pathlib import Path

import pytest

from tests.helpers.astra_record import (
    check_anchor_value,
    extract_anchors,
    load_yaml,
    resolve_anchor,
    value_errors,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def record_with(reference):
    """Put a reference in a scoped decision, as in the real record."""

    return {
        "analyses": {
            "detection": {
                "decisions": {
                    "threshold": {
                        "rationale": f"Affects selection. Anchor: {reference}."
                    }
                }
            }
        }
    }


def test_every_astra_value_matches():
    record = load_yaml(REPO_ROOT / "astra.yaml")
    anchors = extract_anchors(record)
    assert any(" = " in ref for a in anchors for ref in a.references), (
        "astra.yaml contains no value assertions"
    )
    errors = value_errors(REPO_ROOT, record)
    assert not errors, "ASTRA value mismatches:\n - " + "\n - ".join(errors)


def test_anchor_grammar_keeps_values_and_location_only_refs(tmp_path):
    (tmp_path / "image.sex").write_text("THRESH 1.25\n", encoding="utf-8")
    (tmp_path / "code.py").write_text("def fit():\n    pass\n", encoding="utf-8")
    record = record_with("image.sex#THRESH = 1.25; code.py::fit; image.sex")
    anchor, = extract_anchors(record)
    assert anchor.error is None
    assert anchor.references == (
        "image.sex#THRESH = 1.25", "code.py::fit", "image.sex"
    )
    assert all(resolve_anchor(tmp_path, ref) is None for ref in anchor.references)
    assert value_errors(tmp_path, record) == []


@pytest.mark.parametrize(
    "actual, expected",
    [
        ("1.0", "1"),
        ("13.", "13.0"),
        ("5e-4", "0.0005"),
        ("-1e3", "-1000"),
        ("yes", "on"),
        ("off", "No"),
        ("1", "True"),  # ConfigParser.getboolean reads 1/0 as booleans.
        ("0", "false"),
        ("2.5, 3.5", "[2.50, 3.500]"),
        ("XWIN_IMAGE,YWIN_IMAGE", "XWIN_IMAGE, YWIN_IMAGE"),
        ("$ROOT/data{number}.fits", '"$ROOT/data{number}.fits"'),
    ],
)
def test_equivalent_spellings(tmp_path, actual, expected):
    (tmp_path / "config.ini").write_text(f"[SCIENCE]\nKEY = {actual}\n")
    ref = f"config.ini#SCIENCE.KEY = {expected}"
    assert check_anchor_value(tmp_path, ref) is None


@pytest.mark.parametrize("filename", ["config.sex", "model.psfex"])
@pytest.mark.parametrize(
    "actual, expected", [("Y", "True"), ("false", "N"), ("n", "off")]
)
def test_astromatic_formats_accept_y_n(tmp_path, filename, actual, expected):
    (tmp_path / filename).write_text(f"KEY {actual}\n")
    assert check_anchor_value(tmp_path, f"{filename}#KEY = {expected}") is None


@pytest.mark.parametrize(
    "filename, contents, selector, expected",
    [
        # getboolean raises on Y/N, so the record must not call them True/False.
        ("config.ini", "[S]\nKEY = Y\n", "S.KEY", "True"),
        ("config.ini", "[S]\nKEY = True\n", "S.KEY", "Y"),
        ("config.ini", "[S]\nKEY = 2\n", "S.KEY", "True"),
        ("stars.setools", "[RAND_SPLIT:s]\nKEY = Y\n", "RAND_SPLIT:s.KEY", "True"),
        ("constants.py", "KEY = True\n", "KEY", "Y"),
    ],
)
def test_boolean_words_follow_each_reader(
    tmp_path, filename, contents, selector, expected
):
    (tmp_path / filename).write_text(contents)
    separator = "::" if filename.endswith(".py") else "#"
    ref = f"{filename}{separator}{selector} = {expected}"
    assert check_anchor_value(tmp_path, ref) is not None


@pytest.mark.parametrize(
    "actual, expected",
    [
        ("1.000001", "1"),
        ("1", "True"),
        ("0", "False"),
        ("51,51", "51"),  # No scalar broadcasting for ordinary keys.
        ("51,52", "51,51"),
        ("1,2", "2,1"),
        ("1,1,1", "1,1"),
        ("1", "[1]"),
        ("map_weight", "MAP_WEIGHT"),
    ],
)
def test_normalization_does_not_hide_drift(tmp_path, actual, expected):
    (tmp_path / "config.sex").write_text(f"KEY {actual}\n")
    problem = check_anchor_value(tmp_path, f"config.sex#KEY = {expected}")
    assert "expected" in problem and "actual" in problem


@pytest.mark.parametrize(
    "filename, line, selector",
    [
        ("default.sex", "THRESH 0.0005 # comment", "THRESH"),
        ("default.psfex", "PSF_ACCURACY 0.0005", "PSF_ACCURACY"),
        ("default.ww", "WEIGHT_MIN = 0.0005", "WEIGHT_MIN"),
        ("default.conf", "THRESH 0.0005", "THRESH"),
        ("stars.setools", "[RAND_SPLIT:stars]\nRATIO = 0.0005",
         "RAND_SPLIT:stars.RATIO"),
    ],
)
def test_line_config_formats(tmp_path, filename, line, selector):
    (tmp_path / filename).write_text(line + "\n")
    ref = f"{filename}#{selector} = 5e-4"
    assert resolve_anchor(tmp_path, ref) is None
    assert check_anchor_value(tmp_path, ref) is None


@pytest.mark.parametrize(
    "filename, line, key",
    [
        ("columns.param", "VIGNET(51,51)", "VIGNET"),
        ("model.psfex", "PSF_SIZE 51,51", "PSF_SIZE"),
        ("model.psfex", "PSF_SIZE 51", "PSF_SIZE"),
    ],
)
def test_square_stamp_shorthand(tmp_path, filename, line, key):
    target = tmp_path / filename
    target.write_text(line + " # stamp\n")
    for expected in ("51", "51,51", "[51.0, 51]"):
        ref = f"{filename}#{key} = {expected}"
        assert check_anchor_value(tmp_path, ref) is None
    target.write_text(line.replace("51", "53", 1))
    assert check_anchor_value(tmp_path, f"{filename}#{key} = 51") is not None


def test_ini_case_sections_defaults_and_interpolation(tmp_path):
    (tmp_path / "config.ini").write_text(
        "[DEFAULT]\nENABLED = True\n"
        "[SCIENCE]\nKey = 30\nKEY = 31\nPATH = $DATA/%s/file\n"
        "[OTHER]\nKEY = 99\n"
    )
    for selector, value in (
        ("SCIENCE.Key", "30"), ("SCIENCE.KEY", "31"),
        ("OTHER.KEY", "99"), ("SCIENCE.ENABLED", "yes"),
        ("DEFAULT.ENABLED", "True"), ("SCIENCE.PATH", "$DATA/%s/file"),
    ):
        ref = f"config.ini#{selector} = {value}"
        assert resolve_anchor(tmp_path, ref) is None
        assert check_anchor_value(tmp_path, ref) is None
    assert check_anchor_value(tmp_path, "config.ini#SCIENCE.key = 31") is not None
    assert check_anchor_value(tmp_path, "config.ini#KEY = 31") is not None


@pytest.mark.parametrize(
    "filename, contents, selector",
    [
        ("config.sex", "# KEY 1\nKEY_EXTRA 1\n", "KEY"),
        ("config.param", "# VIGNET(51,51)\n", "VIGNET"),
        ("config.setools", "[MASK:stars]\n# FLAGS == 0\n", "MASK:stars.FLAGS"),
    ],
)
def test_commented_keys_resolve_but_cannot_assert_active_values(
    tmp_path, filename, contents, selector
):
    (tmp_path / filename).write_text(contents)
    assert resolve_anchor(tmp_path, f"{filename}#{selector}") is None
    problem = check_anchor_value(tmp_path, f"{filename}#{selector} = 1")
    assert "active" in problem


@pytest.mark.parametrize(
    "filename, contents, selector",
    [
        ("config.sex", "KEY 1\nKEY 2\n", "KEY"),
        ("config.ini", "[SCIENCE]\nKEY = 1\nKEY = 2\n", "SCIENCE.KEY"),
        ("config.setools", "[RAND_SPLIT:s]\nRATIO = 1\nRATIO = 2\n",
         "RAND_SPLIT:s.RATIO"),
    ],
)
def test_duplicate_settings_fail_closed(tmp_path, filename, contents, selector):
    (tmp_path / filename).write_text(contents)
    assert check_anchor_value(tmp_path, f"{filename}#{selector} = 2") is not None


def test_setools_cuts_keep_operators_and_all_bounds(tmp_path):
    config = tmp_path / "stars.setools"
    contents = (
        "[MASK:preselect]\nMAG_AUTO < 21\n"
        "[MASK:stars]\nMAG_AUTO > 18.\nMAG_AUTO < 22.\nFLAGS == 0\n"
    )
    config.write_text(contents)
    ref = 'stars.setools#MASK:stars.MAG_AUTO = ["> 18.", "< 22."]'
    assert check_anchor_value(tmp_path, ref) is None
    flag_ref = 'stars.setools#MASK:stars.FLAGS = "== 0"'
    assert check_anchor_value(tmp_path, flag_ref) is None
    for altered in (
        contents.replace("> 18.", ">= 18."),
        contents.replace("MAG_AUTO < 22.\n", ""),
    ):
        config.write_text(altered)
        assert check_anchor_value(tmp_path, ref) is not None


def test_python_literals_and_dict_paths_without_importing(tmp_path):
    (tmp_path / "constants.py").write_text(
        "raise RuntimeError('must not execute')\n"
        "WIDTH: int = 51\nNOISE = 5e-4\n"
        "COMPLETENESS = {'exp_split': {'split_exp_runner': "
        "dict(expect=121, warn=True)}}\n"
        "class Model:\n    WIDTH = 53\n"
        "def fit():\n"
        "    limits = [-1.0, 1.0e3]\n"
        "    options = {'step': 0.01, 'dynamic': choose_at_runtime()}\n"
    )
    for selector, expected in (
        ("WIDTH", "51.0"), ("NOISE", "0.0005"), ("Model.WIDTH", "53"),
        ("fit.limits", "-1,1000"), ("fit.options[step]", "1e-2"),
        ("COMPLETENESS[exp_split.split_exp_runner.expect]", "121"),
        ("COMPLETENESS[exp_split.split_exp_runner.warn]", "True"),
    ):
        ref = f"constants.py::{selector} = {expected}"
        assert resolve_anchor(tmp_path, ref) is None
        assert check_anchor_value(tmp_path, ref) is None
    missing = "constants.py::fit.options[missing]"
    assert resolve_anchor(tmp_path, missing) is not None


@pytest.mark.parametrize(
    "contents, selector",
    [
        ("X = get_value()", "X"),
        ("X = 1 / 3", "X"),
        ("X = 1\nX = 2", "X"),
        ("X = 1\nX += 1", "X"),
        ("X, Y = 1, 2", "X"),
        ("def X():\n    return 1", "X"),
        ("X = {'a': 1, 'a': 2}", "X[a]"),
        ("X = dict(**other)", "X[a]"),
        ("X = {'a': 1, **other}", "X[a]"),
        ("X = {variable: 1}", "X[a]"),
        ("X = {'a': 1}\nX['a'] = 2", "X[a]"),
        ("X = {'a': 1}\nX['b']['c'] = 2", "X[a]"),
        ("X = {'a': 1}\nX['a'] += 1", "X[a]"),
        ("X = {'a': 1}\ndel X['a']", "X[a]"),
        ("X = 1\nX.attr = 2", "X"),
        ("def f():\n    X = {'a': 1}\n    X['a'] = 2", "f.X[a]"),
    ],
)
def test_python_nonliteral_or_ambiguous_values_fail_closed(
    tmp_path, contents, selector
):
    (tmp_path / "constants.py").write_text(contents + "\n")
    ref = f"constants.py::{selector} = 1"
    assert check_anchor_value(tmp_path, ref) is not None


@pytest.mark.parametrize(
    "reference",
    [
        "config.sex#KEY =", "config.sex#KEY == 1", "config.sex = 1",
        "config.sex#KEY = [1,", "config.sex#KEY = 1,,2",
        "config.sex#KEY = {'value': 1}", "config.sex#KEY = [[1]]",
        "../outside.sex#KEY = 1", "/outside.sex#KEY = 1",
        "config.sex#KEY = 1; config.sex#KEY = 2",
    ],
)
def test_malformed_assertions_are_errors(tmp_path, reference):
    (tmp_path / "config.sex").write_text("KEY 1\n")
    assert check_anchor_value(tmp_path, reference) is not None


def test_injected_drift_reports_decision_ref_expected_actual(tmp_path):
    config = tmp_path / "detect.sex"
    config.write_text("DETECT_THRESH 1.0\n")
    reference = "detect.sex#DETECT_THRESH = 1"
    record = record_with(reference)
    assert value_errors(tmp_path, record) == []

    config.write_text("DETECT_THRESH 1.5\n")
    assert resolve_anchor(tmp_path, reference) is None
    error, = value_errors(tmp_path, record)
    for detail in (
        "analyses.detection.decisions.threshold", "detect.sex#DETECT_THRESH",
        "expected", "1", "actual", "1.5",
    ):
        assert detail in error


def test_python_drift_is_not_hidden_by_a_cached_ast(tmp_path):
    target = tmp_path / "constants.py"
    reference = "constants.py::WIDTH = 51"
    target.write_text("WIDTH = 51\n")
    assert check_anchor_value(tmp_path, reference) is None
    target.write_text("WIDTH = 53\n")
    assert resolve_anchor(tmp_path, reference) is None
    assert check_anchor_value(tmp_path, reference) is not None


def test_malformed_anchor_sentence_cannot_silently_skip_values(tmp_path):
    record = {"decisions": {"width": {"rationale": "Anchor: config.sex#KEY = 1"}}}
    assert value_errors(tmp_path, record)


@pytest.mark.parametrize(
    "filename, contents, selector",
    [
        ("config.sex", "SATUR_KEY SATURATE\n# SATUR_LEVEL 50000\n", "SATUR_LEVEL"),
        ("config.ini", "[DEFAULT]\nA = 1\n[S]\nKEY = 1\n", "S.OTHER"),
        ("stars.setools", "[MASK:other]\nMASK_EXT == 0\n[MASK:stars]\n"
         "IMAFLAGS_ISO == 0\n# MASK_EXT == 0\n", "MASK:stars.MASK_EXT"),
    ],
)
def test_absent_passes_when_no_active_setting(tmp_path, filename, contents, selector):
    (tmp_path / filename).write_text(contents)
    ref = f"{filename}#{selector} = absent"
    assert resolve_anchor(tmp_path, ref) is None
    assert check_anchor_value(tmp_path, ref) is None


@pytest.mark.parametrize(
    "filename, contents, selector",
    [
        ("config.sex", "SATUR_LEVEL 50000\n", "SATUR_LEVEL"),
        ("config.sex", "SATUR_LEVEL\n", "SATUR_LEVEL"),
        ("config.ini", "[S]\nKEY = 1\n", "S.KEY"),
        ("config.ini", "[DEFAULT]\nKEY = 1\n[S]\n", "S.KEY"),
        ("stars.setools", "[MASK:stars]\nIMAFLAGS_ISO == 0\nMASK_EXT == 0\n",
         "MASK:stars.MASK_EXT"),
    ],
)
def test_absent_fails_on_an_active_setting(tmp_path, filename, contents, selector):
    (tmp_path / filename).write_text(contents)
    problem = check_anchor_value(tmp_path, f"{filename}#{selector} = absent")
    assert "expected no active setting" in problem


@pytest.mark.parametrize(
    "filename, contents, reference",
    [
        # A renamed section must not make every absence trivially true.
        ("config.ini", "[S]\nKEY = 1\n", "config.ini#RENAMED.KEY = absent"),
        ("stars.setools", "[MASK:stars]\nFLAGS == 0\n",
         "stars.setools#MASK:renamed.MASK_EXT = absent"),
        ("constants.py", "X = 1\n", "constants.py::Y = absent"),
        ("missing.sex", None, "missing.sex#KEY = absent"),
    ],
)
def test_absent_needs_a_real_file_and_section(tmp_path, filename, contents, reference):
    if contents is not None:
        (tmp_path / filename).write_text(contents)
    assert resolve_anchor(tmp_path, reference) is not None
    assert check_anchor_value(tmp_path, reference) is not None


def test_quoted_absent_is_ordinary_text(tmp_path):
    (tmp_path / "config.sex").write_text("KEY absent\n")
    assert check_anchor_value(tmp_path, 'config.sex#KEY = "absent"') is None
    assert check_anchor_value(tmp_path, "config.sex#KEY = absent") is not None


@pytest.mark.parametrize(
    "path, pattern, replacement",
    [
        ("workflow/config/cfis/config_tile_Sx.ini",
         "WEIGHT_IMAGE = True", "WEIGHT_IMAGE = False"),
        ("workflow/config/cfis/config_tile_Sx.ini",
         "MAKE_POST_PROCESS = True", "MAKE_POST_PROCESS = False"),
        ("workflow/config/cfis/star_selection.setools",
         "[MASK:star_selection]\n", "[MASK:star_selection]\nMASK_EXT == 0\n"),
        ("workflow/config/cfis/default_tile.sex",
         "SATUR_KEY", "SATUR_LEVEL      50000\nSATUR_KEY"),
        ("src/shapepipe/modules/ngmix_package/ngmix.py",
         "\n    boot = ngmix.metacal.",
         "\n    metacal_pars['step'] = 0.02\n    boot = ngmix.metacal."),
    ],
)
def test_real_record_catches_gate_and_absence_drift(
    tmp_path, path, pattern, replacement
):
    """Each mutation once passed the record; each must now be reported."""

    record = load_yaml(REPO_ROOT / "astra.yaml")
    for ref in {
        ref.split(" = ", 1)[0].split("#", 1)[0].split("::", 1)[0]
        for anchor in extract_anchors(record) for ref in anchor.references
    }:
        source = REPO_ROOT / ref
        if source.is_file():
            target = tmp_path / ref
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(source.read_text(encoding="utf-8"))
    assert value_errors(tmp_path, record) == []

    target = tmp_path / path
    text = target.read_text(encoding="utf-8")
    assert text.count(pattern) == 1
    target.write_text(text.replace(pattern, replacement))
    assert value_errors(tmp_path, record)
