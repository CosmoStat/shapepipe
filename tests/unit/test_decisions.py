"""Decision-tag parsing and the static Values resolver."""

from pathlib import Path

import pytest

from tests.helpers.decisions import (
    decision_ids,
    decision_marker_errors,
    decision_markers,
    forbid_rules,
    import_violations,
    load_yaml,
    main,
    scan_tags,
    tag_errors,
    value_errors,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def _record(rationale="Values: THRESH = 1.", decision="choice"):
    return {"decisions": {decision: {"label": "Choice", "rationale": rationale}}}


def _write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _config_tag(path, decision="choice", content="THRESH 1\n"):
    return _write(path, f"# @sc [decision:{decision}]\n{content}")


def test_config_tags_govern_paragraph_and_section(tmp_path):
    config = _write(
        tmp_path / "settings.ini",
        "# @sc [decision:choice]\n[SCIENCE]\nKEY = 1\n# comment\nOTHER = 2\n\n"
        "# @sc [decision:choice]\n[OUTPUT]\nSAVE = True\n",
    )
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert len(tags) == 2
    science, output = tags
    assert science.site.path == "settings.ini"
    assert (science.site.start, science.site.end) == (2, 6)
    assert science.site.section == "SCIENCE"
    assert output.site.section == "OUTPUT"
    assert output.site.scope == "section"


def test_section_tag_covers_the_entire_section_across_nested_key_tags(tmp_path):
    _write(
        tmp_path / "settings.ini",
        "# @sc [decision:section_choice]\n[S]\nA = 1\n"
        "# @sc [decision:key_choice]\nB = 2\n\n"
        "# @sc [decision:next_choice]\n[N]\nC = 3\n",
    )
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    section_tag = next(tag for tag in tags if tag.decisions == ("section_choice",))
    assert section_tag.site.section == "S"
    assert section_tag.site.start <= 6 <= section_tag.site.end


def test_prose_mentions_of_sc_are_not_tags(tmp_path):
    _write(
        tmp_path / "mod.py",
        '"""An @sc citation points to a decision.\n\n'
        "A paragraph that explains the @sc syntax.\n\n"
        "@sc [decision:choice]\n\n\"\"\"\n",
    )
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert len(tags) == 1
    assert tags[0].decisions == ("choice",)


def test_python_declaration_statement_and_snakemake_tags(tmp_path):
    _write(
        tmp_path / "src" / "mod.py",
        '"""Module."""\n\n'
        "# @sc [decision:choice]\nWIDTH = 51\n\n"
        "def fit():\n"
        '    """Fit.\n\n'
        "    @sc [decision:choice,label:coupling] fit-coupling\n"
        "    The local fit constraint is preserved.\n"
        '    """\n'
        "    SCALE = 1.0\n",
    )
    _write(
        tmp_path / "workflow" / "rules.smk",
        "# @sc [decision:choice]\nrule detect:\n    output: 'catalogue'\n\n",
    )
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert {(tag.site.kind, tag.site.symbol) for tag in tags} == {
        ("python_statement", "WIDTH"),
        ("python_declaration", "fit"),
        ("snakemake", ""),
    }
    contract = next(tag for tag in tags if tag.ident)
    assert contract.ident == "fit-coupling"
    assert contract.prose == "The local fit constraint is preserved."


def test_tag_grammar_rejects_bad_metadata_missing_sites_and_duplicate_contract_ids(
    tmp_path,
):
    _write(
        tmp_path / "a.sex",
        "# @sc [decision choice] malformed\nKEY 1\n\n"
        "# @sc [decision:missing,label:x] same-id\n# Local prose.\nKEY 2\n\n"
        "# @sc [decision:choice,label:x] same-id\n# Local prose.\nKEY 3\n\n"
        "# @sc [decision:choice]",
    )
    tags, errors = scan_tags(tmp_path)

    assert len(tags) == 2
    assert any("malformed @sc metadata" in error for error in errors)
    assert any("duplicate local-contract id same-id" in error for error in errors)
    assert any("same-id" in error for error in errors)
    assert any("@sc tag governs no site" in error for error in errors)


def test_decision_citations_are_checked_in_both_directions(tmp_path):
    record = {
        "decisions": {"top": {}, "orphan": {}},
        "analyses": {"stage": {"decisions": {"inner": {}}}},
    }
    _config_tag(tmp_path / "a.sex", "top")
    _config_tag(tmp_path / "b.sex", "stage.inner")
    _config_tag(tmp_path / "c.sex", "not-real")
    tags, parse_errors = scan_tags(tmp_path)

    assert parse_errors == []
    assert decision_ids(record) == {"top", "orphan", "stage.inner"}
    errors = tag_errors(tags, record)
    assert any("unknown decision 'not-real'" in error for error in errors)
    assert any("orphan: decision has no tagged site" in error for error in errors)


def test_repeated_decision_metadata_cites_multiple_real_decisions(tmp_path):
    _write(
        tmp_path / "shared.sex",
        "# @sc [decision:first,decision:second]\nTHRESH 1\n",
    )
    record = {"decisions": {"first": {}, "second": {}}}
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert tags[0].decisions == ("first", "second")
    assert tag_errors(tags, record) == []


def test_value_change_is_reported_with_decision_ref_expected_and_actual(tmp_path):
    _config_tag(tmp_path / "detect.sex", content="THRESH 1.5\n")
    record = _record("Detection. Values: THRESH = 1.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    problems = value_errors(tmp_path, record, tags)
    assert len(problems) == 1
    for detail in ("choice", "THRESH", "expected", "1", "actual", "1.5"):
        assert detail in problems[0]


def test_key_movement_inside_tagged_paragraph_preserves_value(tmp_path):
    config = _config_tag(
        tmp_path / "detect.sex",
        content="# explanatory comment\nOTHER 2\nTHRESH 1\n",
    )
    record = _record("Detection. Values: THRESH = 1.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert value_errors(tmp_path, record, tags) == []
    config.write_text(
        "# @sc [decision:choice]\nTHRESH 1\nOTHER 2\n", encoding="utf-8"
    )
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, record, tags) == []


def test_removed_tag_orphans_decision(tmp_path):
    config = _config_tag(tmp_path / "detect.sex")
    record = _record()
    tags, _ = scan_tags(tmp_path)
    assert tag_errors(tags, record) == []

    config.write_text("THRESH 1\n", encoding="utf-8")
    tags, _ = scan_tags(tmp_path)
    assert tag_errors(tags, record) == ["choice: decision has no tagged site"]


def test_equal_values_at_multiple_sites_still_require_qualified_refs(tmp_path):
    _config_tag(tmp_path / "default.param", content="VIGNET(51,51)\n")
    _config_tag(tmp_path / "default_noimaflags.param", content="VIGNET(51,51)\n")
    record = _record("Stamps. Values: VIGNET = 51.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    problem, = value_errors(tmp_path, record, tags)
    assert "exactly one tagged site" in problem
    assert "2" in problem


def test_deleting_one_of_two_vignet_tagged_sites_fails_its_assertion(tmp_path):
    one = _config_tag(tmp_path / "default.param", content="VIGNET(51,51)\n")
    two = _config_tag(tmp_path / "default_noimaflags.param", content="VIGNET(51,51)\n")
    record = _record(
        "Stamps. Values: default.param#VIGNET = 51; "
        "default_noimaflags.param#VIGNET = 51."
    )
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, record, tags) == []

    two.unlink()
    tags, errors = scan_tags(tmp_path)
    assert not errors
    problems = value_errors(tmp_path, record, tags)
    assert len(problems) == 1
    assert "default_noimaflags.param#VIGNET" in problems[0]
    assert "found 0" in problems[0]
    assert one.exists()


def test_absent_key_requires_tagged_scope_and_detects_added_key(tmp_path):
    config = _write(
        tmp_path / "settings.sex",
        "# @sc [decision:choice,scope:file]\nOTHER 2\n",
    )
    record = _record("No fixed threshold. Values: settings.sex#THRESH = absent.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert value_errors(tmp_path, record, tags) == []

    config.write_text(
        "# @sc [decision:choice,scope:file]\nOTHER 2\nTHRESH 1\n",
        encoding="utf-8",
    )
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert "expected no active setting" in value_errors(tmp_path, record, tags)[0]


def test_absent_assertion_fails_when_its_scope_is_removed(tmp_path):
    config = _write(
        tmp_path / "settings.ini",
        "# @sc [decision:choice]\n[SCIENCE]\nOTHER = 2\n\n",
    )
    record = _record("No fixed threshold. Values: SCIENCE.THRESH = absent.")
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, record, tags) == []

    # The tag no longer governs the named section; it cannot assert absence.
    config.write_text(
        "# @sc [decision:choice]\n[OTHER]\nKEY = 2\n", encoding="utf-8"
    )
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    problem, = value_errors(tmp_path, record, tags)
    assert "SCIENCE.THRESH" in problem
    assert "found 0" in problem


def test_numeric_lists_and_astromatic_boolean_words_keep_reader_semantics(tmp_path):
    _config_tag(
        tmp_path / "values.sex",
        content="THRESH 5e-4\nAPERTURE 2.5, 3.5\nFLAG Y\n",
    )
    record = _record(
        "Values. Values: THRESH = 0.0005; APERTURE = [2.50, 3.500]; FLAG = True."
    )
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert value_errors(tmp_path, record, tags) == []


def test_ini_multiline_value_keeps_continuation_lines(tmp_path):
    _write(
        tmp_path / "config.ini",
        "# @sc [decision:choice]\n[SCIENCE]\nMODULE = alpha, beta,\n"
        "         gamma, delta\n",
    )
    record = _record("Runner chain. Values: SCIENCE.MODULE = alpha,beta,gamma,delta.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert value_errors(tmp_path, record, tags) == []


def test_config_boolean_semantics_and_setools_predicates(tmp_path):
    ini = _write(tmp_path / "config.ini", "# @sc [decision:choice]\n[SCIENCE]\nENABLED = 1\n")
    record = _record("Toggle. Values: SCIENCE.ENABLED = True.")
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, record, tags) == []

    ini.write_text("# @sc [decision:choice]\n[SCIENCE]\nENABLED = Y\n", encoding="utf-8")
    tags, _ = scan_tags(tmp_path)
    assert value_errors(tmp_path, record, tags)
    ini.write_text("# @sc [decision:choice]\n[SCIENCE]\nENABLED = False\n", encoding="utf-8")
    tags, _ = scan_tags(tmp_path)
    assert value_errors(tmp_path, record, tags)

    setools = _write(
        tmp_path / "stars.setools",
        "# @sc [decision:choice]\n[MASK:stars]\nMAG_AUTO > 18.\nMAG_AUTO < 22.\n",
    )
    cuts = _record('Star cut. Values: stars.setools#MASK:stars.MAG_AUTO = ["> 18.", "< 22."].')
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, cuts, tags) == []
    setools.write_text(
        "# @sc [decision:choice]\n[MASK:stars]\nMAG_AUTO >= 18.\nMAG_AUTO < 22.\n",
        encoding="utf-8",
    )
    tags, _ = scan_tags(tmp_path)
    assert value_errors(tmp_path, cuts, tags)


def test_python_value_reassignment_fails_closed(tmp_path):
    _write(
        tmp_path / "constants.py",
        'def fit():\n    """Fit.\n\n    @sc [decision:choice]\n    """\n'
        "    WIDTH = 51\n    WIDTH = 53\n",
    )
    record = _record("Stamp. Values: WIDTH = 51.")
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    problem, = value_errors(tmp_path, record, tags)
    assert "exactly one tagged site" in problem


def test_python_qualified_selector_ignores_other_tagged_scopes(tmp_path):
    _write(
        tmp_path / "constants.py",
        "# @sc [decision:choice]\nOTHER = 5\n\n"
        'def fit():\n    """Fit.\n\n    @sc [decision:choice]\n    """\n'
        "    PARAMS = {'limits': {'T': 1}}\n",
    )
    record = _record("Prior. Values: fit.PARAMS[limits.T] = 1.")
    tags, errors = scan_tags(tmp_path)

    assert errors == []
    assert value_errors(tmp_path, record, tags) == []


def test_python_literal_and_dict_selector_is_static(tmp_path):
    _write(
        tmp_path / "constants.py",
        "raise RuntimeError('must not execute')\n"
        "# @sc [decision:choice]\nCOMPLETENESS = {'run': {'expect': 40, 'warn': True}}\n",
    )
    record = _record("Counts. Values: COMPLETENESS[run.expect] = 40.")
    tags, errors = scan_tags(tmp_path)
    assert errors == []
    assert value_errors(tmp_path, record, tags) == []


def test_decision_markers_keep_unknown_id_check(tmp_path):
    _write(
        tmp_path / "tests" / "test_markers.py",
        "import pytest\npytestmark = [pytest.mark.decision('choice')]\n"
        "@pytest.mark.decision('missing')\ndef test_it():\n    pass\n",
    )
    markers, errors = decision_markers(tmp_path)

    assert errors == []
    assert [(marker.decision, marker.line) for marker in markers] == [
        ("choice", 2), ("missing", 3)
    ]
    assert decision_marker_errors(markers, _record()) == [
        "tests/test_markers.py:3: decision marker cites unknown decision 'missing'"
    ]


def test_cli_reports_decision_and_local_contract_for_a_location(tmp_path, capsys):
    _write(tmp_path / "astra.yaml", 'decisions:\n  choice:\n    label: Choice\n    rationale: "First sentence. Values: THRESH = 1."\n')
    _write(
        tmp_path / "detect.sex",
        "# @sc [decision:choice,label:coupling] threshold-coupling\n"
        "# The threshold must remain coupled.\nTHRESH 1\n",
    )

    assert main(["detect.sex:3", "--root", str(tmp_path)]) == 0
    output = capsys.readouterr().out
    assert "choice — Choice" in output
    assert "First sentence." in output
    assert "Values: THRESH = 1." in output
    assert "threshold-coupling" in output
    assert "The threshold must remain coupled." in output


def test_preserved_utilities_import_rule(tmp_path):
    contracts = _write(
        tmp_path / "pkg" / "utilities" / "CONTRACTS",
        "@cc no-up-imports\nforbid: pkg.utilities.* -> pkg.modules.*\n",
    )
    _write(tmp_path / "pkg" / "utilities" / "bad.py", "from ..modules import runner\n")
    rules = forbid_rules(contracts)
    assert rules == [("no-up-imports", "pkg.utilities.*", "pkg.modules.*")]
    assert len(import_violations(tmp_path, rules)) == 2


def test_decision_ids_cover_nested_analysis_paths():
    assert decision_ids(
        {"decisions": {"outer": {}}, "analyses": {"child": {"decisions": {"inner": {}}}}}
    ) == {"outer", "child.inner"}
