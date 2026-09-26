"""Keep the @sc contracts well-formed and tied to the ASTRA decision record."""

from functools import cache
from pathlib import Path
import textwrap

import pytest

from tests.helpers.astra_record import load_yaml
from tests.helpers.contracts import (
    collect,
    coverage_report,
    decision_errors,
    decision_ids,
    forbid_rules,
    governed_refs,
    governs_errors,
    import_violations,
)

REPO_ROOT = Path(__file__).resolve().parents[2]

RECORD = {
    "decisions": {"top_choice": {}},
    "analyses": {"stage": {"decisions": {"inner_choice": {}}}},
}


@cache
def _repository():
    record = load_yaml(REPO_ROOT / "astra.yaml")
    contracts, errors = collect(REPO_ROOT)
    return record, contracts, errors


def _write(root, relative, text):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text), encoding="utf-8")


def test_parser_reads_docstrings_contracts_files_and_snakemake(tmp_path):
    _write(tmp_path, "src/pkg/mod.py", '''
        """Module."""


        class Thing:
            def method(self):
                """Do it.

                @sc [decision:stage.inner_choice,label:convention] method-rule
                The method keeps its promise.
                Second prose line.

                Returns
                -------
                None
                """
        ''')
    _write(tmp_path, "src/pkg/CONTRACTS", """
        @cc boundary-rule
        forbid: pkg.a.* -> pkg.b.*
        """)
    _write(tmp_path, "workflow/rules/x.smk", """
        # @sc [decision:top_choice] rule-contract
        # The rule keeps its promise.

        rule x:
            output: "a"
        """)

    contracts, errors = collect(tmp_path)

    assert errors == []
    by_id = {c.id: c for c in contracts}
    assert set(by_id) == {"method-rule", "boundary-rule", "rule-contract"}
    method = by_id["method-rule"]
    assert method.scope == "Thing.method"
    assert method.meta == {
        "decision": "stage.inner_choice",
        "label": "convention",
    }
    assert method.prose == "The method keeps its promise. Second prose line."
    assert method.line == 9
    assert by_id["boundary-rule"].scope == "src/pkg"
    assert decision_errors(contracts, RECORD) == []


def test_malformed_meta_and_missing_id_are_errors(tmp_path):
    _write(tmp_path, "src/mod.py", '''
        def f():
            """F.

            @sc [decision stage.inner_choice] bad-meta
            Prose.

            @sc [label:x]
            Prose.

            @sc two words
            Prose.
            """
        ''')

    contracts, errors = collect(tmp_path)

    assert contracts == []
    assert len(errors) == 3
    assert "malformed contract metadata" in errors[0]
    assert "missing contract id" in errors[1]
    assert "malformed contract line" in errors[2]


def test_duplicate_ids_and_comment_contracts_are_errors(tmp_path):
    _write(tmp_path, "src/a.py", '''
        def f():
            """F.

            @sc same-id
            Prose.
            """
        ''')
    _write(tmp_path, "scripts/b.py", '''
        def g():
            """G.

            @sc same-id
            Prose.
            """
            # @sc hidden-id
            return None
        ''')

    _, errors = collect(tmp_path)

    assert any("duplicate contract id same-id" in e for e in errors)
    assert any("contract in a comment" in e for e in errors)


def test_unknown_decision_is_an_error(tmp_path):
    _write(tmp_path, "src/mod.py", '''
        def f():
            """F.

            @sc [decision:inner_choice] undotted
            Sub-analysis decisions need their analysis prefix.

            @sc [decision:no_such_choice] unknown
            Prose.
            """
        ''')

    contracts, errors = collect(tmp_path)

    assert errors == []
    problems = decision_errors(contracts, RECORD)
    assert len(problems) == 2
    assert "'inner_choice'" in problems[0]
    assert "'no_such_choice'" in problems[1]
    assert decision_ids(RECORD) == {"top_choice", "stage.inner_choice"}


def test_governs_resolves_all_refs_relative_to_the_contract_file(tmp_path):
    """A multi-key coupling must not lose refs or resolve them from cwd."""

    base = "workflow/config/cfis"
    _write(tmp_path, f"{base}/default.sex", "DEBLEND_MINCONT 0.002\n")
    _write(tmp_path, f"{base}/default.psfex", "PSF_SIZE 51,51\n")
    _write(tmp_path, f"{base}/stamps.ini", "[STAMP]\nSIZE = 51\n")
    _write(tmp_path, f"{base}/default.param", "VIGNET(51,51)\n")
    _write(tmp_path, f"{base}/stars.setools", "[MASK:stars]\nFLAGS == 0\n")
    _write(tmp_path, f"{base}/kernel.conv", "CONV NORM\n1 2 1\n")
    # Bare-file refs also cover formats the shared resolver cannot select into.
    _write(tmp_path, "workflow/config.yaml", "psf_model: psfex\n")
    refs = (
        "default.sex#DEBLEND_MINCONT",
        "default.psfex#PSF_SIZE",
        "stamps.ini#STAMP.SIZE",
        "default.param#VIGNET",
        "stars.setools#MASK:stars.FLAGS",
        "kernel.conv",
    )
    _write(tmp_path, f"{base}/CONTRACTS", f"""
        @sc [decision:top_choice,governs:{';'.join(refs)}] coupled-config
        Keep the apertures coupled.
        Ordinary prose is not an import rule.
        """)
    _write(tmp_path, "workflow/CONTRACTS", """
        @sc [decision:stage.inner_choice,governs:config.yaml] model-choice
        Select matching exposure and tile models.
        """)

    contracts, errors = collect(tmp_path)
    by_id = {c.id: c for c in contracts}
    coupled = by_id["coupled-config"]

    assert errors == []
    assert coupled.meta["governs"] == ";".join(refs)
    assert coupled.scope == base
    assert coupled.line == 2
    assert coupled.prose == (
        "Keep the apertures coupled. Ordinary prose is not an import rule."
    )
    assert governed_refs(coupled) == tuple(f"{base}/{ref}" for ref in refs)
    assert governs_errors(contracts, tmp_path) == []
    assert decision_errors(contracts, RECORD) == []
    assert forbid_rules(tmp_path / base / "CONTRACTS") == []


@pytest.mark.parametrize("bad_ref", [
    "missing.sex#KEY",
    "default.sex#MISSING",
    "stamps.ini#STAMP.MISSING",
    "stamps.ini#SIZE",  # INI selectors need SECTION.KEY.
    "stars.setools#MASK:missing.FLAGS",
    "../default.sex#KEY",  # No escape from the governing directory.
    "/absolute/default.sex#KEY",
])
@pytest.mark.parametrize("bad_first", [True, False])
def test_every_unresolvable_governs_ref_is_an_error(
    tmp_path, bad_ref, bad_first
):
    """Checking only the first/last ref silently loses part of a coupling."""

    base = "workflow/config"
    _write(tmp_path, f"{base}/default.sex", "KEY 1\n")
    _write(tmp_path, f"{base}/stamps.ini", "[STAMP]\nSIZE = 51\n")
    _write(tmp_path, f"{base}/stars.setools", "[MASK:stars]\nFLAGS == 0\n")
    refs = [bad_ref, "default.sex#KEY"]
    if not bad_first:
        refs.reverse()
    _write(tmp_path, f"{base}/CONTRACTS", f"""
        @sc [governs:{';'.join(refs)}] broken-coupling
        Prose.
        """)

    contracts, errors = collect(tmp_path)
    problems = governs_errors(contracts, tmp_path)

    assert errors == []
    assert len(problems) == 1
    assert f"{base}/CONTRACTS:2" in problems[0]
    assert "broken-coupling" in problems[0]
    assert bad_ref in problems[0]


@pytest.mark.parametrize("metadata, message", [
    ("governs:", "malformed contract metadata"),
    ("governs:a.sex#KEY b.sex#KEY", "malformed contract metadata"),
    ("governs:a.sex#KEY,b.sex#KEY", "malformed contract metadata"),
    ("governs:;a.sex#KEY", "empty governs ref"),
    ("governs:a.sex#KEY;", "empty governs ref"),
    ("governs:a.sex#KEY;;b.sex#KEY", "empty governs ref"),
    ("governs:a.sex#KEY,governs:b.sex#KEY", "duplicate metadata key"),
])
def test_malformed_governs_does_not_silently_drop_refs(
    tmp_path, metadata, message
):
    _write(tmp_path, "workflow/config/CONTRACTS", f"""
        @sc [{metadata}] malformed-coupling
        Prose.
        """)

    contracts, errors = collect(tmp_path)

    assert contracts == []
    assert len(errors) == 1
    assert "workflow/config/CONTRACTS:2" in errors[0]
    assert message in errors[0]


def test_config_coverage_uses_governed_refs_not_the_sidecar_path(tmp_path):
    _write(tmp_path, "workflow/config/CONTRACTS", """
        @sc [decision:top_choice,governs:a.sex#KEY;b.ini#S.SIZE] config-coupling
        Prose.

        @sc [decision:stage.inner_choice,governs:a.sex#OTHER] off-record-key
        A different key in the same file is not the anchored key.

        @sc [governs:config.yaml] whole-file
        Prose.
        """)
    record = {
        "decisions": {
            "top_choice": {
                "rationale": "Anchor: workflow/config/b.ini#S.SIZE."
            },
            "uncovered": {},
        },
        "analyses": {"stage": {"decisions": {"inner_choice": {
            "rationale": "Anchor: workflow/config/a.sex#KEY."
        }}}},
        "description": "Anchor: workflow/config/config.yaml.",
    }

    contracts, errors = collect(tmp_path)
    uncovered, unanchored = coverage_report(contracts, record)

    assert errors == []
    assert uncovered == ["uncovered"]
    assert [c.id for c in unanchored] == ["off-record-key"]


def test_repository_contracts_are_valid_and_cite_real_decisions():
    record, contracts, errors = _repository()
    errors = (
        errors
        + decision_errors(contracts, record)
        + governs_errors(contracts, REPO_ROOT)
    )

    message = "Contract problems:\n - " + "\n - ".join(errors)
    assert not errors, message


def test_contract_coverage_report():
    """Report-only: print record decisions and contracts that lack a partner."""

    record, contracts, _ = _repository()
    uncovered, unanchored = coverage_report(contracts, record)

    print(f"\n{len(contracts)} contracts; "
          f"{len(uncovered)} decisions cited by no contract:")
    for decision in uncovered:
        print(f"  {decision}")
    print(f"{len(unanchored)} @sc contracts off the record's anchors:")
    for contract in unanchored:
        targets = governed_refs(contract)
        where = "; ".join(targets) if targets else (
            f"{contract.path}::{contract.scope}"
        )
        print(f"  {contract.id} at {where}")


def test_forbidden_import_is_found(tmp_path):
    _write(tmp_path, "pkg/utilities/CONTRACTS", """
        @cc no-up-imports
        forbid: pkg.utilities.* -> pkg.modules.*
        """)
    _write(tmp_path, "pkg/utilities/good.py", "import os\n")
    _write(tmp_path, "pkg/utilities/bad.py", "from ..modules import runner\n")
    _write(tmp_path, "pkg/modules/runner.py", "from pkg.utilities import good\n")

    rules = forbid_rules(tmp_path / "pkg/utilities/CONTRACTS")
    violations = import_violations(tmp_path, rules)

    assert rules == [("no-up-imports", "pkg.utilities.*", "pkg.modules.*")]
    assert len(violations) == 2
    assert all("bad.py:1" in v and "no-up-imports" in v for v in violations)


def test_utilities_do_not_import_modules():
    contracts_file = REPO_ROOT / "src/shapepipe/utilities/CONTRACTS"
    rules = forbid_rules(contracts_file)
    assert [rule[0] for rule in rules] == ["utilities-do-not-import-modules"]

    violations = import_violations(REPO_ROOT / "src", rules)

    message = "Forbidden imports:\n - " + "\n - ".join(violations)
    assert not violations, message
