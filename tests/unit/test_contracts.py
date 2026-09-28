"""Preserve the utilities import-boundary contract."""

import textwrap
from pathlib import Path

from tests.helpers.decisions import forbid_rules, import_violations


REPO_ROOT = Path(__file__).resolve().parents[2]


def _write(root, relative, text):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text), encoding="utf-8")


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
    assert all("bad.py:1" in problem and "no-up-imports" in problem
               for problem in violations)


def test_utilities_do_not_import_modules():
    contracts_file = REPO_ROOT / "src/shapepipe/utilities/CONTRACTS"
    rules = forbid_rules(contracts_file)
    assert [rule[0] for rule in rules] == ["utilities-do-not-import-modules"]

    violations = import_violations(REPO_ROOT / "src", rules)
    message = "Forbidden imports:\n - " + "\n - ".join(violations)
    assert not violations, message
