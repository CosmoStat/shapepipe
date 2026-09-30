"""The ``--diff`` report: decisions a branch touches, and their record."""

import subprocess

import pytest

from tests.helpers.decisions import diff_report

_RECORD = """\
decisions:
  choice:
    rationale: Choice rationale.
  untouched:
    rationale: Untouched rationale.
analyses:
  sub:
    decisions:
      inner:
        rationale: Inner rationale.
"""
_CONFIG = """\
# @sc [decision:choice]
THRESH 1

# @sc [decision:untouched]
OTHER 7

# @sc [decision:sub.inner]
INNER 3
"""
_LINE = {"choice": "THRESH 1", "untouched": "OTHER 7", "sub.inner": "INNER 3"}
_RATIONALE = {
    "choice": "Choice rationale.",
    "untouched": "Untouched rationale.",
    "sub.inner": "Inner rationale.",
}
_UNCHANGED = "record unchanged — check the rationale still holds"
CONFIG = "config/settings.sex"


def _git(repo, *args):
    return subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@t", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _edit(repo, relative, old, new):
    path = repo / relative
    path.write_text(path.read_text().replace(old, new))


@pytest.fixture
def repo(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "astra.yaml").write_text(_RECORD)
    (tmp_path / CONFIG).write_text(_CONFIG)
    _git(tmp_path, "init", "--quiet", "-b", "base")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "--quiet", "-m", "base")
    _git(tmp_path, "checkout", "--quiet", "-b", "feature")
    return tmp_path


def _commit(repo):
    _git(repo, "commit", "--quiet", "-am", "change")


def _rows(repo):
    return [
        line
        for line in diff_report(repo, "base").splitlines()
        if line.startswith("| `")
    ]


@pytest.mark.parametrize(
    "changed, amended, verdict",
    [
        ("choice", None, _UNCHANGED),
        ("choice", "untouched", _UNCHANGED),
        ("choice", "choice", "record amended"),
        ("sub.inner", "sub.inner", "record amended"),
        ("sub.inner", "choice", _UNCHANGED),
    ],
)
def test_changed_site_reported_with_its_own_record_verdict(
    repo, changed, amended, verdict
):
    _edit(repo, CONFIG, _LINE[changed], _LINE[changed] + "0")
    if amended:
        _edit(repo, "astra.yaml", _RATIONALE[amended], "Rewritten.")
    _commit(repo)
    line = _CONFIG.splitlines().index(_LINE[changed]) + 1
    assert _rows(repo) == [
        f"| `{changed}` | `{CONFIG}:{line}-{line}` | {verdict} |"
    ]


def test_deleted_site_mapped_against_the_base_file(repo):
    _edit(repo, CONFIG, "# @sc [decision:choice]\nTHRESH 1\n\n", "")
    _commit(repo)
    assert _rows(repo) == [
        f"| `choice` | `{CONFIG}:2-2` (base) | {_UNCHANGED} |"
    ]


def test_deleted_file_reports_every_site(repo):
    _git(repo, "rm", "--quiet", CONFIG)
    _commit(repo)
    assert [row.split(" | ")[0] for row in _rows(repo)] == [
        "| `choice`",
        "| `sub.inner`",
        "| `untouched`",
    ]


def test_base_advanced_past_the_merge_base(repo):
    """Base-branch line shifts and record edits are not the branch's."""
    _edit(repo, CONFIG, "OTHER 7", "OTHER 8")
    _commit(repo)
    _git(repo, "checkout", "--quiet", "base")
    _edit(
        repo,
        CONFIG,
        "# @sc [decision:choice]",
        "# header\n# @sc [decision:choice]",
    )
    _edit(repo, "astra.yaml", "Untouched rationale.", "Rewritten on base.")
    _commit(repo)
    _git(repo, "checkout", "--quiet", "feature")
    assert _rows(repo) == [f"| `untouched` | `{CONFIG}:5-5` | {_UNCHANGED} |"]


def test_rename_and_binary_touch_nothing(repo):
    _git(repo, "mv", CONFIG, "config/renamed.sex")
    (repo / "blob.bin").write_bytes(b"\0\1\2")
    _git(repo, "add", "blob.bin")
    _commit(repo)
    assert "No tagged decision site overlaps" in diff_report(repo, "base")
