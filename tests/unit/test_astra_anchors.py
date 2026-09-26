"""Keep ASTRA decision anchors and the committed universe resolvable."""

from pathlib import Path

from tests.helpers.astra_record import (
    extract_anchors,
    load_yaml,
    resolve_anchor,
    universe_errors,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_snakemake_rule_and_function_anchors_resolve(tmp_path):
    rules_dir = tmp_path / "workflow" / "rules"
    rules_dir.mkdir(parents=True)
    (rules_dir / "example.smk").write_text(
        "def tile_local(tile):\n    return tile\n\n\n"
        "rule tile_detect:\n    input: 'a'\n    output: 'b'\n",
        encoding="utf-8",
    )

    assert resolve_anchor(tmp_path, "workflow/rules/example.smk::tile_detect") is None
    assert resolve_anchor(tmp_path, "workflow/rules/example.smk::tile_local") is None

    problem = resolve_anchor(tmp_path, "workflow/rules/example.smk::no_such_rule")
    assert problem == "no rule/checkpoint/def named 'no_such_rule'"


def test_snakefile_rule_anchor_resolves(tmp_path):
    (tmp_path / "workflow").mkdir()
    (tmp_path / "workflow" / "Snakefile").write_text(
        "checkpoint plan:\n    input: 'a'\n", encoding="utf-8"
    )

    assert resolve_anchor(tmp_path, "workflow/Snakefile::plan") is None


def test_every_astra_anchor_resolves():
    record = load_yaml(REPO_ROOT / "astra.yaml")
    anchors = extract_anchors(record)
    errors = []

    assert anchors, "astra.yaml contains no Anchor: sentences"
    for anchor in anchors:
        if anchor.error:
            errors.append(f"{anchor.location}: {anchor.error}")
            continue
        for reference in anchor.references:
            problem = resolve_anchor(REPO_ROOT, reference)
            if problem:
                errors.append(
                    f"{anchor.location}: {reference}: {problem}"
                )

    message = (
        "Unresolved ASTRA anchors or rationales:\n - "
        + "\n - ".join(errors)
    )
    assert not errors, message


def test_committed_universe_matches_astra_decisions():
    record = load_yaml(REPO_ROOT / "astra.yaml")
    universe = load_yaml(REPO_ROOT / "universes" / "committed.yaml")

    errors = universe_errors(record, universe)

    message = (
        "ASTRA / committed universe mismatch:\n - "
        + "\n - ".join(errors)
    )
    assert not errors, message
