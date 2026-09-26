"""Keep ASTRA decision anchors and the committed universe resolvable."""

from pathlib import Path

from tests.helpers.astra_record import (
    extract_anchors,
    load_yaml,
    resolve_anchor,
    universe_errors,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


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
