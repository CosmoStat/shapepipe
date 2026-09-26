"""Reusable parsing and resolution helpers for ShapePipe's ASTRA record."""

import argparse
import ast
import configparser
from dataclasses import dataclass
import json
from pathlib import Path
import re
import subprocess
import sys

import yaml


@dataclass(frozen=True)
class Anchor:
    """An anchor sentence found in a YAML value."""

    location: str
    references: tuple[str, ...]
    error: str | None = None


def load_yaml(path):
    """Load YAML with PyYAML's safe loader."""

    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def _walk(value, location=""):
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{location}.{key}" if location else str(key)
            yield from _walk(child, path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk(child, f"{location}[{index}]")
    else:
        yield location, value


def _rationales(document):
    for location, value in _walk(document):
        if location.endswith(".rationale"):
            yield location, value


def extract_anchors(document):
    """Parse all ``Anchor:`` sentences and check every rationale has one."""

    anchors = []
    for location, value in _walk(document):
        if not isinstance(value, str) or "Anchor:" not in value:
            continue
        tail = value.split("Anchor:", 1)[1].strip()
        error = None
        refs = ()
        if value.count("Anchor:") != 1:
            count = value.count("Anchor:")
            error = f"expected one Anchor: marker, found {count}"
        elif not tail.endswith("."):
            error = "anchor sentence must end with a period"
        else:
            refs = tuple(part.strip() for part in tail[:-1].split(";"))
            if not refs or any(not ref for ref in refs):
                error = "anchor sentence contains an empty ref"
        anchors.append(Anchor(location, refs, error))

    for location, value in _rationales(document):
        if (
            not isinstance(value, str)
            or value.count("Anchor:") != 1
            or not value.rstrip().endswith(".")
        ):
            anchors.append(
                Anchor(
                    location,
                    (),
                    "rationale must end with exactly one Anchor: sentence",
                )
            )
    return anchors


def _parse_reference(reference):
    if "::" in reference:
        path, symbol = reference.split("::", 1)
        return "code", path, symbol
    if "#" in reference:
        path, key = reference.split("#", 1)
        return "config", path, key
    return "path", reference, ""


_SNAKEMAKE_SUFFIXES = {".smk"}
_SNAKEFILE_NAMES = {"Snakefile"}


def _is_snakemake_file(target):
    return target.suffix in _SNAKEMAKE_SUFFIXES or target.name in _SNAKEFILE_NAMES


def _snakemake_symbol(text, symbol):
    rule_pattern = re.compile(
        rf"^\s*(?:rule|checkpoint)\s+{re.escape(symbol)}\s*:", re.MULTILINE
    )
    if rule_pattern.search(text):
        return None
    def_pattern = re.compile(rf"^\s*def\s+{re.escape(symbol)}\(", re.MULTILINE)
    if def_pattern.search(text):
        return None
    return f"no rule/checkpoint/def named {symbol!r}"


def resolve_anchor(root, reference):
    """Return ``None`` if a reference resolves, otherwise a diagnostic."""

    kind, relative, selector = _parse_reference(reference)
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        return "path must be relative to the repository root"
    target = Path(root) / path
    if not target.exists():
        return "path does not exist"
    if kind == "path":
        return None
    if not target.is_file():
        return "code/config refs must name a file"

    try:
        text = target.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        return f"cannot read file: {error}"

    if kind == "code":
        if _is_snakemake_file(target):
            return _snakemake_symbol(text, selector)
        if target.suffix != ".py":
            return "code-symbol refs must name a .py file"
        try:
            tree = ast.parse(text, filename=str(target))
        except SyntaxError as error:
            return f"cannot parse Python file: {error}"
        if not _has_symbol(tree, selector):
            return f"no def/class/assignment target named {selector!r}"
        return None

    suffix = target.suffix.lower()
    if suffix == ".ini":
        return _ini_key(text, selector)
    if suffix == ".setools":
        return _setools_key(text, selector)
    if suffix in {".sex", ".psfex", ".ww", ".param", ".conf"}:
        key = selector.rsplit(".", 1)[-1]
        pattern = re.compile(rf"^\s*(?:#\s*)?{re.escape(key)}(?=$|\s|=|\()")
        if any(pattern.search(line) for line in text.splitlines()):
            return None
        return f"no line starts with key {key!r} (commented keys are allowed)"
    return f"unsupported config-key file type {suffix or '(no extension)'}"


def _ini_key(text, selector):
    if "." not in selector:
        return "INI config ref needs SECTION.KEY"
    section, key = selector.rsplit(".", 1)
    parser = configparser.ConfigParser(
        interpolation=None, strict=False, allow_no_value=True
    )
    try:
        parser.read_string(text)
    except configparser.Error as error:
        return f"cannot parse INI file: {error}"
    if not parser.has_section(section):
        return f"INI section {section!r} is missing"
    if not parser.has_option(section, key):
        return f"INI key {key!r} is missing from section {section!r}"
    return None


def _setools_key(text, selector):
    if "." not in selector:
        return "SETools config ref needs SECTION.KEY"
    section, key = selector.rsplit(".", 1)
    pattern = re.compile(rf"^\s*(?:#\s*)?{re.escape(key)}(?=$|\s|=|<|>)")
    active = False
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            active = stripped[1:-1].strip() == section
        elif active and pattern.search(line):
            return None
    return f"SETools key {key!r} is missing from section {section!r}"


def _target_names(target):
    if isinstance(target, ast.Name):
        return [target.id]
    if isinstance(target, ast.Attribute):
        return [target.attr]
    if isinstance(target, (ast.Tuple, ast.List)):
        return [name for item in target.elts for name in _target_names(item)]
    if isinstance(target, ast.Starred):
        return _target_names(target.value)
    return []


def _bindings(scope):
    """Collect declarations and assignment targets in one lexical scope."""

    result = {}

    def visit(node):
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            result[node.name] = node
            return
        if isinstance(node, ast.Lambda):
            return
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign, ast.NamedExpr)):
            targets = [node.target]
        elif isinstance(node, (ast.For, ast.AsyncFor)):
            targets = [node.target]
        elif isinstance(node, (ast.With, ast.AsyncWith)):
            targets = [item.optional_vars for item in node.items]
        else:
            targets = []
        for target in targets:
            if target is not None:
                result.update(dict.fromkeys(_target_names(target), node))
        if isinstance(node, ast.ExceptHandler) and node.name:
            result[node.name] = node
        for child in ast.iter_child_nodes(node):
            visit(child)

    for statement in scope.body:
        visit(statement)
    return result


def _has_symbol(tree, symbol):
    scope = tree
    parts = symbol.split(".")
    for index, part in enumerate(parts):
        declaration = _bindings(scope).get(part)
        if declaration is None:
            return False
        if index == len(parts) - 1:
            return True
        if not isinstance(
            declaration,
            (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef),
        ):
            return False
        scope = declaration
    return False


def universe_errors(record, universe):
    """Check scoped decision IDs and options against the ASTRA record."""

    record_decisions = _decisions(record)
    pinned = _decisions(universe)
    errors = []
    for location in sorted(pinned.keys() - record_decisions.keys()):
        errors.append(
            f"{location}: universe decision is absent from astra.yaml"
        )
    for location in sorted(record_decisions.keys() - pinned.keys()):
        errors.append(
            f"{location}: astra.yaml decision is not pinned in the universe"
        )
    for location in sorted(record_decisions.keys() & pinned.keys()):
        definition = record_decisions[location]
        options = (
            definition.get("options", {})
            if isinstance(definition, dict)
            else {}
        )
        if not isinstance(options, dict) or pinned[location] not in options:
            errors.append(
                f"{location}: pinned option {pinned[location]!r} is not in "
                "ASTRA options"
            )
    return errors


def _decisions(document, location=""):
    scope = document if isinstance(document, dict) else {}
    result = {}
    for decision_id, definition in (scope.get("decisions") or {}).items():
        key = (
            f"{location}.decisions.{decision_id}"
            if location
            else f"decisions.{decision_id}"
        )
        result[key] = definition
    for analysis_id, analysis in (scope.get("analyses") or {}).items():
        child = (
            f"{location}.analyses.{analysis_id}"
            if location
            else f"analyses.{analysis_id}"
        )
        if isinstance(analysis, dict):
            result.update(_decisions(analysis, child))
    return result


def _git_sha(root):
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            capture_output=True,
            check=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def build_report(root):
    """Resolve every anchor and universe pin under ``root`` into a report dict."""

    root = Path(root)
    astra_yaml = root / "astra.yaml"
    record = load_yaml(astra_yaml)
    anchors = extract_anchors(record)

    unresolved = []
    for anchor in anchors:
        if anchor.error:
            unresolved.append(
                {"location": anchor.location, "ref": None, "problem": anchor.error}
            )
            continue
        for reference in anchor.references:
            problem = resolve_anchor(root, reference)
            if problem:
                unresolved.append(
                    {
                        "location": anchor.location,
                        "ref": reference,
                        "problem": problem,
                    }
                )

    universe = load_yaml(root / "universes" / "committed.yaml")
    errors = universe_errors(record, universe)

    return {
        "astra_yaml": str(astra_yaml),
        "git_sha": _git_sha(root),
        "anchors_total": len(anchors),
        "unresolved": unresolved,
        "universe_errors": errors,
        "ok": not unresolved and not errors,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report", required=True, help="path to write the JSON report to"
    )
    parser.add_argument(
        "--root",
        default=Path(__file__).resolve().parents[2],
        help="repository root (default: repo root inferred from this file)",
    )
    args = parser.parse_args(argv)

    report = build_report(args.root)
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
