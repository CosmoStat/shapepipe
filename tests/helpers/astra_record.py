"""Reusable parsing and resolution helpers for ShapePipe's ASTRA record."""

import argparse
import ast
import configparser
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
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


def _split_assertion(reference):
    """Split the reserved, whitespace-delimited `` = `` (never a cut's ==)."""

    if ";" in reference:
        raise ValueError("semicolon is reserved for separating anchor refs")
    parts = re.split(r"\s+=\s*", reference.strip(), maxsplit=1)
    locator = parts[0]
    expected = parts[1].strip() if len(parts) == 2 else None
    if not locator or re.search(r"\s|=", locator):
        raise ValueError("expected a locator optionally followed by ' = value'")
    if expected is not None and (not expected or expected.startswith("=")):
        raise ValueError("expected a nonempty value after ' = '")
    return locator, expected


def _parse_reference(reference):
    reference, _ = _split_assertion(reference)
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

    try:
        kind, relative, selector = _parse_reference(reference)
        _, expected = _split_assertion(reference)
    except ValueError as error:
        return str(error)
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

    if kind == "code" and expected == ABSENT:
        return "absent assertions need a config key, not a code symbol"
    if kind == "code":
        if _is_snakemake_file(target):
            return _snakemake_symbol(text, selector)
        if target.suffix != ".py":
            return "code-symbol refs must name a .py file"
        try:
            tree = _python_tree(text)
            symbol, keys = _code_selector(selector)
            if not _has_symbol(tree, symbol):
                return f"no def/class/assignment target named {symbol!r}"
            if keys:
                _selected_python_node(tree, selector)
        except (SyntaxError, ValueError) as error:
            return f"cannot resolve Python selector: {error}"
        return None

    suffix = target.suffix.lower()
    if expected == ABSENT:
        return _absent_scope(text, selector, suffix)
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


ABSENT = "absent"
_LINE_SUFFIXES = {".sex", ".psfex", ".ww", ".param", ".conf", ".setools"}
_SECTIONED = {".ini", ".setools"}


def _absent_scope(text, selector, suffix):
    """An absent key still needs a real file type and, if sectioned, section.

    Without the section check a renamed section would make every absence
    trivially true.
    """

    if suffix not in _LINE_SUFFIXES | {".ini"}:
        return f"unsupported config-key file type {suffix or '(no extension)'}"
    if suffix not in _SECTIONED:
        return None
    if "." not in selector:
        return "sectioned config ref needs SECTION.KEY"
    section = selector.rsplit(".", 1)[0]
    if suffix == ".ini":
        try:
            parser = _ini_parser(text, strict=False)
        except configparser.Error as error:
            return f"cannot parse INI file: {error}"
        if section == parser.default_section or parser.has_section(section):
            return None
    elif any(
        line.strip() == f"[{section}]" for line in text.splitlines()
    ):
        return None
    return f"section {section!r} is missing"


def _ini_key(text, selector):
    if "." not in selector:
        return "INI config ref needs SECTION.KEY"
    section, key = selector.rsplit(".", 1)
    try:
        parser = _ini_parser(text, strict=False)
    except configparser.Error as error:
        return f"cannot parse INI file: {error}"
    if section != parser.default_section and not parser.has_section(section):
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


@lru_cache(maxsize=128)
def _bindings(scope):
    """Collect all bindings per name; value reads must not pick one silently."""

    result = {}

    def visit(node):
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            result.setdefault(node.name, []).append(node)
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
                for name in _target_names(target):
                    result.setdefault(name, []).append(node)
        if isinstance(node, ast.ExceptHandler) and node.name:
            result.setdefault(node.name, []).append(node)
        for child in ast.iter_child_nodes(node):
            visit(child)

    for statement in scope.body:
        visit(statement)
    return result


def _mutated_name(target):
    """Base name of ``NAME[...] =`` / ``NAME.attr =`` (nested included)."""

    while isinstance(target, (ast.Subscript, ast.Attribute)):
        target = target.value
        if isinstance(target, ast.Name):
            return target.id
    return None


@lru_cache(maxsize=128)
def _mutations(scope):
    """Names whose bound object is item- or attribute-assigned in ``scope``.

    Same lexical scope only, like ``_bindings``; method calls such as
    ``.update()`` and mutation from other scopes are not seen.
    """

    names = set()

    def visit(node):
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)
        ):
            return
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            targets = [node.target]
        elif isinstance(node, ast.Delete):
            targets = node.targets
        elif isinstance(node, (ast.For, ast.AsyncFor)):
            targets = [node.target]
        elif isinstance(node, (ast.With, ast.AsyncWith)):
            targets = [item.optional_vars for item in node.items]
        else:
            targets = []
        stack = [target for target in targets if target is not None]
        while stack:
            target = stack.pop()
            if isinstance(target, (ast.Tuple, ast.List)):
                stack.extend(target.elts)
            elif isinstance(target, ast.Starred):
                stack.append(target.value)
            else:
                name = _mutated_name(target)
                if name:
                    names.add(name)
        for child in ast.iter_child_nodes(node):
            visit(child)

    for statement in scope.body:
        visit(statement)
    return frozenset(names)


def _has_symbol(tree, symbol):
    scope = tree
    parts = symbol.split(".")
    for index, part in enumerate(parts):
        declarations = _bindings(scope).get(part)
        if not declarations:
            return False
        declaration = declarations[-1]
        if index == len(parts) - 1:
            return True
        if not isinstance(
            declaration,
            (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef),
        ):
            return False
        scope = declaration
    return False


def _ini_parser(text, *, strict=True):
    parser = configparser.ConfigParser(
        interpolation=None, strict=strict, allow_no_value=True
    )
    parser.optionxform = str
    parser.read_string(text)
    return parser


@lru_cache(maxsize=16)
def _python_tree(text):
    # Cache by source, not path: editing a file must invalidate the read.
    return ast.parse(text)


def _code_selector(selector):
    match = re.fullmatch(r"([\w.]+)(?:\[([\w.]+)\])?", selector)
    if not match or any(not p.isidentifier() for p in match[1].split(".")):
        raise ValueError(f"invalid Python selector {selector!r}")
    keys = tuple(match[2].split(".")) if match[2] else ()
    if any(not key.isidentifier() for key in keys):
        raise ValueError("dict paths need dot-separated identifier keys")
    return match[1], keys


def _dict_entry(node, key):
    """Select syntax, not a runtime value; never execute a dict() call."""

    if isinstance(node, ast.Dict):
        if any(
            not isinstance(k, ast.Constant) or not isinstance(k.value, str)
            for k in node.keys
        ):
            raise ValueError("dict selectors need literal string keys, no **")
        items = [(k.value, v) for k, v in zip(node.keys, node.values)]
    elif (
        isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        and node.func.id == "dict" and not node.args
        and all(k.arg is not None for k in node.keywords)
    ):
        items = [(k.arg, k.value) for k in node.keywords]
    else:
        raise ValueError("dict selectors need {...} or dict(key=value) syntax")
    names = [name for name, _ in items]
    if len(names) != len(set(names)):
        raise ValueError("ambiguous duplicate dict keys")
    if key not in names:
        raise ValueError(f"dict key {key!r} is missing")
    return dict(items)[key]


def _selected_python_node(tree, selector):
    symbol, keys = _code_selector(selector)
    scope = tree
    parts = symbol.split(".")
    for index, part in enumerate(parts):
        declarations = _bindings(scope).get(part, [])
        if len(declarations) != 1:
            raise ValueError(
                f"{symbol!r} needs one binding; found {len(declarations)} "
                f"for {part!r}"
            )
        node = declarations[0]
        if index == len(parts) - 1 and part in _mutations(scope):
            raise ValueError(
                f"{symbol!r} is item- or attribute-assigned after binding"
            )
        if index < len(parts) - 1:
            if not isinstance(
                node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            ):
                raise ValueError(f"{part!r} is not a lexical scope")
            scope = node
    if not isinstance(node, (ast.Assign, ast.AnnAssign)):
        raise ValueError(f"{symbol!r} is not a literal assignment")
    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
    if any(not isinstance(target, ast.Name) for target in targets):
        raise ValueError("value assertions need simple named assignment targets")
    node = node.value
    for key in keys:
        node = _dict_entry(node, key)
    return node


def _active_lines(text, selector, suffix):
    """Every active (uncommented) setting of the key, as (value, predicate)."""

    section = None
    if suffix == ".setools":
        if "." not in selector:
            raise ValueError("SETools config ref needs SECTION.KEY")
        section, key = selector.rsplit(".", 1)
    else:
        key = selector.rsplit(".", 1)[-1]
    pattern = re.compile(rf"^{re.escape(key)}(?=$|\s|=|\(|<|>)(.*)$")
    active = section is None
    values = []
    predicates = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if section is not None and line.startswith("[") and line.endswith("]"):
            active = line[1:-1].strip() == section
            continue
        match = pattern.fullmatch(line) if active else None
        if not match:
            continue
        value = match[1].strip()
        predicate = suffix == ".setools" and value.startswith(
            ("==", "!=", "<", ">")
        )
        if suffix == ".param" and value.startswith("(") and value.endswith(")"):
            value = value[1:-1]
        elif value.startswith("=") and not predicate:
            value = value[1:].strip()
        values.append(value)
        predicates.append(predicate)
    return values, predicates


def _line_value(text, selector, suffix):
    """Read active lines; SETools repeated predicates form an ordered list."""

    values, predicates = _active_lines(text, selector, suffix)
    if not values or any(not value for value in values):
        raise ValueError(f"no active value for {selector!r}")
    if len(values) > 1 and not all(predicates):
        raise ValueError(f"ambiguous active values for {selector!r}: {values!r}")
    return values if len(values) > 1 else values[0]


_NUMBER = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z")
# Boolean words each reader accepts. INI follows ConfigParser.getboolean,
# whose 1/0 spellings are handled at comparison (see _ini_bool); SExtractor
# and PSFEx add Y/N; Python literals are already bool, and the record spells
# them True/False. Elsewhere words stay text.
_INI_BOOLEANS = {
    "yes": True, "true": True, "on": True,
    "no": False, "false": False, "off": False,
}
_ASTROMATIC_BOOLEANS = {**_INI_BOOLEANS, "y": True, "n": False}
_PYTHON_BOOLEANS = {"true": True, "false": False}
_NO_BOOLEANS = {}


def _booleans_for(suffix):
    if suffix == ".ini":
        return _INI_BOOLEANS
    if suffix in {".sex", ".psfex"}:
        return _ASTROMATIC_BOOLEANS
    if suffix == ".py":
        return _PYTHON_BOOLEANS
    return _NO_BOOLEANS


def _ini_bool(want, got):
    """getboolean also reads 1/0; accept them only against a bool expectation."""

    if want[0] == "bool" and got[0] == "number" and got[1] in (0, 1):
        return "bool", got[1] == 1
    return got


def _normalise_value(value, booleans=_ASTROMATIC_BOOLEANS):
    """Use tagged atoms so boolean True cannot compare equal to number 1."""

    if isinstance(value, bool):
        return "bool", value
    if isinstance(value, (int, float)):
        number = Decimal(str(value))
        if not number.is_finite():
            raise ValueError("numeric values must be finite")
        return "number", number
    if isinstance(value, (list, tuple)):
        elements = tuple(_normalise_value(item, booleans) for item in value)
        if any(kind == "list" for kind, _ in elements):
            raise ValueError("only flat lists are supported")
        return "list", elements
    if not isinstance(value, str):
        raise ValueError("expected a number, boolean, string or flat list")
    value = value.strip()
    if not value:
        raise ValueError("empty values/list elements are not supported")
    if value[0] in "[{'\"":
        # BaseLoader keeps even 5e-4 and Y as strings, avoiding YAML 1.1's
        # inconsistent numeric/boolean coercions. It constructs no objects.
        parsed = yaml.load(value, Loader=yaml.BaseLoader)
        if isinstance(parsed, list):
            return _normalise_value(parsed, booleans)
        if not isinstance(parsed, str):
            raise ValueError("expected a scalar or flat list, not a mapping")
        # Quotes protect commas/operators; their contents are a single atom.
        value = parsed.strip()
    elif "," in value:
        return _normalise_value(value.split(","), booleans)
    if _NUMBER.fullmatch(value):
        return "number", Decimal(value)
    if value.lower() in booleans:
        return "bool", booleans[value.lower()]
    return "text", value


def _square_stamp(value):
    if value[0] == "number":
        return "list", (value, value)
    return value


def check_anchor_value(root, reference):
    """Return a diagnostic for a mismatched/unreadable assertion, else None.

    A reference without `` = value`` is location-only. Numbers compare
    exactly after decimal normalization, not with a tolerance. No imported
    code, environment expansion, function calls or expressions are evaluated.
    """

    expected = None
    actual = "<unreadable>"
    try:
        _, expected = _split_assertion(reference)
        if expected is None:
            return None
        kind, relative, selector = _parse_reference(reference)
        problem = resolve_anchor(root, reference)
        if problem:
            raise ValueError(problem)
        target = Path(root) / relative
        text = target.read_text(encoding="utf-8")
        suffix = target.suffix.lower()
        if expected == ABSENT:
            if suffix == ".ini":
                section, key = selector.rsplit(".", 1)
                parser = _ini_parser(text)
                if not parser.has_option(section, key):
                    return None
                actual = parser.get(section, key)
            else:
                active, _ = _active_lines(text, selector, suffix)
                if not active:
                    return None
                actual = active if len(active) > 1 else active[0]
            return f"expected no active setting, actual {actual!r}"
        if kind == "code" and suffix == ".py":
            node = _selected_python_node(_python_tree(text), selector)
            try:
                actual = ast.literal_eval(node)
            except (ValueError, TypeError) as error:
                raise ValueError(
                    "selected Python value is not a literal"
                ) from error
        elif kind == "config" and suffix == ".ini":
            section, key = selector.rsplit(".", 1)
            actual = _ini_parser(text).get(section, key)
            if actual is None:
                raise ValueError("no active value for INI key")
        elif kind == "config" and suffix in {
            ".sex", ".psfex", ".ww", ".param", ".conf", ".setools"
        }:
            actual = _line_value(text, selector, suffix)
        else:
            raise ValueError(
                "value assertions need a config key or Python assignment"
            )
        booleans = _booleans_for(suffix)
        want = _normalise_value(expected, booleans)
        got = _normalise_value(actual, booleans)
        if suffix == ".ini":
            got = _ini_bool(want, got)
        if (suffix, selector) in {(".param", "VIGNET"), (".psfex", "PSF_SIZE")}:
            want, got = _square_stamp(want), _square_stamp(got)
        if want == got:
            return None
        return f"expected {expected!r}, actual {actual!r}"
    except (
        ValueError, OSError, SyntaxError, configparser.Error, yaml.YAMLError
    ) as error:
        return f"expected {expected!r}, actual {actual!r}: {error}"


def value_errors(root, record):
    """Check assertions, naming decision, ref, expected and actual in errors."""

    errors = []
    for anchor in extract_anchors(record):
        if anchor.error:
            errors.append(f"{anchor.location}: {anchor.error}")
            continue
        for reference in anchor.references:
            problem = check_anchor_value(root, reference)
            if problem:
                errors.append(f"{anchor.location}: {reference}: {problem}")
    return errors


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
