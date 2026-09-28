"""Parse ShapePipe's decision tags, ASTRA Values and local contracts.

The value grammar behind ``Values: ref = value`` is:

* Values are numbers, boolean words, strings (quote expressions), or flat
  comma lists, optionally bracketed. Semicolons separate refs. Decimal
  equality is exact (1 = 1.0, 5e-4 = 0.0005), without rounding.
* Boolean words follow the file's reader, case-insensitively: INI as
  ConfigParser.getboolean (yes/true/on/1, no/false/off/0; Y/N are text),
  .sex/.psfex also Y/N, Python True/False; elsewhere words are text and
  1/0 are numbers. Outer whitespace is trimmed; other strings are
  case-sensitive. Lists preserve order and length. Only .param VIGNET and
  .psfex PSF_SIZE accept square-size shorthand: 51 = 51,51.
* ``= absent`` asserts that a config key has no active line anywhere in the
  file or named section. The same decision must tag at least one site in that
  file (and, for INI/SETools, in the named section or at file scope). INI
  absence checks include values inherited from ``[DEFAULT]``. Quote it
  ("absent") to mean the text.
* .setools refs use ``SECTION.KEY``. Predicates keep their operators as
  quoted text; repeated cuts on one key are an ordered list, e.g.
  ``MAG_AUTO = ["> 18.", "< 22."]``. Expressions compare as text.
* Python selectors may append ``[key.subkey]`` to a named assignment
  (identifier-like string dict keys); only the selected literal is read,
  including ``dict(key=value)`` syntax. Ambiguous bindings fail, including
  ``NAME[...] =`` or ``NAME.attr =`` in the same scope; ``.update()`` calls
  and mutation from other scopes are not seen. No imports, calls,
  arithmetic, argument defaults, environment expansion or implicit tool
  defaults are evaluated: the check stays static rather than becoming a
  second pipeline runtime.
* Option ids are stable references and are never parsed for values.
"""

import argparse
import ast
import configparser
import fnmatch
import re
import sys
import tokenize
import warnings
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
from io import StringIO
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Site:
    """One code declaration, statement, config paragraph, or config scope."""

    path: str
    start: int
    end: int
    kind: str
    symbol: str = ""
    section: str = ""
    scope: str = ""

    @property
    def identity(self):
        return (self.path, self.start, self.end, self.kind, self.symbol, self.section)


@dataclass(frozen=True)
class Tag:
    """One parsed ``@sc`` tag and the site it governs."""

    path: str
    line: int
    decisions: tuple[str, ...]
    ident: str | None
    meta: dict
    prose: str
    site: Site | None


@dataclass(frozen=True)
class DecisionMarker:
    """One pytest marker linking a test to an ASTRA decision."""

    decision: str
    path: str
    line: int


def load_yaml(path):
    """Load YAML with PyYAML's safe loader."""

    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


ABSENT = "absent"
_NO_SETTING = object()
_TAG_LINE = re.compile(r"^\s*@sc(?:\s+\[([^\]]*)\])?(?:\s+(\S+))?\s*$")
_META_KEY = re.compile(r"[A-Za-z_][\w.-]*:[^,\s]+\Z")
_ID = re.compile(r"[A-Za-z_][\w.-]*\Z")
_CONFIG_SUFFIXES = {
    ".ini", ".sex", ".param", ".conv", ".psfex", ".setools",
    ".ww", ".conf", ".yaml", ".yml",
}
_SKIP = {".git", ".venv", "venv", "__pycache__", "node_modules", ".felt"}


def _tag_line(text):
    """Parse one ``@sc`` line into metadata and an optional local id."""

    match = _TAG_LINE.fullmatch(text.strip())
    if not match:
        return None, None, f"malformed @sc tag: {text.strip()}"
    meta_text, ident = match.groups()
    meta = {}
    if meta_text is not None:
        pairs = [part.strip() for part in meta_text.split(",")]
        if not pairs or any(not _META_KEY.fullmatch(part) for part in pairs):
            return None, None, f"malformed @sc metadata: {meta_text}"
        for pair in pairs:
            key, value = pair.split(":", 1)
            if key != "decision" and key in meta:
                return None, None, f"duplicate @sc metadata key: {key}"
            meta.setdefault(key, []).append(value)
    if ident is not None and not _ID.fullmatch(ident):
        return None, None, f"malformed local-contract id {ident!r}"
    if "scope" in meta and meta["scope"] != ["file"]:
        return None, None, "scope metadata must be scope:file"
    if "decision" not in meta and ident is None:
        return None, None, "@sc needs a decision citation or local-contract id"
    return meta, ident, None


def _comment_body(line):
    stripped = line.lstrip()
    if not stripped.startswith("#"):
        return None
    return stripped[1:].lstrip()


def _comment_prose(lines, index, end=None):
    stop = len(lines) if end is None else end
    prose = []
    for line in lines[index + 1:stop]:
        if not line.strip():
            break
        body = _comment_body(line)
        if body is None or body.startswith("@sc"):
            break
        if body:
            prose.append(body.strip())
    return " ".join(prose)


def _ast_symbol(node, parents):
    parts = [node.name]
    parent = parents.get(node)
    while parent is not None and not isinstance(parent, ast.Module):
        if isinstance(parent, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            parts.append(parent.name)
        parent = parents.get(parent)
    return ".".join(reversed(parts))


def _python_docstring_tags(path, source, tree):
    lines = source.splitlines()
    parents = {child: parent for parent in ast.walk(tree)
               for child in ast.iter_child_nodes(parent)}
    owners = [(tree, "")]
    owners.extend((node, _ast_symbol(node, parents)) for node in ast.walk(tree)
                  if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)))
    tags, errors, consumed = [], [], set()
    for owner, symbol in owners:
        body = getattr(owner, "body", [])
        if not body or not isinstance(body[0], ast.Expr) or not isinstance(
            body[0].value, ast.Constant
        ) or not isinstance(body[0].value.value, str):
            continue
        doc_node = body[0].value
        for offset, doc_line in enumerate(doc_node.value.splitlines()):
            if not doc_line.strip().startswith("@sc"):
                continue
            line_no = doc_node.lineno + offset
            consumed.add(line_no)
            meta, ident, error = _tag_line(doc_line)
            if error:
                errors.append(f"{path}:{line_no}: {error}")
                continue
            prose = []
            for part in doc_node.value.splitlines()[offset + 1:]:
                if not part.strip() or part.strip().startswith("@sc"):
                    break
                prose.append(part.strip())
            prose_text = " ".join(prose)
            if ident is not None and not prose_text:
                errors.append(f"{path}:{line_no}: local contract {ident} has no prose")
            site = Site(path, getattr(owner, "lineno", 1),
                        getattr(owner, "end_lineno", len(lines)),
                        "python_declaration", symbol=symbol)
            tags.append(Tag(path, line_no, tuple(meta.get("decision", ())), ident,
                            {k: v for k, v in meta.items() if k != "decision"},
                            prose_text, site))
    return tags, errors, consumed


def _comment_site(path, lines, index, meta, tree=None, *, snakemake=False):
    line_no = index + 1
    if meta.get("scope") == ["file"]:
        if tree is not None or snakemake:
            return None
        has_content = any(
            line.strip() and _comment_body(line) is None for line in lines
        )
        if not has_content:
            return None
        return Site(path, 1, len(lines), "config", scope="file")
    if snakemake:
        for n in range(index + 1, len(lines)):
            line = lines[n]
            if not line.strip():
                return None
            if line.lstrip().startswith("#"):
                continue
            start = n + 1
            if re.match(r"^\s*(?:rule|checkpoint)\s+[\w.-]+\s*:", line):
                end = len(lines)
                for j in range(n + 1, len(lines)):
                    if lines[j].strip() and not lines[j].startswith((" ", "\t", "#")):
                        end = j
                        break
                return Site(path, start, end, "snakemake")
            return Site(path, start, start, "snakemake")
        return None
    if tree is not None:
        node = next((n for n in tree.body if getattr(n, "lineno", 0) > line_no), None)
        if node is None:
            return None
        for between in lines[index + 1:node.lineno - 1]:
            if not between.strip() or not between.lstrip().startswith("#"):
                return None
        names = []
        for target in getattr(node, "targets", []):
            names.extend(_target_names(target))
        if hasattr(node, "target"):
            names.extend(_target_names(node.target))
        return Site(path, node.lineno, getattr(node, "end_lineno", node.lineno),
                    "python_statement", symbol=names[0] if names else "")

    # Config paragraphs end at blank lines; comments before/inside a run are
    # skipped, but a blank before the first setting leaves an empty site.
    for n in range(index + 1, len(lines)):
        line = lines[n]
        if not line.strip():
            return None
        comment = _comment_body(line)
        if comment is not None:
            if comment.startswith("@sc"):
                return None
            continue
        start = n + 1
        header = re.match(r"^\s*\[([^]]+)\]\s*(?:[#;].*)?$", line)
        if header:
            end = len(lines)
            for j in range(n + 1, len(lines)):
                if re.match(r"^\s*\[[^]]+\]\s*(?:[#;].*)?$", lines[j]):
                    end = j
                    k = j - 1
                    while k > n and _comment_body(lines[k]) is not None:
                        if _comment_body(lines[k]).startswith("@sc"):
                            end = k
                            break
                        k -= 1
                    break
            return Site(path, start, end, "config", section=header.group(1).strip(), scope="section")
        end = len(lines)
        for j in range(n + 1, len(lines)):
            if not lines[j].strip():
                end = j
                break
            comment = _comment_body(lines[j])
            if comment is not None and comment.startswith("@sc"):
                end = j
                break
        suffix = Path(path).suffix.lower()
        section = ""
        if suffix in {".ini", ".setools"}:
            for prior in lines[:index]:
                header = re.match(r"^\s*\[([^]]+)\]\s*(?:[#;].*)?$", prior)
                if header:
                    section = header.group(1).strip()
            if suffix == ".ini" and not section:
                section = "DEFAULT"
        return Site(path, start, end, "config", section=section,
                    scope="paragraph")
    return None


def _parse_file_tags(path, root):
    relative = Path(path).relative_to(root).as_posix()
    source = Path(path).read_text(encoding="utf-8")
    lines = source.splitlines()
    suffix = Path(path).suffix.lower()
    tags, errors = [], []
    is_snakemake = suffix == ".smk" or Path(path).name == "Snakefile"
    tree = None
    consumed = set()
    comment_tokens = {}
    if suffix == ".py":
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(source, filename=relative)
        except SyntaxError as error:
            return [], [f"{relative}: cannot parse tagged Python: {error}"]
        found, issues, consumed = _python_docstring_tags(relative, source, tree)
        tags.extend(found)
        errors.extend(issues)
        try:
            comment_tokens = {
                token.start[0]: (token.start[1], token.string)
                for token in tokenize.generate_tokens(StringIO(source).readline)
                if token.type == tokenize.COMMENT
            }
        except tokenize.TokenError as error:
            return tags, errors + [f"{relative}: cannot tokenize Python comments: {error}"]
    if suffix == ".py" or is_snakemake or suffix in _CONFIG_SUFFIXES:
        if suffix == ".py":
            tag_comments = [
                (line_no, column, comment[1:].lstrip())
                for line_no, (column, comment) in comment_tokens.items()
                if comment[1:].lstrip().startswith("@sc")
            ]
        else:
            tag_comments = [
                (index + 1, len(line) - len(line.lstrip()), body)
                for index, line in enumerate(lines)
                if (body := _comment_body(line)) is not None and body.startswith("@sc")
            ]
        for line_no, column, body in tag_comments:
            index = line_no - 1
            if suffix == ".py" and line_no in consumed:
                continue
            if suffix == ".py" and column != 0:
                errors.append(f"{relative}:{line_no}: Python statement tags must be module-level comments")
                continue
            meta, ident, error = _tag_line(body)
            if error:
                errors.append(f"{relative}:{line_no}: {error}")
                continue
            site = _comment_site(relative, lines, index, meta, tree,
                                 snakemake=is_snakemake)
            if site is None:
                errors.append(f"{relative}:{line_no}: @sc tag governs no site")
                continue
            prose = _comment_prose(lines, index)
            if ident is not None and not prose:
                errors.append(f"{relative}:{line_no}: local contract {ident} has no prose")
            tags.append(Tag(relative, line_no, tuple(meta.get("decision", ())), ident,
                            {k: v for k, v in meta.items() if k != "decision"},
                            prose, site))
    return tags, errors


def scan_tags(root):
    """Collect all supported-site tags and report malformed/duplicate tags."""

    root = Path(root)
    tags, errors = [], []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or any(part in _SKIP for part in path.parts):
            continue
        if path.name == "CONTRACTS":
            continue
        if path.suffix.lower() not in _CONFIG_SUFFIXES | {".py", ".smk"} and path.name != "Snakefile":
            continue
        found, issues = _parse_file_tags(path, root)
        tags.extend(found)
        errors.extend(issues)
    by_id = {}
    for tag in tags:
        if tag.ident:
            by_id.setdefault(tag.ident, []).append(tag)
    for ident, found in sorted(by_id.items()):
        if len(found) > 1:
            where = ", ".join(f"{tag.path}:{tag.line}" for tag in found)
            errors.append(f"duplicate local-contract id {ident}: {where}")
    return tags, errors


def decision_ids(record, prefix=""):
    """Scoped decision IDs: bare at top level, dotted in sub-analyses."""

    if not isinstance(record, dict):
        return set()
    result = {f"{prefix}{key}" for key in (record.get("decisions") or {})}
    for name, analysis in (record.get("analyses") or {}).items():
        if isinstance(analysis, dict):
            result.update(decision_ids(analysis, f"{prefix}{name}."))
    return result


def tag_errors(tags, record):
    """Check citations and require every ASTRA decision to have a site."""

    known = decision_ids(record)
    cited = set()
    errors = []
    for tag in tags:
        for decision in tag.decisions:
            cited.add(decision)
            if decision not in known:
                errors.append(f"{tag.path}:{tag.line}: @sc cites unknown decision {decision!r}")
    errors.extend(f"{decision}: decision has no tagged site" for decision in sorted(known - cited))
    return errors


def _parse_values(rationale):
    if not isinstance(rationale, str) or "Values:" not in rationale:
        return (), None
    if rationale.count("Values:") != 1:
        return (), "rationale must contain at most one Values: sentence"
    tail = rationale.rsplit("Values:", 1)[1].strip()
    if not tail.endswith("."):
        return (), "Values: sentence must end with a period"
    entries = tuple(part.strip() for part in tail[:-1].split(";"))
    if not entries or any(not entry for entry in entries):
        return (), "Values: sentence contains an empty ref"
    return entries, None


def _split_assertion(reference):
    """Split one ``ref = value`` entry (not a quoted cut's ``==``)."""

    parts = re.split(r"\s+=\s*", reference.strip(), maxsplit=1)
    ref = parts[0]
    expected = parts[1].strip() if len(parts) == 2 else None
    if not ref or re.search(r"\s|=", ref):
        raise ValueError("expected a ref optionally followed by ' = value'")
    if expected is None or not expected or expected.startswith("="):
        raise ValueError("Values entries require a nonempty ' = value'")
    return ref, expected


def _parse_reference(reference):
    if "::" in reference:
        path, selector = reference.split("::", 1)
        return "python", path, selector
    if "#" in reference:
        path, selector = reference.split("#", 1)
        return "config", path, selector
    return "bare", "", reference


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


def _ini_parser(text, *, strict=True):
    parser = configparser.ConfigParser(
        interpolation=None, strict=strict, allow_no_value=True
    )
    parser.optionxform = str.lower
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
_YAML_BOOLEANS = {"true": True, "false": False, "yes": True, "no": False, "on": True, "off": False}
_NO_BOOLEANS = {}


def _booleans_for(suffix):
    if suffix == ".ini":
        return _INI_BOOLEANS
    if suffix in {".sex", ".psfex"}:
        return _ASTROMATIC_BOOLEANS
    if suffix == ".py":
        return _PYTHON_BOOLEANS
    if suffix in {".yaml", ".yml"}:
        return _YAML_BOOLEANS
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


def _iter_decision_defs(document, prefix=""):
    if not isinstance(document, dict):
        return
    for ident, definition in (document.get("decisions") or {}).items():
        yield f"{prefix}{ident}", definition
    for name, analysis in (document.get("analyses") or {}).items():
        if isinstance(analysis, dict):
            yield from _iter_decision_defs(analysis, f"{prefix}{name}.")


def _sites_for(tags, decision):
    unique = {}
    for tag in tags:
        if decision in tag.decisions and tag.site is not None:
            unique[tag.site.identity] = tag.site
    return list(unique.values())


def _path_matches(site_path, qualifier):
    path = Path(qualifier)
    if path.is_absolute() or ".." in path.parts or not qualifier:
        return False
    normalized = path.as_posix().lstrip("./")
    return site_path == normalized or site_path.endswith("/" + normalized)


def _split_config_key(selector, suffix, site):
    if suffix in {".ini", ".setools"} and "." in selector:
        return selector.rsplit(".", 1)
    section = site.section if suffix in {".ini", ".setools"} else ""
    return section, selector.rsplit(".", 1)[-1]


def _strip_config_comment(line, suffix):
    stripped = line.strip()
    if not stripped or stripped.startswith(("#", ";")):
        return ""
    if suffix != ".ini" and "#" in stripped:
        stripped = stripped.split("#", 1)[0].strip()
    return stripped


def _config_values_in_site(root, site, selector):
    """Return this config site's value for selector, or None if it doesn't."""

    target = Path(root) / site.path
    text = target.read_text(encoding="utf-8")
    suffix = target.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        data = yaml.safe_load(text)
        parts = selector.split(".")
        node = data
        for part in parts:
            if not isinstance(node, dict) or part not in node:
                return None
            node = node[part]
        key = parts[-1]
        line_matches = any(
            re.match(rf"^\s*{re.escape(key)}\s*:", line)
            for line in text.splitlines()[site.start - 1:site.end]
        )
        return node if line_matches else None

    lines = text.splitlines()
    section, key = _split_config_key(selector, suffix, site)
    if suffix == ".ini":
        key = key.lower()
    current_section = "DEFAULT" if suffix == ".ini" else ""
    values, predicates = [], []
    for number, raw in enumerate(lines, 1):
        stripped = _strip_config_comment(raw, suffix)
        header = re.match(r"^\[([^]]+)\]$", stripped)
        if header:
            current_section = header.group(1).strip()
            continue
        if number < site.start or number > site.end or not stripped:
            continue
        if suffix in {".ini", ".setools"}:
            wanted_section = section or (site.section if site.scope == "section" else "")
            if wanted_section and current_section != wanted_section:
                continue
        if suffix == ".ini":
            match = re.match(r"^([^:=\s][^:=]*?)\s*[:=]\s*(.*)$", stripped)
            if not match or match.group(1).strip().lower() != key:
                continue
            value = match.group(2).strip()
            base_indent = len(raw) - len(raw.lstrip())
            for continuation in lines[number:site.end]:
                if not continuation.strip():
                    break
                if continuation.lstrip().startswith(("#", ";")):
                    continue
                indent = len(continuation) - len(continuation.lstrip())
                if indent <= base_indent:
                    break
                value += " " + continuation.strip()
            predicate = False
        elif suffix == ".setools":
            match = re.match(rf"^{re.escape(key)}(?=$|\s|=|<|>)(.*)$", stripped)
            if not match:
                continue
            value = match.group(1).strip()
            predicate = value.startswith(("==", "!=", "<", ">"))
            if value.startswith("=") and not predicate:
                value = value[1:].strip()
        elif suffix == ".param":
            match = re.match(rf"^{re.escape(key)}(?:\s*\(([^)]*)\)|\s+(.*))?$", stripped)
            if not match:
                continue
            value = match.group(1) if match.group(1) is not None else (match.group(2) or "")
            predicate = False
        else:
            match = re.match(rf"^{re.escape(key)}(?=$|\s|=|\()(.*)$", stripped)
            if not match:
                continue
            value = match.group(1).strip()
            predicate = False
            if value.startswith("="):
                value = value[1:].strip()
        if not value:
            raise ValueError(f"no active value for {selector!r}")
        values.append(value)
        predicates.append(predicate)
    if not values:
        return None
    if len(values) > 1 and not all(predicates):
        raise ValueError(f"ambiguous active values for {selector!r}: {values!r}")
    return values if len(values) > 1 else values[0]


def _section_exists(text, suffix, section):
    if suffix == ".ini":
        try:
            parser = _ini_parser(text, strict=False)
        except configparser.Error as error:
            raise ValueError(f"cannot parse INI file: {error}") from error
        return section == parser.default_section or parser.has_section(section)
    return any(line.strip() == f"[{section}]" for line in text.splitlines())


def _yaml_key_lines(lines, selector):
    """Return lines for a simple nested YAML mapping selector."""

    wanted = selector.split(".")
    stack, found = [], []
    for number, raw in enumerate(lines, 1):
        match = re.match(r"^(\s*)([^:#][^:]*?)\s*:\s*(?:.*)?$", raw)
        if not match:
            continue
        indent = len(match.group(1))
        while stack and stack[-1][0] >= indent:
            stack.pop()
        key = match.group(2).strip().strip("'\"")
        stack.append((indent, key))
        if [part for _, part in stack] == wanted:
            found.append(number)
    return found


def _active_config_lines(root, path, selector, site=None):
    """Find active line settings for a key throughout one config file."""

    target = Path(root) / path
    suffix = target.suffix.lower()
    lines = target.read_text(encoding="utf-8").splitlines()
    if suffix in {".yaml", ".yml"}:
        return _yaml_key_lines(lines, selector)

    file_site = site or Site(path, 1, len(lines), "config", scope="file")
    section, key = _split_config_key(selector, suffix, file_site)
    if suffix == ".ini":
        key = key.lower()
    current_section = "DEFAULT" if suffix == ".ini" else ""
    found = []
    for number, raw in enumerate(lines, 1):
        stripped = _strip_config_comment(raw, suffix)
        if not stripped:
            continue
        header = re.match(r"^\[([^]]+)\]$", stripped)
        if header:
            current_section = header.group(1).strip()
            continue
        if suffix in {".ini", ".setools"} and section:
            if current_section != section:
                continue
        if suffix == ".ini":
            match = re.match(r"^([^:=\s][^:=]*?)\s*(?:[:=]\s*(.*))?$", stripped)
            active_key = match.group(1).strip().lower() if match else None
        elif suffix == ".setools":
            match = re.match(rf"^{re.escape(key)}(?=$|\s|=|<|>)(.*)$", stripped)
            active_key = key if match else None
        elif suffix == ".param":
            match = re.match(
                rf"^{re.escape(key)}(?:\s*\(([^)]*)\)|\s+(.*))?$", stripped
            )
            active_key = key if match else None
        else:
            match = re.match(rf"^{re.escape(key)}(?=$|\s|=|\()(.*)$", stripped)
            active_key = key if match else None
        if active_key == key:
            found.append(number)
    return found


def _config_absent_actual(root, path, selector):
    """Return whether a config setting is active anywhere in its target scope."""

    target = Path(root) / path
    suffix = target.suffix.lower()
    text = target.read_text(encoding="utf-8")
    file_site = Site(path, 1, len(text.splitlines()), "config", scope="file")
    section, key = _split_config_key(selector, suffix, file_site)
    if suffix == ".ini":
        key = key.lower()
        try:
            parser = _ini_parser(text, strict=False)
        except configparser.Error as error:
            raise ValueError(f"cannot parse INI file: {error}") from error
        if section and not _section_exists(text, suffix, section):
            return None
        sections = [section] if section else [parser.default_section, *parser.sections()]
        active_sections = [
            name for name in sections if parser.has_option(name, key)
        ]
        if section and section != parser.default_section:
            explicit = parser._sections.get(section, {})
            if key not in explicit and parser.has_option(parser.default_section, key):
                return f"active setting inherited from [{parser.default_section}]"
        return (
            f"active setting in section(s) {active_sections!r}"
            if active_sections else _NO_SETTING
        )
    if suffix == ".setools" and section and not _section_exists(
        text, suffix, section
    ):
        return None
    active = _active_config_lines(root, path, selector)
    return f"active setting on line(s) {active!r}" if active else _NO_SETTING


def _check_absent_entry(root, decision, ref, kind, qualifier, selector, tags):
    """Check absence against a same-decision tag in the file or section."""

    if kind == "python":
        return (
            f"{decision}: ref {ref!r}: expected no active setting, actual "
            "<unreadable> (absence refs apply only to config files)"
        )
    matching_paths = set()
    for tag in tags:
        if decision not in tag.decisions or tag.site is None:
            continue
        site = tag.site
        if qualifier and not _path_matches(site.path, qualifier):
            continue
        suffix = Path(site.path).suffix.lower()
        if suffix not in _CONFIG_SUFFIXES:
            continue
        section, _ = _split_config_key(selector, suffix, site)
        if suffix in {".ini", ".setools"} and section:
            if site.scope != "file" and site.section != section:
                continue
            try:
                text = (Path(root) / site.path).read_text(encoding="utf-8")
            except OSError:
                continue
            if not _section_exists(text, suffix, section):
                continue
        matching_paths.add(site.path)

    if len(matching_paths) != 1:
        actual = (
            "<unreadable>" if not matching_paths
            else f"<ambiguous: {len(matching_paths)} tagged files>"
        )
        return (
            f"{decision}: ref {ref!r}: expected no active setting, actual "
            f"{actual} (absence requires a same-decision tag in the file or "
            f"named section; found {len(matching_paths)})"
        )

    path, = matching_paths
    try:
        actual = _config_absent_actual(root, path, selector)
    except (ValueError, OSError, configparser.Error, yaml.YAMLError) as error:
        return (
            f"{decision}: ref {ref!r}: expected no active setting, "
            f"actual <unreadable>: {error}"
        )
    if actual is None:
        section, _ = _split_config_key(
            selector, Path(path).suffix.lower(),
            Site(path, 1, 1, "config", scope="file"),
        )
        return (
            f"{decision}: ref {ref!r}: expected no active setting, actual "
            f"<unreadable> (section {section!r} does not exist)"
        )
    if actual is _NO_SETTING:
        return None
    return f"{decision}: ref {ref!r}: expected no active setting, actual {actual}"


def _python_value_in_site(root, site, selector):
    source = (Path(root) / site.path).read_text(encoding="utf-8")
    tree = _python_tree(source)
    try:
        symbol, _ = _code_selector(selector)
    except ValueError:
        return None
    direct = bool(
        site.symbol
        and (symbol == site.symbol or symbol.startswith(site.symbol + "."))
    )
    full_selector = (
        selector if direct else f"{site.symbol}.{selector}"
        if site.symbol else selector
    )
    try:
        node = _selected_python_node(tree, full_selector)
    except ValueError as error:
        message = str(error)
        if "needs one binding" in message:
            found = re.search(r"found (\d+)", message)
            if found and int(found.group(1)) > 1:
                raise ValueError(f"ambiguous binding: {message}") from error
            return None
        if not direct:
            return None
        if "missing" in message:
            return None
        raise
    try:
        return ast.literal_eval(node)
    except (ValueError, TypeError) as error:
        raise ValueError("selected Python value is not a literal") from error


def _site_actual(root, site, kind, selector):
    suffix = Path(site.path).suffix.lower()
    if kind == "python":
        if suffix != ".py":
            return None
        return _python_value_in_site(root, site, selector)
    if suffix not in _CONFIG_SUFFIXES:
        return None
    return _config_values_in_site(root, site, selector)


def _value_comparison(expected, actual, suffix, selector):
    booleans = _booleans_for(suffix)
    want = _normalise_value(expected, booleans)
    got = _normalise_value(actual, booleans)
    if suffix == ".ini":
        got = _ini_bool(want, got)
    if (suffix, selector.rsplit(".", 1)[-1]) in {
        (".param", "VIGNET"), (".psfex", "PSF_SIZE")
    }:
        want, got = _square_stamp(want), _square_stamp(got)
    return want == got


def _check_value_entry(root, decision, reference, tags):
    try:
        ref, expected = _split_assertion(reference)
        kind, qualifier, selector = _parse_reference(ref)
        if kind == "bare":
            if not selector or re.search(r"\s", selector):
                raise ValueError("invalid unqualified ref")
        elif not qualifier or Path(qualifier).is_absolute() or ".." in Path(qualifier).parts:
            raise ValueError("path qualifier must be a relative path suffix")
    except ValueError as error:
        return f"{decision}: ref {reference!r}: expected <valid ref>, actual <unreadable>: {error}"

    if expected == ABSENT:
        return _check_absent_entry(
            root, decision, ref, kind, qualifier, selector, tags
        )

    candidates = []
    for site in _sites_for(tags, decision):
        suffix = Path(site.path).suffix.lower()
        expected_kind = "python" if suffix == ".py" else "config"
        if kind != "bare" and kind != expected_kind:
            continue
        if qualifier and not _path_matches(site.path, qualifier):
            continue
        try:
            actual = _site_actual(root, site, expected_kind, selector)
        except (ValueError, OSError, SyntaxError, configparser.Error, yaml.YAMLError) as error:
            return f"{decision}: ref {ref!r}: expected {expected!r}, actual <unreadable>: {error}"
        if actual is not None:
            if expected_kind == "config":
                try:
                    active_lines = _active_config_lines(
                        root, site.path, selector, site
                    )
                except (ValueError, OSError, configparser.Error, yaml.YAMLError) as error:
                    return f"{decision}: ref {ref!r}: expected {expected!r}, actual <unreadable>: {error}"
                outside = [
                    number for number in active_lines
                    if not site.start <= number <= site.end
                ]
                if outside:
                    return (
                        f"{decision}: ref {ref!r}: expected {expected!r}, "
                        f"actual {actual!r}; duplicate active setting outside "
                        f"the governed paragraph on line(s) {outside!r}"
                    )
            candidates.append((site, actual))

    if len(candidates) != 1:
        actual = "<unreadable>" if not candidates else f"<ambiguous: {len(candidates)} tagged sites>"
        return f"{decision}: ref {ref!r}: expected {expected!r}, actual {actual} (ref must resolve to exactly one tagged site; found {len(candidates)})"

    site, actual = candidates[0]
    if expected == ABSENT:
        if actual is _NO_SETTING:
            return None
        return f"{decision}: ref {ref!r}: expected no active setting, actual {actual!r}"
    try:
        suffix = Path(site.path).suffix.lower()
        if _value_comparison(expected, actual, suffix, selector):
            return None
        return f"{decision}: ref {ref!r}: expected {expected!r}, actual {actual!r}"
    except (ValueError, OSError, SyntaxError, configparser.Error, yaml.YAMLError) as error:
        return f"{decision}: ref {ref!r}: expected {expected!r}, actual {actual!r}: {error}"


def value_errors(root, record, tags=None):
    """Check every Values entry against exactly one site tagged for its decision."""

    root = Path(root)
    if tags is None:
        tags, _ = scan_tags(root)
    errors = []
    for decision, definition in _iter_decision_defs(record):
        rationale = definition.get("rationale") if isinstance(definition, dict) else None
        if not isinstance(rationale, str):
            continue
        entries, problem = _parse_values(rationale)
        if problem:
            errors.append(f"{decision}: {problem}")
            continue
        for entry in entries:
            try:
                _split_assertion(entry)
            except ValueError as error:
                errors.append(f"{decision}: Values ref {entry!r}: {error}")
                continue
            result = _check_value_entry(root, decision, entry, tags)
            if result:
                errors.append(result)
    return errors


def _is_decision_marker_call(node):
    function = node.func
    return (
        isinstance(function, ast.Attribute) and function.attr == "decision"
        and isinstance(function.value, ast.Attribute) and function.value.attr == "mark"
        and isinstance(function.value.value, ast.Name)
        and function.value.value.id == "pytest"
    )


def decision_markers(root):
    root = Path(root)
    markers, errors = [], []
    test_root = root / "tests"
    if not test_root.is_dir():
        return markers, errors
    for path in sorted(test_root.rglob("*.py")):
        if any(part in _SKIP for part in path.parts):
            continue
        relative = path.relative_to(root).as_posix()
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=relative)
        except (OSError, SyntaxError) as error:
            errors.append(f"{relative}: cannot parse decision markers: {error}")
            continue
        calls = sorted(
            (node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and _is_decision_marker_call(node)),
            key=lambda node: (node.lineno, node.col_offset),
        )
        for call in calls:
            where = f"{relative}:{call.lineno}"
            if not call.args:
                errors.append(f"{where}: decision marker needs literal string ids")
            if call.keywords:
                errors.append(f"{where}: decision markers accept positional ids only")
            for argument in call.args:
                if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                    markers.append(DecisionMarker(argument.value, relative, call.lineno))
                else:
                    errors.append(f"{where}: decision marker ids must be literal strings")
    return markers, errors


def decision_marker_errors(markers, record):
    known = decision_ids(record)
    return [f"{marker.path}:{marker.line}: decision marker cites unknown decision {marker.decision!r}"
            for marker in markers if marker.decision not in known]


def _rationale_sentence(definition):
    text = str(definition.get("rationale", "")).strip()
    text = re.sub(r"\s+", " ", text)
    parts = re.split(r"(?<=[.!?])\s+", text, maxsplit=1)
    return parts[0] if parts else ""


def _value_sentence(definition):
    entries, problem = _parse_values(definition.get("rationale", ""))
    return "Values: " + "; ".join(entries) + "." if entries and not problem else ""


def repository_errors(root):
    """Run the bidirectional tag/value and marker checks for a repository."""

    root = Path(root)
    record = load_yaml(root / "astra.yaml")
    tags, errors = scan_tags(root)
    markers, marker_parse_errors = decision_markers(root)
    errors = list(errors)
    errors.extend(tag_errors(tags, record))
    errors.extend(value_errors(root, record, tags))
    errors.extend(marker_parse_errors)
    errors.extend(decision_marker_errors(markers, record))
    universe_path = root / "universes" / "committed.yaml"
    if universe_path.is_file():
        errors.extend(universe_errors(record, load_yaml(universe_path)))
    return errors


_FORBID = re.compile(r"^\s*forbid:\s*(\S+)\s*->\s*(\S+)\s*$")
_CONTRACT_TAG = re.compile(r"^\s*@(sc|cc)(?:\s+\[[^]]*\])?\s+([\w][\w.-]*)\s*$")


def forbid_rules(contracts_file):
    """Return ``(contract_id, source, target)`` import-boundary rules."""

    rules, ident = [], None
    for line in Path(contracts_file).read_text(encoding="utf-8").splitlines():
        match = _CONTRACT_TAG.match(line)
        if match:
            ident = match.group(2)
            continue
        match = _FORBID.match(line)
        if match and ident:
            rules.append((ident, *match.groups()))
    return rules


def module_matches(name, pattern):
    return fnmatch.fnmatchcase(name, pattern) or (
        pattern.endswith(".*") and name == pattern[:-2]
    )


def module_name(path, src_root):
    parts = list(Path(path).relative_to(src_root).with_suffix("").parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def imported_names(path, module):
    """Yield ``(line, dotted_name)`` for imports in a Python module."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)
        tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    package = module.split(".")
    if Path(path).name != "__init__.py":
        package = package[:-1]
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield node.lineno, alias.name
        elif isinstance(node, ast.ImportFrom):
            base = node.module or ""
            if node.level:
                parent = package[:len(package) - (node.level - 1)]
                base = ".".join([*parent, *([base] if base else [])])
            yield node.lineno, base
            for alias in node.names:
                if alias.name != "*":
                    yield node.lineno, f"{base}.{alias.name}"


def import_violations(src_root, rules):
    """Find imports forbidden by ``CONTRACTS`` rules."""

    src_root = Path(src_root)
    found = []
    for path in sorted(src_root.rglob("*.py")):
        if any(part in _SKIP for part in path.parts):
            continue
        module = module_name(path, src_root)
        for ident, source, target in rules:
            if not module_matches(module, source):
                continue
            for line, name in imported_names(path, module):
                if module_matches(name, target):
                    found.append(
                        f"{path.relative_to(src_root)}:{line}: {module} imports "
                        f"{name} (contract {ident})"
                    )
    return found


def _decision_description(record, decision):
    for ident, definition in _iter_decision_defs(record):
        if ident == decision:
            label = definition.get("label", decision) if isinstance(definition, dict) else decision
            return label, _rationale_sentence(definition), _value_sentence(definition)
    return decision, "", ""


def _location_path(root, value):
    raw = Path(value)
    if raw.is_absolute():
        try:
            return raw.resolve().relative_to(Path(root).resolve()).as_posix()
        except ValueError:
            return None
    return raw.as_posix().lstrip("./")


def _print_location(root, record, tags, location):
    line = None
    match = re.match(r"^(.*):(\d+)$", location)
    if match:
        location, line = match.group(1), int(match.group(2))
    relative = _location_path(root, location)
    if relative is None:
        print(f"{location}: outside repository")
        return
    matches = []
    for tag in tags:
        site = tag.site
        if site is None or site.path != relative:
            continue
        if line is not None and not (site.start <= line <= site.end or tag.line == line):
            continue
        matches.append(tag)
    by_decision = sorted({decision for tag in matches for decision in tag.decisions})
    print(f"{relative}" + (f":{line}" if line is not None else ""))
    for decision in by_decision:
        label, rationale, values = _decision_description(record, decision)
        print(f"  {decision} — {label}")
        if rationale:
            print(f"    {rationale}")
        if values:
            print(f"    {values}")
    for tag in matches:
        if tag.ident:
            label = tag.meta.get("label", [""])[0]
            suffix = f" [{label}]" if label else ""
            print(f"  Local contract {tag.ident}{suffix}: {tag.prose}")
    if not matches:
        print("  No tagged site governs this location.")


def _print_decision_sites(record, tags, decision):
    label, rationale, values = _decision_description(record, decision)
    print(f"{decision} — {label}")
    if rationale:
        print(f"  {rationale}")
    if values:
        print(f"  {values}")
    sites = _sites_for(tags, decision)
    if not sites:
        print("  No tagged sites.")
    for site in sorted(sites, key=lambda item: (item.path, item.start)):
        print(f"  {site.path}:{site.start}-{site.end} ({site.kind})")
        for tag in tags:
            if decision in tag.decisions and tag.site and tag.site.identity == site.identity and tag.ident:
                print(f"    Local contract {tag.ident}: {tag.prose}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("location", nargs="?", help="path[:line] to inspect")
    parser.add_argument("--decision", help="list all sites for one decision")
    parser.add_argument("--root", default=Path(__file__).resolve().parents[2])
    args = parser.parse_args(argv)
    root = Path(args.root)
    record = load_yaml(root / "astra.yaml")
    tags, errors = scan_tags(root)
    if args.decision:
        _print_decision_sites(record, tags, args.decision)
    elif args.location:
        _print_location(root, record, tags, args.location)
    else:
        parser.error("supply a path[:line] or --decision <id>")
    if errors:
        print(f"\n{len(errors)} tag parse error(s):")
        for error in errors:
            print(f"  {error}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
