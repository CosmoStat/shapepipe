"""Parse and validate the repository's scientific contracts.

A contract is a tagged block colocated with the code it governs: a tag line
(``@sc`` or ``@cc``, an optional bracketed ``key:value`` meta list, then a
stable id) followed by prose up to the next blank line. The line grammar is
the loom ``sc-list`` reference parser's, regex for regex, so a contract that
parses here parses there.

Where contracts live:

* Python: in the docstring of the module, class or function they govern.
  ``sc-list`` reads docstrings only, so a tag line in a ``#`` comment of a
  ``.py`` file is reported as an error rather than silently ignored.
* Snakemake (``.smk``, ``Snakefile``): in ``#`` comment blocks.
* ``CONTRACTS`` files: anywhere in the file; they govern their directory.

A ``decision:<id>`` meta names a decision in ``astra.yaml``: a top-level
decision by its bare id, a sub-analysis decision as ``<analysis>.<id>``.
"""

from dataclasses import dataclass, field
import ast
import io
from pathlib import Path
import re
import tokenize
import warnings

from tests.helpers.astra_record import extract_anchors

TAG = re.compile(r"^\s*@(sc|cc)\b.*$")
VALID = re.compile(r"^\s*@(sc|cc)(?:\s+\[([^\]]*)\])?\s+([\w][\w.-]*)\s*$")
META_PAIR = re.compile(r"[\w.-]+:[^,\s]+")
MISSING_ID = re.compile(r"\s*@(sc|cc)(?:\s+\[[^\]]*\])?\s*")

SCAN_ROOTS = ("src", "workflow", "scripts")
SKIP = {".git", ".venv", "venv", "__pycache__", "node_modules", ".felt"}


@dataclass
class Contract:
    """One parsed contract."""

    tag: str
    id: str
    meta: dict = field(default_factory=dict)
    prose: str = ""
    path: str = ""
    line: int = 0
    scope: str = ""


def parse_block(text, path, offset=0, scope=""):
    """Parse every contract in ``text``; return ``(contracts, errors)``."""

    lines = text.splitlines()
    found, errors = [], []
    for index, line in enumerate(lines):
        if not TAG.match(line):
            continue
        where = f"{path}:{index + 1 + offset}"
        match = VALID.match(line)
        if not match:
            detail = (
                "missing contract id"
                if MISSING_ID.fullmatch(line)
                else "malformed contract line"
            )
            errors.append(f"{where}: {detail}: {line.strip()}")
            continue
        tag, meta_text, ident = match.groups()
        meta = {}
        if meta_text is not None:
            pairs = [part.strip() for part in meta_text.split(",")]
            if any(not META_PAIR.fullmatch(part) for part in pairs):
                errors.append(f"{where}: malformed contract metadata: {meta_text}")
                continue
            meta = dict(part.split(":", 1) for part in pairs)
        prose = []
        for body in lines[index + 1:]:
            if not body.strip():
                break
            prose.append(body.strip())
        found.append(
            Contract(tag, ident, meta, " ".join(prose), str(path),
                     index + 1 + offset, scope)
        )
    return found, errors


def _declarations(tree):
    parents = {
        child: parent
        for parent in ast.walk(tree)
        for child in ast.iter_child_nodes(parent)
    }
    yield tree, "module"
    for node in ast.walk(tree):
        if not isinstance(
            node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            continue
        names, parent = [node.name], parents.get(node)
        while parent is not None and not isinstance(parent, ast.Module):
            if isinstance(
                parent, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
            ):
                names.append(parent.name)
            parent = parents.get(parent)
        yield node, ".".join(reversed(names))


def python_contracts(path, source):
    """Contracts in a Python file's docstrings; tags in comments are errors."""

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse(source, filename=str(path))
    except SyntaxError as error:
        return [], [f"{path}: cannot parse: {error}"]
    found, errors = [], []
    for node, scope in _declarations(tree):
        doc = ast.get_docstring(node, clean=False)
        if not doc:
            continue
        first = 1 if node is tree else node.body[0].lineno
        records, issues = parse_block(doc, path, first - 1, scope)
        found.extend(records)
        errors.extend(issues)
    try:
        tokens = tokenize.generate_tokens(io.StringIO(source).readline)
        for token in tokens:
            if token.type == tokenize.COMMENT and TAG.match(
                token.string.lstrip("#")
            ):
                errors.append(
                    f"{path}:{token.start[0]}: contract in a comment; "
                    "move it into the governing docstring"
                )
    except (tokenize.TokenError, SyntaxError):
        pass
    return found, errors


def snakemake_contracts(path, source):
    """Contracts in a Snakemake file's ``#`` comment blocks."""

    stripped = []
    for line in source.splitlines():
        text = line.strip()
        stripped.append(text[1:] if text.startswith("#") else "")
    return parse_block("\n".join(stripped), path, 0, "file")


def _is_snakemake(path):
    return path.suffix == ".smk" or path.name == "Snakefile"


def contract_files(root, scan_roots=SCAN_ROOTS):
    """Yield every file under ``scan_roots`` that can carry contracts."""

    root = Path(root)
    for top in scan_roots:
        base = root / top
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if any(part in SKIP for part in path.parts) or not path.is_file():
                continue
            if (
                path.name == "CONTRACTS"
                or path.suffix == ".py"
                or _is_snakemake(path)
            ):
                yield path


def collect(root, scan_roots=SCAN_ROOTS):
    """Parse every contract under ``root``; return ``(contracts, errors)``.

    Errors cover malformed tag lines, malformed meta, missing ids and ids
    used more than once. Paths are relative to ``root``.
    """

    root = Path(root)
    contracts, errors = [], []
    for path in contract_files(root, scan_roots):
        relative = path.relative_to(root)
        text = path.read_text(encoding="utf-8")
        if path.name == "CONTRACTS":
            found, issues = parse_block(
                text, relative, 0, str(relative.parent)
            )
        elif _is_snakemake(path):
            found, issues = snakemake_contracts(relative, text)
        else:
            found, issues = python_contracts(relative, text)
        contracts.extend(found)
        errors.extend(issues)
    seen = {}
    for contract in contracts:
        seen.setdefault(contract.id, []).append(contract)
    for ident, items in seen.items():
        if len(items) > 1:
            where = ", ".join(f"{c.path}:{c.line}" for c in items)
            errors.append(f"duplicate contract id {ident}: {where}")
    return contracts, errors


def decision_ids(record, prefix=""):
    """Scoped decision ids: bare at top level, dotted inside sub-analyses."""

    ids = set()
    for decision in (record.get("decisions") or {}):
        ids.add(f"{prefix}{decision}")
    for name, analysis in (record.get("analyses") or {}).items():
        if isinstance(analysis, dict):
            ids |= decision_ids(analysis, f"{prefix}{name}.")
    return ids


def decision_errors(contracts, record):
    """Contracts whose ``decision:`` meta names no decision in the record."""

    known = decision_ids(record)
    return [
        f"{c.path}:{c.line}: contract {c.id} cites unknown decision "
        f"{c.meta['decision']!r}"
        for c in contracts
        if "decision" in c.meta and c.meta["decision"] not in known
    ]


def anchored_refs(record):
    """``(code_symbols, anchored_paths)`` named by the record's anchors.

    ``code_symbols`` holds ``(path, Symbol)`` pairs from ``path::Symbol``
    refs; ``anchored_paths`` holds every path any ref names.
    """

    symbols, paths = set(), set()
    for anchor in extract_anchors(record):
        for reference in anchor.references:
            if "::" in reference:
                path, symbol = reference.split("::", 1)
                symbols.add((path, symbol))
            else:
                path = reference.split("#", 1)[0]
            paths.add(path)
    return symbols, paths


def coverage_report(contracts, record):
    """Report-only gaps between the contracts and the record.

    Returns ``(uncovered, unanchored)``: decision ids no contract cites, and
    ``@sc`` contracts whose declaration is not an anchored symbol. A
    module-docstring contract counts as anchored when an anchor names its
    file or a symbol in it.
    """

    cited = {c.meta["decision"] for c in contracts if "decision" in c.meta}
    uncovered = sorted(decision_ids(record) - cited)
    symbols, paths = anchored_refs(record)
    unanchored = []
    for contract in contracts:
        if contract.tag != "sc":
            continue
        if contract.scope == "module":
            if contract.path in paths:
                continue
        elif (contract.path, contract.scope) in symbols:
            continue
        unanchored.append(contract)
    return uncovered, unanchored
