"""Config lineage of the campaign name and the campaign product paths.

Enforces workflow/CONTRACTS ``campaign-name-is-run``: the campaign's name has
one source, the run config's ``run:``, and every campaign product path is
rooted in ``PRODUCTS_DIR``.

The skill-level lineage checker scans ``.py`` only, and the name is bound in
the Snakefile, so this reads the workflow's sources directly. The path helpers
are lifted out of the Snakefile by name and evaluated against sentinel roots,
which tests what they RETURN rather than how they are spelled.
"""

import re
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / "workflow"
SNAKEFILE = WORKFLOW / "Snakefile"
RULES = sorted((WORKFLOW / "rules").glob("*.smk"))
SCRIPTS = sorted((WORKFLOW / "scripts").glob("*.py"))
LAUNCHER = WORKFLOW / "bin" / "sp"

SNAKEMAKE_SOURCES = [SNAKEFILE, *RULES]
ALL_SOURCES = [*SNAKEMAKE_SOURCES, *SCRIPTS, LAUNCHER]

PRODUCTS = "/sentinel/products"
SCRATCH = "/sentinel/scratch"
CAMPAIGN = "campaign-sentinel"

# Paths that hold the campaign's durable products, and whether each one
# carries the campaign name. Called with a tile/exposure ID where they take one.
PRODUCT_HELPERS = {
    "final_cat": (("210.282",), False),
    "final_cat_hdf5": ((), True),
    "full_starcat": ((), True),
    "prod_exp_dir": (("2605805",), False),
    "prod_exp_manifest": (("2605805", "exp_persist"), False),
    "prod_exp_tar": (("2605805",), False),
}
PRODUCT_TEMPLATES = ("PROD_TILE_DIR", "PROD_EXP_DIR")
COVERAGE_TEMPLATES = {
    "COVERAGE_DIR": f"{PRODUCTS}/coverage",
    "COVERAGE_HSP": f"{PRODUCTS}/coverage/coverage_{CAMPAIGN}.hsp",
    "COVERAGE_MANIFEST": f"{PRODUCTS}/coverage/manifests/coverage_map.json",
}


def _code_lines(path):
    """Source lines with comment-only lines dropped."""
    return [(n, line) for n, line in
            enumerate(path.read_text().splitlines(), 1)
            if not line.lstrip().startswith("#")]


def _snakefile_def(name):
    """The text of one top-level ``def`` in the Snakefile."""
    text = SNAKEFILE.read_text()
    m = re.search(rf"^def {name}\(.*?(?=^\S)", text, re.M | re.S)
    assert m, f"Snakefile no longer defines {name}(); update PRODUCT_HELPERS"
    return m.group(0)


def _assignment(path, name):
    """The expression bound to a top-level workflow constant."""
    m = re.search(rf"^{name}\s*=\s*(.+)$", path.read_text(), re.M)
    assert m, f"{path.name} must assign {name}"
    return m.group(1)


@pytest.fixture(scope="module")
def helpers():
    """The Snakefile's product-path helpers, bound to sentinel roots."""
    ns = {"Path": Path, "PRODUCTS_DIR": Path(PRODUCTS),
          "RUN_DIR": Path(SCRATCH), "CAMPAIGN": CAMPAIGN}
    for name in PRODUCT_HELPERS:
        exec(_snakefile_def(name), ns)
    for name in PRODUCT_TEMPLATES:
        ns[name] = eval(_assignment(SNAKEFILE, name), ns)
    for name in COVERAGE_TEMPLATES:
        ns[name] = eval(_assignment(WORKFLOW / "rules" / "coverage.smk",
                                    name), ns)
    return ns


def test_campaign_is_bound_once_from_run():
    """``CAMPAIGN`` has exactly one binding, and it reads ``config["run"]``."""
    bindings = [(p.name, n, line.strip()) for p in SNAKEMAKE_SOURCES
                for n, line in _code_lines(p)
                if re.match(r"\s*CAMPAIGN\s*=", line)]
    assert [b[2] for b in bindings] == ['CAMPAIGN = config["run"]'], bindings


def test_no_second_campaign_source():
    """No source reads a ``campaign`` config key or names a campaign after a
    root directory; either could disagree with ``run:``."""
    forbidden = re.compile(
        r"""\[\s*['"]campaign['"]\s*\]"""
        r"""|\.get\(\s*['"]campaign['"]"""
        r"|\b(?:PRODUCTS_DIR|RUN_DIR|products_dir)\b[\w)\]]*"
        r"(?:\.parent)*\.(?:name|stem)\b")
    hits = [f"{p.relative_to(REPO_ROOT)}:{n}: {line.strip()}"
            for p in ALL_SOURCES for n, line in _code_lines(p)
            if forbidden.search(line)]
    assert not hits, "second campaign-name source:\n" + "\n".join(hits)


def test_merge_rules_pass_campaign_from_campaign():
    """Every rule param named ``campaign`` is ``CAMPAIGN``, and every
    ``--campaign`` on a shell line is that param."""
    params, flags = [], []
    for p in RULES:
        for n, line in _code_lines(p):
            m = re.match(r"\s*campaign\s*=\s*(.+?),?\s*$", line)
            if m:
                params.append((p.name, n, m.group(1)))
            for arg in re.findall(r"--campaign\s+(\S+)", line):
                flags.append((p.name, n, arg))
    assert params, "no rule passes a campaign param; update this test"
    assert all(v == "CAMPAIGN" for *_, v in params), params
    assert flags and all(a.strip("'\"") == "{params.campaign}"
                         for *_, a in flags), flags


@pytest.mark.parametrize("name", sorted(PRODUCT_HELPERS))
def test_product_paths_are_rooted_in_products_dir(helpers, name):
    args, named = PRODUCT_HELPERS[name]
    path = str(helpers[name](*args))
    assert path.startswith(PRODUCTS + "/"), path
    assert (CAMPAIGN in path) == named, path


@pytest.mark.parametrize("name", PRODUCT_TEMPLATES)
def test_product_templates_are_rooted_in_products_dir(helpers, name):
    assert str(helpers[name]).startswith(PRODUCTS + "/"), helpers[name]


@pytest.mark.parametrize("name,expected", COVERAGE_TEMPLATES.items())
def test_coverage_paths_use_products_dir_and_run(helpers, name, expected):
    """The portable map carries the run name; its manifest stays beside it."""
    assert str(helpers[name]) == expected


def _machine_outputs():
    config = yaml.safe_load((WORKFLOW / "config.yaml").read_text())
    for machine, entry in (config.get("machines") or {}).items():
        for input_type, defaults in entry.items():
            if isinstance(defaults, dict) and "outputs" in defaults:
                yield f"{machine}.{input_type}", defaults["outputs"]


@pytest.mark.parametrize("where,outputs", list(_machine_outputs()))
def test_machine_defaults_name_the_campaign_by_run(where, outputs):
    """Where a machine default sets the persistent root, the campaign in it
    is ``$run``, and the index lives beneath it."""
    products = outputs.get("products_dir")
    if products is None:
        pytest.skip(f"{where} sets no products_dir")
    assert products.rstrip("/").endswith("/$run"), products
    index = outputs.get("index_db")
    if index is not None:
        assert index.startswith(products.rstrip("/") + "/"), (index, products)
