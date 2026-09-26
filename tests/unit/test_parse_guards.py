"""Parse-time refusals in the Snakefile for campaigns that would otherwise fail
late or destroy products.

Each guard is a top-level ``def`` in the Snakefile, called once at parse time.
The tests lift the ``def`` by name (as ``test_campaign_lineage.py`` does for the
path helpers), exercise it against a stub ``WorkflowError``, and check that the
Snakefile calls it on the live values.
"""

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"


class WorkflowError(Exception):
    pass


def _lift(name):
    text = SNAKEFILE.read_text()
    m = re.search(rf"^def {name}\(.*?(?=^\S)", text, re.M | re.S)
    assert m, f"Snakefile no longer defines {name}()"
    ns = {"Path": Path, "WorkflowError": WorkflowError}
    exec(m.group(0), ns)
    return ns[name]


def _called(call):
    """The Snakefile calls ``call`` at top level (not only inside a def)."""
    return re.search(rf"^{re.escape(call)}", SNAKEFILE.read_text(), re.M)


# --- MCCD: persistence and the star-catalogue merge read PSFEx products only --

def test_mccd_is_refused_naming_the_psfex_only_readers():
    guard = _lift("refuse_unpersistable_psf")
    with pytest.raises(WorkflowError) as exc:
        guard("mccd")
    msg = str(exc.value)
    assert "persist_exp.py" in msg and "merge_star_cat.py" in msg
    assert "PSFEx" in msg


@pytest.mark.parametrize("model", ["psfex", "fake"])
def test_psfex_and_fake_pass(model):
    _lift("refuse_unpersistable_psf")(model)


def test_psf_guard_runs_on_the_parsed_model():
    assert _called("refuse_unpersistable_psf(PSF_MODEL)")


# --- one root: the cleaners would reclaim the products ------------------------

@pytest.mark.parametrize("clean,clean_tiles",
                         [(True, False), (False, True), (True, True)])
def test_one_root_with_a_cleaner_is_refused(tmp_path, clean, clean_tiles):
    guard = _lift("refuse_one_root_cleaners")
    with pytest.raises(WorkflowError) as exc:
        guard(tmp_path, tmp_path, clean, clean_tiles)
    msg = str(exc.value)
    assert "clean_tile" in msg and "clean_exposure" in msg
    assert "clean: false" in msg and "clean_tiles: false" in msg


def test_one_root_is_compared_resolved(tmp_path):
    (tmp_path / "run").mkdir()
    (tmp_path / "alias").symlink_to(tmp_path / "run")
    with pytest.raises(WorkflowError):
        _lift("refuse_one_root_cleaners")(
            tmp_path / "alias", tmp_path / "run", True, False)


def test_one_root_without_cleaners_passes(tmp_path):
    _lift("refuse_one_root_cleaners")(tmp_path, tmp_path, False, False)


def test_two_roots_with_cleaners_pass(tmp_path):
    _lift("refuse_one_root_cleaners")(
        tmp_path / "products", tmp_path / "run", True, True)


def test_root_guard_runs_on_the_raw_flags():
    """The raw config flags, not CLEAN/CLEAN_TILES: those are off outside the
    compute phase, and the prepare parse should already refuse."""
    assert _called('refuse_one_root_cleaners(PRODUCTS_DIR, RUN_DIR, '
                   'flag(config.get("clean")),')
