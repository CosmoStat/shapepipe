"""Every name for a ShapePipe run directory must be the same name.

A ShapePipe stage's run directory is named once, by ``RUN_NAME`` in the stage's
``.ini``, and then referred to FOUR more times by things that must find it:

  * later ``INPUT_DIR`` lines in that same ``.ini``, which chain one module's
    output into the next module's input;
  * ``completeness.STAGE_DIR``, which the Snakefile's ``unit_pre`` uses to
    ``rm -rf`` the stage's run dir before the run, and ``completeness.py`` uses
    to count the products after it;
  * ``persist_exp.RUN_NAME``, which is where ``exp_persist`` looks for the PSF
    products it tars onto /project.

Nothing enforces that agreement at run time, and each way of breaking it fails
LATE and in a different voice. The smk-g7 campaign lost every exposure to
exactly this: renaming ``run_sp_exp_SxSePsfPi`` to ``run_sp_exp_SxSePsf`` (so
that the psfex and mccd chains share one downstream path) updated the config
line that existed at the time and left the other four behind. The first ran as
``ERROR: Invalid INPUT_DIR``; the retry ran as ``ERROR: Directory ... already
exists``, because the ``rm -rf`` had been clearing a directory that no longer
had that name.

So the invariant is asserted here, statically, on the committed files.
"""

import configparser
import importlib.util
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = REPO_ROOT / "workflow" / "scripts"
CONFIG_DIR = REPO_ROOT / "workflow" / "config" / "cfis"

# The `$SP_RUN/output/<run>/<module>/output` shape that a chained INPUT_DIR has.
CHAINED = re.compile(r"\$SP_RUN/output/(run_sp_[A-Za-z0-9_]+)/")

# The exposure PSF stage is the one place two configs must agree with EACH
# OTHER: `psf:` in config.yaml picks between them, and everything downstream
# reads one path precisely because neither name mentions the model.
PSF_CONFIGS = ["config_exp_psfex.ini", "config_exp_mccd.ini"]


def _load(name, *, required=True):
    path = SCRIPTS / f"{name}.py"
    if not path.exists():
        if required:
            raise AssertionError(f"{path} not found; the workflow calls it by path")
        return None
    sys.path.insert(0, str(SCRIPTS))
    try:
        spec = importlib.util.spec_from_file_location(f"_{name}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS))
    return module


def _run_name(config_path):
    parser = configparser.ConfigParser()
    assert parser.read(config_path) == [str(config_path)]
    return parser["DEFAULT"]["RUN_NAME"].strip()


CONFIGS = sorted(CONFIG_DIR.glob("config_*.ini"))


@pytest.mark.parametrize("config_path", CONFIGS, ids=lambda p: p.name)
def test_chained_input_dirs_name_this_config_s_own_run(config_path):
    """An INPUT_DIR naming a run dir must name the run THIS config writes.

    A config's modules chain through its own output. Naming another config's
    run dir is always a stale rename, never a real cross-stage read: a genuine
    one goes through a run dir this config never writes, and there are none.
    """
    parser = configparser.ConfigParser()
    parser.read(config_path)
    if "RUN_NAME" not in parser["DEFAULT"]:
        pytest.skip(f"{config_path.name} sets no RUN_NAME")
    run_name = _run_name(config_path)

    # Run dirs written by an EARLIER stage are legitimately read across configs
    # (a tile config reads run_sp_tile_Fe, say), so this does not forbid every
    # foreign reference — only a reference to a run dir the exposure PSF stage
    # owns, which is the one name that has drifted and the one no config may
    # spell for itself in two ways.
    offenders = []
    for section in parser.sections():
        value = parser[section].get("INPUT_DIR", "")
        for referenced in CHAINED.findall(value):
            # Only judge references to run dirs the PSF stage owns: those are
            # the ones a config must not spell with a stale name.
            if referenced.startswith("run_sp_exp_SxSe") and referenced != run_name:
                offenders.append((section, referenced))

    assert not offenders, (
        f"{config_path.name} declares RUN_NAME = {run_name} but chains "
        f"INPUT_DIR through a different run dir: {offenders}. "
        "ShapePipe fails this as 'Invalid INPUT_DIR' on every unit."
    )


def test_both_psf_configs_write_the_same_run_dir():
    """psfex and mccd share a run-dir name so downstream never branches."""
    names = {name: _run_name(CONFIG_DIR / name) for name in PSF_CONFIGS}
    assert len(set(names.values())) == 1, (
        f"the exposure PSF configs disagree on RUN_NAME: {names}. "
        "Downstream (completeness, persist_exp, clean_exposure) reads one path "
        "for both, so a per-model name silently breaks the model it is not."
    )


def test_stage_dir_and_persist_agree_with_the_configs():
    """completeness.STAGE_DIR and persist_exp.RUN_NAME name the real dir.

    STAGE_DIR is what ``unit_pre`` clears before a run, so a stale entry here
    is the 'Directory already exists' failure on the FIRST retry, after a
    first attempt that looked like something else entirely.
    """
    expected = _run_name(CONFIG_DIR / "config_exp_psfex.ini")

    completeness = _load("completeness")
    level, subdir = completeness.STAGE_DIR["exp_psf"]
    assert level == "exp"
    assert subdir == expected, (
        f"completeness.STAGE_DIR['exp_psf'] is {subdir!r} but the config "
        f"writes {expected!r}; unit_pre would rm -rf the wrong directory and "
        "the retry would die on 'Directory already exists'."
    )

    # persist_exp.py arrives with feat/persist-exp-products; on a branch
    # without it there is simply no third copy of the name to disagree.
    persist = _load("persist_exp", required=False)
    if persist is None:
        return
    assert persist.RUN_NAME == expected, (
        f"persist_exp.RUN_NAME is {persist.RUN_NAME!r} but the config writes "
        f"{expected!r}; exp_persist would tar an empty product set."
    )
