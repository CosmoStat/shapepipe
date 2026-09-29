"""``workflow/scripts/run_config.py``: `run:` is a required key, and every
machine-key path expands fully or is reported.

`run:` names the campaign's merged catalogues (``final_cat_<run>.hdf5``,
``full_starcat_<run>.hdf5``). ``unresolved()`` already reports a ``$run`` left
unexpanded in a path; these tests pin that a run config whose paths never
mention ``$run`` is refused too, and that the shipped ``config.yaml`` leaves the
name to the run config. A ``$base_dir`` whose value holds ``$run`` expands
through both, and an optional path left holding a ``$`` is reported.
"""

import importlib.util
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "workflow" / "scripts" / "run_config.py"
CONFIG_YAML = REPO_ROOT / "workflow" / "config.yaml"

_spec = importlib.util.spec_from_file_location("_run_config", SCRIPT)
run_config = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_config)

# Every other REQUIRED key, as literal paths with no `$run` in them.
LITERAL = {
    "tile_list": "/data/tiles.txt",
    "inputs": {"tiles": "/data/tiles", "exposures": "/data/exp"},
    "outputs": {"run_dir": "/scratch/run", "index_db": "/data/index.sqlite"},
}


def test_run_is_required():
    assert "run" in run_config.REQUIRED


def test_literal_paths_without_run_are_refused():
    assert run_config.unresolved(dict(LITERAL)) == ["run"]


def test_literal_paths_with_run_resolve():
    assert run_config.unresolved({**LITERAL, "run": "smk-g6"}) == []


def test_shipped_config_leaves_run_to_the_run_config(tmp_path, monkeypatch):
    monkeypatch.setenv("SP_PROFILE", "nibi")
    assert "run" not in (yaml.safe_load(CONFIG_YAML.read_text()) or {})

    no_run = tmp_path / "no_run.yaml"
    no_run.write_text(yaml.safe_dump(LITERAL))
    assert "run" in run_config.unresolved(
        run_config.load(CONFIG_YAML, no_run))

    with_run = tmp_path / "with_run.yaml"
    with_run.write_text(yaml.safe_dump({"run": "smk-test"}))
    cfg = run_config.load(CONFIG_YAML, with_run)
    assert run_config.unresolved(cfg) == []
    assert cfg["outputs"]["run_dir"].endswith("/smk-test")


def _machine_config(base_dir, **over):
    return {"run": "smk-g6", "machine": "m", "input_type": "data",
            "machines": {"m": {"base_dir": base_dir, "data": {
                "tile_list": "$base_dir/tiles.txt",
                "inputs": {"tiles": "$base_dir/tiles",
                           "exposures": "$base_dir/exp"},
                "outputs": {"run_dir": "$base_dir/run",
                            "index_db": "$base_dir/index.sqlite"}}}},
            **over}


def test_base_dir_holding_run_expands_fully():
    cfg = run_config.apply_machine_defaults(_machine_config("/b/$run"))
    assert cfg["outputs"]["run_dir"] == "/b/smk-g6/run"
    assert run_config.unresolved(cfg) == []


def test_run_config_shorthands_nest_and_run_is_written_back():
    """Top-level scalars are variables, ${name} delimits, nesting resolves."""
    cfg = run_config.apply_machine_defaults(_machine_config(
        "/b", shear="1p2z", grid="grid_2", run="${shear}_${grid}",
        outputs={"run_dir": "/o/$run/scratch"}))
    assert cfg["run"] == "1p2z_grid_2"
    assert cfg["outputs"]["run_dir"] == "/o/1p2z_grid_2/scratch"
    assert run_config.unresolved(cfg) == []


def test_run_holding_an_unknown_variable_is_reported():
    cfg = run_config.apply_machine_defaults(_machine_config("/b", run="$nope"))
    assert "run" in run_config.unresolved(cfg)


def test_optional_path_with_an_unknown_variable_is_reported():
    cfg = run_config.apply_machine_defaults(_machine_config(
        "/b", outputs={"products_dir": "/p/$nope/products"},
        inputs={"masks": "/m/$nope"}, container="/c/$nope.sif"))
    assert set(run_config.unresolved(cfg)) == {
        "outputs.products_dir", "inputs.masks", "container"}


def test_shipped_config_carries_no_retired_key():
    assert run_config.retired(yaml.safe_load(CONFIG_YAML.read_text())) == []


def test_retired_keys_are_found_wherever_they_sit():
    """`coverage:` / `defect_map:` at the top level, on a machine, or on a
    machine's input_type block are each reported with their replacement."""
    config = {
        "coverage": {"enabled": False},
        "exposure_maps": {"nexp": {"enabled": True}},
        "machines": {
            "nibi": {"defect_map": {"oversample": 3},
                     "data": {"coverage": {"nside": 131072}}},
            "candide": {"image_sims": {"retrieve": "symlink"}},
        },
    }
    found = dict(run_config.retired(config))
    assert set(found) == {"coverage", "machines.nibi.defect_map",
                          "machines.nibi.data.coverage"}
    assert "exposure_maps.nexp" in found["coverage"]
    assert "exposure_maps.defect" in found["machines.nibi.defect_map"]
    assert "exposure_maps.nside" in found["machines.nibi.data.coverage"]


def test_snakefile_refuses_retired_keys_at_parse_time():
    """The Snakefile raises on retired() before anything else reads config."""
    text = (REPO_ROOT / "workflow" / "Snakefile").read_text()
    check = text.index("run_config.retired(config)")
    assert check < text.index("run_config.apply_machine_defaults(config)")
    assert "raise WorkflowError" in text[check:check + 400]
