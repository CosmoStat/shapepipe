"""``workflow/scripts/run_config.py``: the resolver's layering (config.yaml <
input_types < machines < run config, with $variables expanded across every
layer), `run:` as a required key, and every path expanding fully or being
reported.

`run:` names the campaign's merged catalogues (``final_cat_<run>.hdf5``,
``full_starcat_<run>.hdf5``), so a run config whose paths never mention
``$run`` is refused too, and the shipped ``config.yaml`` leaves the name to
the run config.

Container-free: run_config.py needs only PyYAML.
"""

import importlib.util
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_YAML = REPO_ROOT / "workflow" / "config.yaml"

_spec = importlib.util.spec_from_file_location(
    "_run_config", REPO_ROOT / "workflow" / "scripts" / "run_config.py")
run_config = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(run_config)

BASE = {
    "input_type": "data",
    "run": "r1",
    "input_types": {
        "data": {"psf_model": "psfex", "tile_detection": "unions_catalogue"},
        "image_sims": {"psf_model": "fake", "tile_detection": "sextractor"},
    },
    "machines": {
        "m": {
            "base_dir": "/base",
            "data": {
                "tile_list": "$base_dir/$run/tiles.txt",
                "inputs": {"tiles": "$base_dir/tiles",
                           "exposures": "$base_dir/exp"},
                "outputs": {"run_dir": "/scratch/$run",
                            "index_db": "/idx/$run.sqlite"},
            },
            "image_sims": {"psf_model": "mccd",
                           "inputs": {"tiles": "/sims/$run"}},
        },
    },
}


def _resolve(tmp_path, monkeypatch, over, base=BASE):
    monkeypatch.setenv("SP_PROFILE", "m")
    base_path = tmp_path / "config.yaml"
    base_path.write_text(yaml.safe_dump(base))
    run_path = tmp_path / "run.yaml"
    run_path.write_text(yaml.safe_dump(over))
    return run_config.load(str(base_path), str(run_path))


def test_input_type_table_supplies_its_defaults(tmp_path, monkeypatch):
    config = _resolve(tmp_path, monkeypatch, {})
    assert config["psf_model"] == "psfex"
    assert config["tile_detection"] == "unions_catalogue"


def test_machine_entry_beats_input_type_table(tmp_path, monkeypatch):
    config = _resolve(tmp_path, monkeypatch, {"input_type": "image_sims"})
    assert config["psf_model"] == "mccd"
    assert config["tile_detection"] == "sextractor"


def test_run_config_beats_both_tables(tmp_path, monkeypatch):
    config = _resolve(tmp_path, monkeypatch, {
        "input_type": "image_sims", "psf_model": "psfex",
        "tile_detection": "unions_catalogue"})
    assert config["psf_model"] == "psfex"
    assert config["tile_detection"] == "unions_catalogue"


def test_dicts_merge_per_subkey_and_expand(tmp_path, monkeypatch):
    config = _resolve(tmp_path, monkeypatch, {
        "run": "r2", "inputs": {"exposures": "$base_dir/own_exp"}})
    assert config["inputs"] == {"tiles": "/base/tiles",
                                "exposures": "/base/own_exp"}
    assert config["tile_list"] == "/base/r2/tiles.txt"
    assert run_config.unresolved(config) == []


def test_input_type_table_values_expand(tmp_path, monkeypatch):
    base = dict(BASE, input_types={"data": {"psf_dict": "$base_dir/$run.pkl"}})
    config = _resolve(tmp_path, monkeypatch, {}, base=base)
    assert config["psf_dict"] == "/base/r1.pkl"


def test_unexpanded_variable_is_unresolved(tmp_path, monkeypatch):
    config = _resolve(tmp_path, monkeypatch,
                      {"outputs": {"run_dir": "$nowhere/x"}})
    assert run_config.unresolved(config) == ["outputs.run_dir"]


def test_committed_top_level_sets_no_table_key():
    """config.yaml's top level sits beneath both tables only because it sets
    none of their keys: one that did would shadow them, since the tables fill
    only what is unset."""
    config = yaml.safe_load(CONFIG_YAML.read_text())
    table_keys = set()
    for entry in config["input_types"].values():
        table_keys |= set(entry)
    for machine in config["machines"].values():
        for input_type in config["input_types"]:
            table_keys |= set(machine.get(input_type) or {})
    assert not table_keys & set(config)


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
    cfg = run_config.apply_defaults(_machine_config("/b/$run"))
    assert cfg["outputs"]["run_dir"] == "/b/smk-g6/run"
    assert run_config.unresolved(cfg) == []


def test_run_config_shorthands_nest_and_run_is_written_back():
    """Top-level scalars are variables, ${name} delimits, nesting resolves."""
    cfg = run_config.apply_defaults(_machine_config(
        "/b", shear="1p2z", grid="grid_2", run="${shear}_${grid}",
        outputs={"run_dir": "/o/$run/scratch"}))
    assert cfg["run"] == "1p2z_grid_2"
    assert cfg["outputs"]["run_dir"] == "/o/1p2z_grid_2/scratch"
    assert run_config.unresolved(cfg) == []


def test_run_holding_an_unknown_variable_is_reported():
    cfg = run_config.apply_defaults(_machine_config("/b", run="$nope"))
    assert "run" in run_config.unresolved(cfg)


def test_optional_path_with_an_unknown_variable_is_reported():
    cfg = run_config.apply_defaults(_machine_config(
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
