"""The run-config resolver's layering: config.yaml < input_types < machines <
run config, with $variables expanded across every layer.

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
