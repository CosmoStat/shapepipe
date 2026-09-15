#!/usr/bin/env python3
"""Resolve a run config: config.yaml, then SP_RUN_CONFIG merged on top, then
the `machines:` entry for (machine, input_type) filling whatever is still
unset. One definition shared by the Snakefile, bin/sp and container.py.

CLI (used by bin/sp): run_config.py CONFIG_YAML RUN_CONFIG KEY[.SUBKEY]
prints the resolved value, or an empty line if unset. RUN_CONFIG may be "".
"""

import sys

import yaml

PLACEHOLDER = "TBD"
MACHINE_KEYS = ("tile_list", "retrieve", "container", "inputs", "outputs")
REQUIRED = ("tile_list", "inputs.tiles", "inputs.exposures",
            "outputs.run_dir", "outputs.index_db")


def merge(base, over):
    """Recursive dict merge; `over` wins (as snakemake's update_config)."""
    out = dict(base)
    for key, value in over.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], value)
        else:
            out[key] = value
    return out


def _expand(value, base_dir):
    if isinstance(value, str):
        return value.replace("$base_dir", base_dir)
    if isinstance(value, dict):
        return {k: _expand(v, base_dir) for k, v in value.items()}
    return value


def apply_machine_defaults(config):
    """Fill unset MACHINE_KEYS from machines[machine][input_type], in place.

    A key already in `config` wins; for `inputs`/`outputs` the merge is per
    sub-key. `$base_dir` expands to machines[machine].base_dir.
    """
    entry = (config.get("machines") or {}).get(config.get("machine")) or {}
    defaults = entry.get(config.get("input_type", "data")) or {}
    base_dir = str(entry.get("base_dir", ""))
    for key in MACHINE_KEYS:
        if key not in defaults:
            continue
        value = _expand(defaults[key], base_dir)
        if isinstance(value, dict):
            config[key] = merge(value, config.get(key) or {})
        else:
            config.setdefault(key, value)
    return config


def get(config, dotted):
    value = config
    for part in dotted.split("."):
        value = value.get(part) if isinstance(value, dict) else None
    return value


def unresolved(config):
    """REQUIRED keys that are unset or still the placeholder."""
    return [k for k in REQUIRED if get(config, k) in (None, "", PLACEHOLDER)]


def load(config_yaml, run_config=None):
    with open(config_yaml) as f:
        config = yaml.safe_load(f) or {}
    if run_config:
        with open(run_config) as f:
            config = merge(config, yaml.safe_load(f) or {})
    return apply_machine_defaults(config)


if __name__ == "__main__":
    config_yaml, run_config, key = sys.argv[1:4]
    value = get(load(config_yaml, run_config or None), key)
    print("" if value is None else value)
