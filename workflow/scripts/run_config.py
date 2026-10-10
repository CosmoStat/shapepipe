#!/usr/bin/env python3
"""Resolve a run config: config.yaml with SP_RUN_CONFIG merged on top, then
two tables of defaults filling whatever is still unset -- `input_types:`
entry for input_type, overridden by the `machines:` entry for
(machine, input_type). One definition shared by the Snakefile, bin/sp and
container.py.

Precedence, lowest first: config.yaml's top level < input_types[input_type]
< machines[machine][input_type] < the run config. config.yaml's top level
sits beneath both tables only because it sets none of the keys they carry;
tests/unit/test_run_config.py holds it to that.

CLI (used by bin/sp): run_config.py CONFIG_YAML RUN_CONFIG KEY[.SUBKEY]
prints the resolved value, or an empty line if unset. RUN_CONFIG may be "".
"""

import os
import re
import sys

import yaml

PLACEHOLDER = "TBD"
# The tables themselves: read here, never defaults or $-expanded.
TABLES = ("input_types", "machines")
# `run` is required in its own right, not only through `$run` in the paths:
# it names the campaign's merged catalogues, so a run config whose paths
# never mention `$run` must still set it.
REQUIRED = ("run", "tile_list", "inputs.tiles", "inputs.exposures",
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


# $name or ${name}. The brace form is not decoration: `_` is a word character,
# so `$shear_$grid` parses its first name as `shear_` and silently fails to
# match a `shear` key. Anything adjacent to more word characters needs braces.
_VAR_RE = re.compile(r"\$\{(\w+)\}|\$(\w+)")


def _expand(value, variables):
    """Replace $name / ${name} for each set name in `variables`.

    Names come from every top-level scalar in the config plus `base_dir` from
    the machines: entry, so a run config can name its own shorthands -- a long
    output root written once and reused, rather than repeated per path. An
    unknown name is left as-is, which is what makes the "unresolved" check
    below able to see it.
    """
    if isinstance(value, str):
        return _VAR_RE.sub(
            lambda m: str(variables.get(m.group(1) or m.group(2))
                          or m.group(0)),
            value,
        )
    if isinstance(value, dict):
        return {k: _expand(v, variables) for k, v in value.items()}
    return value


def apply_defaults(config):
    """Fill unset keys from the input_types: and machines: tables, in place.

    The machine is `machine:` when the run config states one, else SP_PROFILE
    (default nibi) -- the same value bin/sp picks the SLURM profile with.
    input_types[input_type] holds what follows from the kind of input alone
    (the PSF model, the tile detection); machines[machine][input_type] holds
    where things are on that machine, and wins where the two overlap.

    A key already in `config` wins; for dict values (`inputs`, `outputs`) the
    merge is per sub-key. Every resolved key then has `$name` / `${name}`
    expanded: `$base_dir` is machines[machine].base_dir, and any top-level
    scalar is a variable too (`$run` is the top-level `run:`).
    """
    input_type = config.get("input_type", "data")
    machine = config.get("machine") or os.environ.get("SP_PROFILE", "nibi")
    entry = (config.get("machines") or {}).get(machine) or {}
    defaults = merge((config.get("input_types") or {}).get(input_type) or {},
                     entry.get(input_type) or {})
    # Resolve top-level scalar shorthands against each other before expanding
    # paths: re.sub does not rescan substitutions in one pass. Stop at a
    # fixpoint or ten passes; unresolved() catches remaining $ references,
    # including cycles.
    variables = {k: v for k, v in config.items()
                 if isinstance(v, (str, int, float))}
    variables["base_dir"] = entry.get("base_dir")
    for _ in range(10):
        resolved = {k: (_expand(v, variables) if isinstance(v, str) else v)
                    for k, v in variables.items()}
        if resolved == variables:
            break
        variables = resolved
    for key, default in defaults.items():
        if isinstance(default, dict):
            config[key] = merge(default, config.get(key) or {})
        else:
            config.setdefault(key, default)
    for key in config:
        if key not in TABLES:
            config[key] = _expand(config[key], variables)
    return config


def get(config, dotted):
    value = config
    for part in dotted.split("."):
        value = value.get(part) if isinstance(value, dict) else None
    return value


def _dollar_keys(value, prefix):
    """Dotted keys under `value` whose string still holds a `$`."""
    if isinstance(value, dict):
        return [k for key, sub in value.items()
                for k in _dollar_keys(sub, f"{prefix}.{key}")]
    return [prefix] if isinstance(value, str) and "$" in value else []


def unresolved(config):
    """REQUIRED keys that are unset, the placeholder, or hold an unexpanded
    $variable, then every other resolved value (recursively through
    inputs/outputs; the tables themselves excluded) that still holds one
    (e.g. `$run` with no `run:` set, or a misspelt name)."""
    missing = [k for k in REQUIRED
               if get(config, k) in (None, "", PLACEHOLDER)
               or "$" in str(get(config, k))]
    dollar = [k for key in config if key not in TABLES
              for k in _dollar_keys(config[key], key)]
    return missing + [k for k in dollar if k not in missing]


def catalogue_source(config):
    """`inputs.catalogues` and the retrieve mode its prefix implies.

    The source is a local directory (symlinked in) or a vos: URL (downloaded).
    Returns ("", "symlink") when unset or the placeholder; raises ValueError
    if `tile_detection: unions_catalogue` then has nothing to fetch.
    """
    source = get(config, "inputs.catalogues") or ""
    if source == PLACEHOLDER:
        source = ""
    if config.get("tile_detection") == "unions_catalogue" and not source:
        raise ValueError(
            "tile_detection=unions_catalogue needs `inputs.catalogues:` (a local "
            "directory or vos: URL holding the per-tile CFIS.<tile>.r.cat files).")
    return source, "vos" if source.startswith("vos:") else "symlink"


def load(config_yaml, run_config=None):
    with open(config_yaml) as f:
        config = yaml.safe_load(f) or {}
    if run_config:
        with open(run_config) as f:
            config = merge(config, yaml.safe_load(f) or {})
    return apply_defaults(config)


if __name__ == "__main__":
    config_yaml, run_config, key = sys.argv[1:4]
    value = get(load(config_yaml, run_config or None), key)
    print("" if value is None else value)
