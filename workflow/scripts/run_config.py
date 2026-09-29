#!/usr/bin/env python3
"""Resolve a run config: config.yaml, then SP_RUN_CONFIG merged on top, then
the `machines:` entry for (machine, input_type) filling whatever is still
unset. One definition shared by the Snakefile, bin/sp and container.py.

CLI (used by bin/sp): run_config.py CONFIG_YAML RUN_CONFIG KEY[.SUBKEY]
prints the resolved value, or an empty line if unset. RUN_CONFIG may be "".
"""

import os
import re
import sys

import yaml

PLACEHOLDER = "TBD"
# Keys the machines: table may default, per (machine, input_type).
# psf_model/psf_dict belong here because they are per-input_type facts,
# not per-run ones: psf_model=fake is only legal with
# input_type=image_sims, and psf_dict is the sim PSF it reads. A key
# also present at the TOP level of config.yaml shadows the table (the
# setdefault below only fires when the key is absent), so a key listed
# here must not carry a top-level default as well.
MACHINE_KEYS = ("tile_list", "retrieve", "container", "inputs", "outputs",
                "psf_model", "psf_dict")
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


def apply_machine_defaults(config):
    """Fill unset MACHINE_KEYS from machines[machine][input_type], in place.

    The machine is `machine:` when the run config states one, else SP_PROFILE
    (default nibi) -- the same value bin/sp picks the SLURM profile with.

    A key already in `config` wins; for `inputs`/`outputs` the merge is per
    sub-key. In all of these, `$base_dir` expands to machines[machine].base_dir
    and `$name` / `${name}` to any top-level scalar (`$run` to `run:`).
    """
    machine = config.get("machine") or os.environ.get("SP_PROFILE", "nibi")
    entry = (config.get("machines") or {}).get(machine) or {}
    defaults = entry.get(config.get("input_type", "data")) or {}
    # Every top-level scalar is a variable, so a run config can define its own
    # shorthands. They are resolved AGAINST EACH OTHER first, to a fixpoint, so
    # one shorthand may be written in terms of another
    # (grid: grid_2 / base: .../${grid} / run: ${shear}_${grid}). One pass is
    # not enough: re.sub does not rescan what it substitutes, so a nested name
    # would survive into the paths and only surface at the unresolved check.
    # The loop is capped, so a self-reference (a: $b, b: $a) stops rather than
    # spinning -- the leftover $ is then caught by unresolved() by design.
    variables = {k: v for k, v in config.items()
                 if isinstance(v, (str, int, float))}
    variables["base_dir"] = entry.get("base_dir")
    for _ in range(10):
        resolved = {k: (_expand(v, variables) if isinstance(v, str) else v)
                    for k, v in variables.items()}
        if resolved == variables:
            break
        variables = resolved
    # `run` is consumed downstream (paths, the hdf5 group name), so the
    # resolved value has to go back into the config, not just the table.
    if isinstance(config.get("run"), str):
        config["run"] = variables.get("run", config["run"])
    for key in MACHINE_KEYS:
        default = defaults.get(key)
        if isinstance(default, dict):
            config[key] = merge(default, config.get(key) or {})
        elif default is not None:
            config.setdefault(key, default)
        if key in config:
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


# Config blocks the workflow no longer reads, and what replaced each. A key
# nothing reads is a switch that silently does nothing (`coverage.enabled:
# false` would leave the exposure-count map on), so retired() makes it an error.
RETIRED = {
    "coverage": "exposure_maps.nexp (on/off: exposure_maps.nexp.enabled; "
                "resolution: the shared exposure_maps.nside / "
                "exposure_maps.nside_coverage)",
    "defect_map": "exposure_maps.defect (oversample: "
                  "exposure_maps.defect.oversample; resolution: the shared "
                  "exposure_maps.nside / exposure_maps.nside_coverage)",
}


def retired(config):
    """Every place a RETIRED key appears — top level, a `machines:` entry, or
    one of its per-input_type blocks — as `(path, replacement)` pairs."""
    found = [(key, RETIRED[key]) for key in RETIRED if key in config]
    for machine, entry in (config.get("machines") or {}).items():
        if not isinstance(entry, dict):
            continue
        found += [(f"machines.{machine}.{key}", RETIRED[key])
                  for key in RETIRED if key in entry]
        for input_type, block in entry.items():
            if isinstance(block, dict):
                found += [(f"machines.{machine}.{input_type}.{key}",
                           RETIRED[key])
                          for key in RETIRED if key in block]
    return found


def unresolved(config):
    """REQUIRED keys that are unset, the placeholder, or hold an unexpanded
    $variable, then every MACHINE_KEYS value (recursively through
    inputs/outputs) that still holds one (e.g. `$run` with no `run:` set, or
    a misspelt name)."""
    missing = [k for k in REQUIRED
               if get(config, k) in (None, "", PLACEHOLDER)
               or "$" in str(get(config, k))]
    dollar = [k for key in MACHINE_KEYS for k in _dollar_keys(config.get(key), key)]
    return missing + [k for k in dollar if k not in missing]


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
