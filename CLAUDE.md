# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is ShapePipe

ShapePipe is a galaxy shape measurement pipeline for weak gravitational lensing, developed at CosmoStat (CEA Paris-Saclay). It orchestrates astronomical processing modules (SExtractor, PSFEx, MCCD, ngmix, etc.) into configurable pipeline runs. Source: `src/shapepipe/`.

## Build & Install

```bash
# Development install (editable, with all extras)
pip install -e ".[dev]"

# Full install with external dependencies (Conda + external tools)
./install_shapepipe --develop

# Container (recommended for full pipeline)
# ghcr.io/cosmostat/shapepipe — see Dockerfile
```

Requires Python >= 3.11. Core deps: joblib, modopt, numpy. Many modules need external executables (source-extractor, psfex, etc.) installed separately.

## Running

```bash
# Run pipeline with config file
shapepipe_run -c <config.ini>

# MPI mode
mpiexec -n <cores> shapepipe_run -c <config.ini>
```

Config files are INI format with sections: `[DEFAULT]`, `[EXECUTION]`, `[FILE]`, `[JOB]`, `[WORKER]`, plus per-module sections like `[SEXTRACTOR_RUNNER]`. See `example/config.ini` for annotated reference, `example/cfis/` for production configs.

## Testing

```bash
# Run all tests (pytest configured in pyproject.toml)
pytest

# Run a single test file
pytest src/shapepipe/tests/test_pipeline.py

# Run a specific test
pytest src/shapepipe/tests/test_pipeline.py::ExecuteTestCase::test_execute
```

pytest is configured with `--verbose --cov=shapepipe`, testpaths = `["shapepipe"]`. Tests use `unittest.TestCase` with `numpy.testing`.

## Linting

```bash
black src/shapepipe/
isort src/shapepipe/
```

## Architecture

### Pipeline engine (`pipeline/`)

`ShapePipe` class in `run.py` is the main orchestrator:
1. `args.py` parses CLI args → `config.py` (`CustomParser`) reads INI config
2. `file_handler.py` resolves input/output paths, discovers files by pattern/numbering
3. `dependency_handler.py` checks Python + executable dependencies
4. `job_handler.py` dispatches work via joblib (SMP) or mpi4py (MPI)
5. `worker_handler.py` + `execute.py` run each module, with `timeout.py` support

`file_handler.py` and `file_io.py` are the largest/most complex files — they manage the numbering schemes, inter-module file routing (`last:MODULE`, `all:MODULE` patterns), and output directory structure.

### Module system (`modules/`)

Each processing step is a module with:
- `<name>_runner.py` — decorated with `@module_runner` from `module_decorator.py`
- `<name>_package/` — implementation classes/functions

The `@module_runner` decorator declares: `version`, `input_module`, `file_pattern`, `file_ext`, `depends`, `executes`, `numbering_scheme`, `run_method` (parallel/serial).

Runner signature: `def runner(input_file_list, run_dirs, file_number_string, config, module_config_sec, w_log)` → returns `(stdout, stderr)`.

When a module appears multiple times in `[EXECUTION] MODULE`, its config section gets a `_RUN_N` suffix (e.g., `PYTHON_EXAMPLE_RUNNER_RUN_1`, `PYTHON_EXAMPLE_RUNNER_RUN_2`).

### Utilities (`utilities/`)

- `cfis.py` — CFIS survey-specific helpers (largest utility, ~35KB)
- `summary.py` — pipeline run summary generation
- `file_system.py`, `galaxy.py` — general helpers

### Scripts (`scripts/`)

Python and shell scripts symlinked into the environment during install. Used for catalog creation, merging, and CANFAR cluster operations.

## Style conventions

- Single quotes preferred over double quotes
- Explicit floats: `1.0` not `1.`
- f-strings for formatting
- Line length: 79 characters max; break with `()` not `\`
- Numpydoc-style docstrings
- String concatenation in f-strings: use `+` at start of continuation line
