---
name: 'ShapePipe test suite: tiered, in-image, property-based'
status: open
tags:
    - shapepipe
    - ci
    - testing
    - constitution
created-at: 2026-05-29T15:24:03.630476192+02:00
outcome: 'Tier 0-2 coverage is locally green under the deterministic Hypothesis CI profile: PYTHONPATH=src uv run --extra test pytest -rX reports 249 passed, 5 skipped, no XPASS, and Hypothesis loads derandomize=True/max_examples=50 from conftest.py. Docker/in-image verification remains blocked in this worker because docker is not installed.'
shuttle:
    enabled: true
    kind: oneshot
    host: c03
    project_dir: /automnt/n17data/cdaley/unions/shapepipe
    agent: codex
---

## Desired State

ShapePipe has a test suite that a contributor trusts: green means the code
actually works, and the suite grows with the pipeline instead of rotting behind
it. Concretely:

- **One CI path.** Tests run *inside the dev image* via the existing
  `pytest` step in `deploy-image.yml` — the image is the environment, so the
  suite exercises exactly what ships. There is no second test workflow.
- **Tiered coverage** (below), each tier cheap relative to the bugs it catches.
- **Property-based tests** (via `hypothesis`) wherever a function has a clear
  invariant — round-trip, idempotence, involution, shape/units preservation —
  rather than only example-based cases. Run them in CI under a **deterministic
  profile**: a registered `hypothesis` profile with a fixed seed
  (`derandomize=True`) and a capped `max_examples` (~50) so CI is fast and
  reproducible; a looser local/dev profile can explore more widely.
- **Coverage is reported, not gated.** Run with `--cov-report=term-missing` so
  thin spots are visible, but do **not** fail CI on a coverage threshold — the
  suite is too young to hard-gate; revisit once it's mature.
- **An honest xfail registry.** Known-broken cases are `xfail(strict=True)`
  with a reason and an issue link; nothing is xfailed that actually passes, and
  nothing silently skipped that should run. A strict xfail that starts passing
  is a signal to delete the marker, not a failure to suppress.
- **Fast feedback.** Tiers 0–2 need no external binaries or survey data and run
  in seconds; the binary/data-dependent integration tier is `skipif`-gated on
  tool/data availability so it never blocks a PR for an install reason.

Done is not a number. Done is: the analytic surface in Tier 1–2 below is
covered, property tests guard the invariant-bearing functions, the suite is the
single gate on `develop`, and the xfail registry reflects reality.

## The coverage surface

The pipeline is 29 module runners (`src/shapepipe/modules/*_package/*_runner.py`)
over a shared contract, plus pipeline plumbing (`pipeline/`) and utilities. Most
*scientific* modules can only be exercised end-to-end (they shell out to
`source-extractor`/`psfex` or need survey data); but a real analytic surface
underneath is pure Python and currently almost untested. Tier by tier:

**Tier 0 — structural (from #708, reconcile & keep).** Parametrized, one case
per file/module: `bash -n` over `scripts/sh/*`, configparser over
`example/**/*.ini`, import every `shapepipe.*` submodule, `-h` on every
`[project.scripts]` entry, `@module_runner` metadata present on every
`*_runner.py`. These caught real recent bugs (shell syntax in #706, unreachable
`except ImportError`, broken entry points). **Keep #708's five test files and
its `testpaths` fix; drop its `ci-dev.yml`** (folded into the in-image run).

**Tier 1 — pure functions + property tests (new, no I/O).** The analytic core,
unit- and property-tested:
- `utilities/galaxy.py::sigma_to_fwhm` — round-trip inverse, positivity,
  array broadcast (has example tests; add properties).
- `modules/get_images_package/get_images.py::in2out_pattern` — dots→dashes,
  idempotence on already-transformed input (has example tests; add properties).
- `modules/ngmix_package/ngmix.py::MegaCamFlip` — shape preservation; involution
  (apply twice = identity for the rot180 CCDs); CCD-set selection.
- `pipeline/file_io.py::get_unit_from_fits_header` — TCUNI/TTYPE lookup
  correctness; IndexError on missing column.
- `pipeline/config.py::CustomParser.getlist`/`getexpanded` — list parsing,
  `$SP_*` expansion, defaults.

**Tier 2 — geometry & I/O with synthetic FITS (new, no binaries).** Construct
astropy WCS/FITS in-memory; no survey data:
- `modules/vignetmaker_package/vignetmaker.py::convert_pos` — world→pixel→world
  round-trip within tolerance; `_get_stamp` shape `(n, 2r+1, 2r+1)` and pixel
  rounding; `make_mask` replaces `< -1e29` with the mask value, leaves the rest.
- `modules/psfex_interp_package/psfex_interp.py::_get_galaxy_positions` —
  `(n_gal, 2)` shape, KeyError on missing column.
- `pipeline/file_io.py` catalogue classes — FITS write/read round-trip.
- `modules/module_decorator.py` — attribute injection + the validation rules
  (file_pattern/file_ext length match, run_method ∈ {parallel,serial}).

**Tier 3 — integration (external binaries / data, `skipif`-gated).** Extend
`src/shapepipe/tests/test_tpv_external_tools.py` (already gates on
`source-extractor`/`psfex` on PATH). The full module runners (sextractor,
psfex, mccd, ngmix, mask) live here; the example pipeline (`example/config.ini`,
already a CI smoke) is the end-to-end case. Grow this opportunistically — it is
not where the bulk of the value is.

## Context

- **Read first:** #708's five files in `tests/unit/` (the Tier-0 spec, with a
  worked `xfail`-registry pattern), `src/shapepipe/tests/` (existing examples,
  esp. `test_split_exp.py` for the synthetic-FITS style), and
  `modules/module_decorator.py` (the runner contract).
- **The runner contract:** `@module_runner(version, input_module, file_pattern,
  file_ext, depends, executes, numbering_scheme, run_method)` wrapping a function
  of `(input_file_list, run_dirs, file_number_string, config, module_config_sec,
  w_log)` returning `(stdout, stderr)`.
- **CI today:** `deploy-image.yml` builds the dev image and runs `pytest` inside
  it. `pyproject.toml` has `addopts = --verbose --cov=shapepipe`. #708 corrects
  `testpaths` to `["tests", "src/shapepipe/tests"]` — adopt that so the in-image
  run discovers both trees.
- **`hypothesis`** is not yet a dependency — add it to the `test` extra in
  `pyproject.toml` and regenerate `uv.lock`.

## First task — reconcile #708's stale xfails

#708 was written before recent fixes. Several `strict=True` xfails now PASS,
which would fail the suite. Before extending, make the registry honest against
current `develop`:
- `summary_run` `-h` and the `canfar_monitor`/`canfar_*` `IndentationError` —
  **fixed by #714**; remove those xfails (they should pass now).
- `mask`/`ngmix`/`uncompress_fits` import failures (astroquery/numba/fitsio) —
  **resolved by the container migration** (deps now in the image); remove.
- `stile`/`treecorr.corr2` failures (`mccd_plots`, `random_cat`) — likely still
  real; keep, and link the stile-removal thread.
Run the suite in the image, see what XPASSes, and clean accordingly.

## Skills

`felt`. Verify empirically — run `pytest` inside the dev image (build `--target
dev`, `docker run … pytest`); the GH `deploy-image` run is ground truth.

## Evidence

```bash
# Build the dev image and run the suite (ground truth for "does it pass in CI")
docker build --target dev -t sp-dev . && docker run --rm sp-dev pytest tests/unit src/shapepipe/tests

# What's covered / where coverage is thin
docker run --rm sp-dev pytest --cov=shapepipe --cov-report=term-missing

# Any strict xfail now XPASSing? (registry is dishonest if so)
docker run --rm sp-dev pytest -rX tests/unit
```

Done = the analytic surface (Tier 1–2) is covered, property tests guard the
invariant functions, `deploy-image.yml` is the single gate, and `pytest -rX`
shows no unexpected XPASS.

## Scope

- **In:** `tests/` + `src/shapepipe/tests/`, the in-image CI step, the `test`
  extra in `pyproject.toml`/`uv.lock`, reworking #708 freely.
- **Out:** the conda removal (done), scientific-algorithm changes (Axel's
  centroid fix #725 is Martin's to handle), and a second CI workflow — there is
  one test path.

## Open Questions

- **stile/treecorr** — the `mccd_plots`/`random_cat` import failures come from
  `stile v0.1` importing a `treecorr.corr2` that newer treecorr removed. Is this
  being fixed by the stile-removal work already in flight, or should those
  modules' tests stay xfailed pending it? Check the stile-removal thread; keep
  the xfails honest either way.
