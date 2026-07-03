# ShapePipe test suite

One discovery root, tiered by subdirectory. Everything lives under `tests/` and
is driven by `pytest` from the repo root (in the dev container — see the project
`CLAUDE.md`); `pyproject.toml` `[tool.pytest.ini_options]` carries the config
(`testpaths = ["tests"]`).

## Where tests live

| Location | Holds | Why here |
|----------|-------|----------|
| `tests/module/` | **module-unit tests** — the fitter, file handler, split-exp, vignetmaker, ngmix internals, the GalSim weight-validation suite | per-module unit/property/integration tests; import package internals directly. (Relocated from `src/shapepipe/tests/` so the suite has one home.) |
| `tests/unit/` | **structural tests** — every submodule imports, configs parse, shell scripts lint, runner metadata is well-formed, console entry points respond to `-h` | suite-level checks on the *tree*, not any one module |
| `tests/science/` | **fast scientific guardrails** — controlled simulations with a known answer, runnable in the inner loop with nothing from the cluster | scientific correctness that must stay green on every commit |
| `tests/cluster/` | **candide guardrails** — read real on-disk catalogs / submit cluster jobs | need the cluster + real data; marked and auto-skipped off it |
| `tests/helpers/` | shared, non-test library code (cluster submission, artifact emission, the star-response R-function) | imported by tests as `tests.helpers.*`; not collected as tests |
| `tests/_artifacts/` | plots + status JSON/markdown emitted by guardrail tests | the seam a later GitHub Pages step publishes from |

`tests/` is a Python package (each tier has an `__init__.py`), so the shared
helpers import as `tests.helpers.*` from any tier. A bare `pytest` discovers the
whole tree from the single `testpaths` root. `conftest.py` at the repo root is
the single source of markers, environment detection, and the candide skip
policy — it applies everywhere.

## Markers

```
slow      heavy compute (minutes); excluded from the fast inner loop
candide   needs the candide cluster and/or its real data; auto-skipped elsewhere
```

`--strict-markers` is on, so a typo'd marker is an error, not a silent no-op.

A `candide`-marked test is **collected everywhere** (so `--collect-only` shows
it exists) but **skipped off-cluster** with a clear reason. Candide is detected
by hostname (`c0x` login, `nXX` compute nodes) via `on_candide()` in
`conftest.py`; override with `SHAPEPIPE_ON_CANDIDE=0/1`.

## Running

```bash
pytest                                  # full suite (cluster tests skip off-candide)
pytest -m "not slow"                    # fast inner loop
pytest tests/science                    # just the fast scientific guardrails
pytest tests/science/test_mbias.py      # the m-bias guardrail alone

# cluster guardrails (on candide, or forced):
SHAPEPIPE_ON_CANDIDE=1 pytest tests/cluster
```

## The guardrail tests

### m-bias (`tests/science/test_mbias.py`) — fast, local

Injects a known shear `g1 = 0.02` through `make_ngmix_observation` /
`do_ngmix_metacal`, builds the per-object metacal response
`R11 = (g1_1p - g1_1m) / (2·step)`, and asserts the response-corrected
multiplicative bias `|m| < 5e-3` in the ideal (low-noise) limit. The controlled
path recovers `g1 ≈ 0.01996` vs 0.02 injected (`m ≈ 2e-3`); a broken response
correction or dropped deconvolution blows `|m|` past tolerance in seconds.

### Star shear-response (`tests/cluster/test_star_shear_response.py`) — slow + candide

Faithful reproduction of Fabian's R-function (handoff 2026-06-12). Reads ngmix
metacal catalogs at
`output_*/run_sp_tile_ngmix_Ng1u_*/ngmix_runner/output/ngmix-*.fits`
(HDU map `[2]=1p [1]=1m [4]=2p [3]=2m`), masks `20 < mag < 26`, computes
`R1 = (g1_1p − g1_1m)/0.02`, `R2 = (g2_2p − g2_2m)/0.02` with a 100-group
jackknife error, and asserts `<R1>`, `<R2>` within `0 ± 0.03` — deconvolution
should leave stars no net shear response.

The outputs `base_dir` is parametrized (`SHAPEPIPE_STAR_GRID_OUTPUTS`, default
`/home/hervas/n25/SP_simu_fab/SP_1z2z_star_grid/outputs/`). By default it
evaluates an **existing** outputs dir (seconds). Regenerating those outputs is
a multi-hour SLURM job chain
(`tile_launcher → job_per_tile_newversion → job_sp_14`, through the container)
— wired via `tests/helpers/cluster.py` but **opt-in only**
(`SHAPEPIPE_REGENERATE_STAR_GRID=1` plus `SHAPEPIPE_STAR_GRID_TILES`). The test
never launches a heavy job on its own.

Both R1/R2 assertions share a module-scoped fixture that also **emits
artifacts** (plot + JSON + markdown) to `tests/_artifacts/` — on pass *and*
fail — so the status surface always reflects the latest run.

## Helpers

- `tests/helpers/star_response.py` — Fabian's R-function as a library
  (`compute_star_response`, `find_ngmix_files`, `tile_responses`,
  `jackknife_error`).
- `tests/helpers/cluster.py` — `srun` / `sbatch` / `submit_star_grid_chain` /
  `wait_for_jobs`, the SLURM submission seam. Import-safe off-cluster.
- `tests/helpers/artifacts.py` — `emit_star_response_artifacts`, the GitHub
  Pages publish seam (emit only; the deploy is intentionally elsewhere).
