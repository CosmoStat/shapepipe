# Testing

ShapePipe has two complementary ways to check that things work: the automated
**test suite** (run by CI, and by you during development) and a **smoke run** of
the bundled example pipeline.

## The automated test suite

The test suite runs with [pytest](https://docs.pytest.org/) under a single
discovery root, `tests/`. See `tests/README.md` for the tiers, markers, and
machine-specific requirements.

`.github/workflows/deploy-image.yml` builds the image and runs tests inside it
on branch pushes and pull requests targeting the integration or release branches.
See that workflow for the suite, science checks, workflow-DAG checks, and binary
and example smoke runs.

### Running it yourself

The supported way is inside the image (it carries the `test` extra and the
astromatic binaries), exactly as CI does:

```bash
docker run --rm -e HYPOTHESIS_PROFILE=ci \
    ghcr.io/cosmostat/shapepipe:develop pytest -rX
```

From a development checkout (e.g. inside an Apptainer sandbox with the repo
bind-mounted), run pytest directly:

```bash
pytest                       # the whole suite, with coverage reported
pytest tests/unit            # just the fast structural tier
```

The `ci` Hypothesis profile is deterministic (fixed seed, capped examples) so CI
is fast and reproducible; running without it locally explores more widely.
Coverage is reported (`--cov-report=term-missing`) but not gated.

## Smoke test — the example pipeline

To confirm a fresh install runs end to end, the
[example](https://github.com/CosmoStat/shapepipe/tree/develop/example) directory
contains a `config.ini` that chains a few dummy modules:

- `python_example_runner` — a module written entirely in Python
- `serial_example_runner` — a module run in serial mode
- `execute_example_runner` — a module that calls a system executable

None do anything scientifically interesting, but running them confirms ShapePipe
is up:

```bash
shapepipe_run -c ./example/config.ini
```

Output is saved to `./example/output`. Inside the container, the
`shapepipe_run_example` wrapper does the same against a writable copy of the
example tree, so it works even on a read-only (apptainer/SIF) filesystem.
