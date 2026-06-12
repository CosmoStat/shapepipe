# tests/_artifacts — the publish seam

Guardrail tests emit status artifacts here at run time:

- `<name>.png`  — distribution plot (e.g. R1/R2 with jackknife band)
- `<name>.json` — scalar status (metric values, errors, `PASS`/`FAIL`)
- `<name>.md`   — short human-readable summary

A later GitHub Pages (or other status) step reads this directory and
publishes. **That deploy is intentionally not built here** — this is the emit
side of the seam only. The star shear-response test
(`tests/cluster/test_star_shear_response.py`) is the first producer, via
`tests/helpers/artifacts.py::emit_star_response_artifacts`.

Generated files are git-ignored (see `.gitignore`); only this README and the
ignore rule are tracked, so the seam exists on a clean checkout with nothing
stale committed.
