# Container Workflow

ShapePipe ships as a container image. This page covers the two image
targets and the three configuration files (`pyproject.toml`, `uv.lock`,
`Dockerfile`) that determine what's inside.

## Two image targets

The Dockerfile builds two flavours of the image from a shared base:

| Tag | Target | Use case |
|-----|--------|----------|
| `:<branch>` (e.g. `:develop`, `:latest`) | `dev` | Interactive work, sandboxed apptainers, CI test runs. Includes everyday CLI tools (`vim`, `tmux`, `htop`, `rg`, `fd`, `jq`, `bat`, `less`, `git-lfs`, …) and **all** Python extras (test, lint, doc, jupyter, fitsio, release). |
| `:<branch>-runtime` | `runtime` | Canfar batch jobs, downstream `FROM` clauses. Slim — no interactive tools, only the `jupyter` and `fitsio` Python extras on top of core deps. |

Both share the `base` stage (system libraries + uv + lockfile copy), so
the heavy work happens once during the build.

```bash
# Default = dev (everyday image)
apptainer build --sandbox shapepipe docker://ghcr.io/cosmostat/shapepipe:develop

# Slim runtime image for batch jobs
docker run --rm ghcr.io/cosmostat/shapepipe:develop-runtime shapepipe_run -c /app/example/config.ini
```

## Two ways to use the image

Apptainer gives you two filesystem postures. Pick based on what you're
doing — and the image is designed to behave correctly in both.

### Read-only SIF — running, testing, batch

```bash
apptainer pull shapepipe.sif docker://ghcr.io/cosmostat/shapepipe:develop
apptainer shell shapepipe.sif

# inside the container — every path under / is read-only except /tmp,
# $HOME, and a few system mounts. This is the right posture for *running*
# things: pipelines, tests, scripts.
cd /app && uv run pytest         # works — see "How read-only mode works"
shapepipe_run -c /app/example/config.ini
```

The image bakes three env vars (`UV_NO_SYNC=1`, `UV_CACHE_DIR=/tmp/uv-cache`,
`COVERAGE_FILE=/tmp/.coverage`) so `uv run`, the uv cache, and pytest-cov
all redirect their writes to `/tmp` (which is tmpfs, always writable).
You shouldn't need to set anything yourself.

What you cannot do in this mode: install new packages, edit `/app`, run
`uv sync`. For any of that, switch to a writable sandbox.

### Writable sandbox — development, exploration

```bash
# Build once. The sandbox is a directory tree, not a SIF.
apptainer build --sandbox sp-dev docker://ghcr.io/cosmostat/shapepipe:develop

# Bind-mount your host clone of the repo so edits-on-host are visible
# inside; install editable against the mount so it shadows /app.
apptainer shell --writable --bind /path/to/local/shapepipe:/mnt/shapepipe sp-dev
pip install -e /mnt/shapepipe

# Ad-hoc additions during exploration:
uv add some-package        # mutates pyproject.toml + uv.lock inside container
uv pip install foo         # throwaway install, doesn't touch pyproject

# When you decide an addition should stick, copy the mutated pyproject /
# uv.lock back to the host repo, commit, and rebuild the image.
```

The dev image's `/app` filesystem is `chmod -R go+rwX` so non-root users
in a writable sandbox can mutate the venv freely. The default
`UV_NO_SYNC=1` still prevents `uv run` from auto-syncing, which keeps
the venv stable until you explicitly call `uv sync` or `uv add`. Run
`unset UV_NO_SYNC` if you want auto-sync on `uv run` instead.

### How read-only mode works

A SIF is a single immutable file; apptainer mounts it read-only. A
non-`--writable` `apptainer shell` of a sandbox behaves the same way.
Three failure modes appear without the env-var defaults the image bakes
in, all because something tries to write under `/app` or `$HOME`:

| Symptom | Cause | Fix (already baked in) |
|---|---|---|
| `failed to remove .../INSTALLER: Read-only file system` from `uv sync` | uv tries to mutate `/app/.venv` | `UV_NO_SYNC=1` (use `pytest` directly, or `uv run` skips sync) |
| `Failed to initialize cache at $HOME/.cache/uv` | uv cache wants to write under host home | `UV_CACHE_DIR=/tmp/uv-cache` |
| `OSError: Read-only file system: '/app/.coverage'` | pytest-cov erases `/app/.coverage` at startup | `COVERAGE_FILE=/tmp/.coverage` |

If you bypass `/tmp` (e.g. with a custom apptainer profile) you may need
to override these manually.

## Three configuration layers

Three files determine what the image contains. Each has a clear role; the
trick is knowing which one to edit when something changes.

### `pyproject.toml` — what shapepipe needs

Abstract requirements: package metadata, the runtime dependencies of
shapepipe itself, and optional extras (`test`, `jupyter`, `doc`, …). Uses
*minimum* version constraints, not exact pins:

```toml
dependencies = [
    "astropy>=7.0",
    "numpy>=2.0",
    "galsim>=2.8",
    ...
]
```

This is the **single source of truth** for *what kinds of environments
shapepipe is compatible with*. Edit when:

- adding a new dependency (`uv add foo` will edit it for you)
- bumping a minimum after crossing a major version line
- adding/restructuring extras
- changing project metadata (scripts, version, etc.)

### `uv.lock` — which exact versions

Concrete materialization: of all the version combinations that satisfy
`pyproject.toml`, here is the one we've chosen, with hashes for every
package and every transitive dependency. Generated by uv; never
hand-edited.

```bash
uv lock                              # regenerate from pyproject
uv lock --upgrade                    # bump every pin to latest compatible
uv lock --upgrade-package astropy    # bump only one package
uv add 'foo>=1.2'                    # adds to pyproject AND uv.lock
```

Commit `uv.lock` to the repo — that's how reproducibility transfers
between machines. The Dockerfile uses `uv sync --frozen`, which fails if
`pyproject.toml` and `uv.lock` have drifted. That's the discipline
mechanism: a stale lockfile cannot ship.

### `Dockerfile` — the build recipe

System-level dependencies (apt packages, binaries, non-Python tools),
multi-stage layout, and the invocations that turn pyproject + lockfile
into a working venv. Edit when:

- adding a system package (apt) or non-Python tool (`vim`, `htop`, …)
- changing the build process (new target, new stage)
- adjusting which extras get pre-installed in which target

The Dockerfile does **not** duplicate Python deps — those come from
`pyproject.toml` via `uv sync --frozen`.

## Common workflows

| Change | What to edit | Then |
|---|---|---|
| Add a Python dep | `uv add foo` | Commit `pyproject.toml` + `uv.lock`; rebuild image |
| Bump a Python pin | `uv lock --upgrade-package foo` | Commit `uv.lock`; rebuild |
| Refresh all pins | `uv lock --upgrade` | Commit `uv.lock` (and any pyproject minimums you bumped); rebuild |
| Add a system tool | Edit Dockerfile apt block | Rebuild — pyproject/lockfile untouched |
| Promote a sandbox install | `uv add foo` *inside* sandbox; copy `pyproject.toml` + `uv.lock` back to host | Commit; rebuild |

The asymmetry is deliberate: Python deps go through pyproject + lockfile
(reproducible, auditable), system deps go through Dockerfile (Debian's
versioning). Don't `apt install` something that has a Python wheel; don't
`pip install` something Debian packages directly (e.g. `weightwatcher`).

## Why this shape

- **No conda.** Earlier images double-installed packages via conda *and*
  pip, ballooning the image and creating environment-resolution bugs.
  Dropped in favour of plain Python + uv.
- **`uv sync --frozen`** at build time means the image is bit-exactly
  reproducible from a tagged commit, and impossible to ship with a stale
  lockfile.
- **Astromatic binaries from Debian** (`psfex`, `source-extractor`,
  `weightwatcher`) instead of source builds — Debian carries the
  GCC-compatibility patches that the previous Dockerfile had to apply
  inline with `sed`.
- **Two targets** so canfar batch deployments stay slim while interactive
  users get a working environment with `vim`, `pytest`, etc. on the
  default tag.
