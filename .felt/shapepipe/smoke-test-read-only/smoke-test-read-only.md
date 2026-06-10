---
id: 01KTCHWZZ923H1H5DR6AS3Q5H6
name: Smoke test must work in read-only mode
status: closed
tags:
    - shapepipe
    - docker
    - infra
created-at: 2026-05-28T10:32:25.53742271+02:00
closed-at: 2026-06-10T17:13:38.485047103+02:00
outcome: 'shapepipe_run_example wrapper (copy example to /tmp) landed via PR #731, merged to develop; CI smoke steps use it. Read-only SIF runs work.'
---

## The gap

The recently-landed dc13582e (uv/pytest on read-only filesystems) closed
the `uv run pytest` story but missed the example pipeline smoke test.
Reproducing on the current `:develop` image:

```
$ apptainer exec --pwd /app ./shapepipe-sandbox \
      shapepipe_run -c /app/example/config.ini
ERROR: [Errno 30] Read-only file system: './example/output/shapepipe_runs.txt'
```

The config's `OUTPUT_DIR = ./example/output` is read-relative-to-cwd; cwd
is `/app` (the image WORKDIR); `/app` is read-only under SIF. CI passes
today because docker mounts the container filesystem writable. So the
gap is silent in CI but blocks anyone using the same command under
apptainer.

## The fix

A wrapper script that prepares a writable copy of the example tree in
`/tmp` and runs from there. `/tmp` is tmpfs and is writable in both
read-only SIFs and writable sandboxes.

```bash
#!/usr/bin/env bash
# scripts/sh/shapepipe_run_example.sh
set -euo pipefail
WORK="$(mktemp -d -t shapepipe-example-XXXXXX)"
cp -r /app/example "$WORK/"
cd "$WORK"
exec shapepipe_run -c example/config.ini "$@"
```

Lives under `scripts/sh/` so the Dockerfile's existing auto-symlink rule
makes it `shapepipe_run_example` on `$PATH`. Tested locally by copying
the script into an apptainer sandbox built from the `:develop` image and
running under `apptainer exec` (no `--writable`).

CI smoke steps in `.github/workflows/deploy-image.yml` switch from
`shapepipe_run -c /app/example/config.ini` to `shapepipe_run_example` in
both the runtime and dev target blocks.

## Connections

Sits in the same family as [[shapepipe/docker-multistage]] (which
introduced the runtime/dev split) and [[shapepipe/docker-uv-revert]]
(which moved uv writable targets to `/tmp` via env vars). [[shapepipe/prs-in-flight]]
gets a new "in-flight" entry once the PR is up.
