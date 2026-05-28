#!/usr/bin/env bash
# Run the bundled example pipeline against a writable copy of /app/example.
#
# Works under both docker (writable /app) and apptainer/SIF (read-only /app):
# the example tree is copied into a fresh tmp workdir so the config's
# relative paths (INPUT_DIR=./example/data/..., OUTPUT_DIR=./example/output)
# resolve to writable locations.
#
# Extra args are forwarded to shapepipe_run, e.g.:
#   shapepipe_run_example -v
set -euo pipefail

WORK="$(mktemp -d -t shapepipe-example-XXXXXX)"
cp -r /app/example "$WORK/"
cd "$WORK"
exec shapepipe_run -c example/config.ini "$@"
