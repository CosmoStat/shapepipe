#!/bin/bash
# run_bias_test.sh <config_name>
#
# Run a centroid-bias test under the conda environment and shapepipe branch
# defined in configs/<config_name>.yaml.
#
# The shapepipe repo root is derived automatically from this script's
# location (scripts/validation/centroid/ -> root).
#
# Results are written to ${SP_RESULTS_DIR}/<config_name>/ where
# SP_RESULTS_DIR defaults to ~/astro/Runs/shapepipe/CFIS/centroid_bug/results
# and can be overridden by setting the environment variable.
#
# Conda envs are looked for at ${SCRIPT_DIR}/envs/<conda_env>.
# They are not committed; create them with:
#   mamba env create -f envs/<config_name>.yml --prefix envs/<conda_env>
# i.e.:
#  mamba env create -f envs/centroid_bug.yml --prefix envs/sp_centroid_bug
#  mamba env create -f envs/centroid_fix.yml --prefix envs/sp_centroid_fix
#
# Usage:
#   ./run_bias_test.sh centroid_bug
#   ./run_bias_test.sh centroid_fix

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Repo root is 3 levels up from scripts/validation/centroid/
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

CONFIG_NAME="${1:?Usage: $0 <config_name>}"
CONFIG_FILE="${SCRIPT_DIR}/configs/${CONFIG_NAME}.yaml"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: config not found: $CONFIG_FILE"
    exit 1
fi

# Parse YAML config with Python (avoids yq dependency)
eval "$(python3 - <<EOF
import yaml, shlex
with open('${CONFIG_FILE}') as f:
    c = yaml.safe_load(f)
for k, v in c.items():
    print(f"{k.upper()}={shlex.quote(str(v))}")
EOF
)"

# shapepipe_path defaults to repo root if not set in config
SHAPEPIPE_PATH="${SHAPEPIPE_PATH:-${REPO_ROOT}}"

# Results directory: env var override > config value > default
DEFAULT_RESULTS="${HOME}/astro/Runs/shapepipe/CFIS/centroid_bug/results/${CONFIG_NAME}"
RESULTS_DIR="${SP_RESULTS_DIR:-${RESULTS_DIR:-${DEFAULT_RESULTS}}}"

echo "=== Bias test: ${NAME} ==="
echo "    ${DESCRIPTION}"
echo "    conda env   : ${CONDA_ENV}"
echo "    branch      : ${SHAPEPIPE_BRANCH}"
echo "    shapepipe   : ${SHAPEPIPE_PATH}"
echo "    test script : ${TEST_SCRIPT}"
echo "    results     : ${RESULTS_DIR}"
echo

# Verify branch
CURRENT_BRANCH=$(git -C "${SHAPEPIPE_PATH}" branch --show-current 2>/dev/null || echo "unknown")
if [[ "${CURRENT_BRANCH}" != "${SHAPEPIPE_BRANCH}" ]]; then
    echo "ERROR: ${SHAPEPIPE_PATH} is on branch '${CURRENT_BRANCH}', expected '${SHAPEPIPE_BRANCH}'."
    echo "       Check out the correct branch or add a git worktree."
    exit 1
fi

mkdir -p "${RESULTS_DIR}"
LOG="${RESULTS_DIR}/run_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to ${LOG}"
echo

# Conda envs live next to this script under envs/<conda_env_name>
ENV_PREFIX="${SCRIPT_DIR}/envs/${CONDA_ENV}"
if [[ ! -d "${ENV_PREFIX}" ]]; then
    echo "ERROR: conda env not found at ${ENV_PREFIX}. Create it with:"
    echo "  mamba env create -f ${SCRIPT_DIR}/envs/${CONFIG_NAME}.yml --prefix ${ENV_PREFIX}"
    exit 1
fi
echo "    env prefix  : ${ENV_PREFIX}"
echo

# PYTHONNOUSERSITE=1 prevents ~/.local site-packages from overriding the
# conda env packages (e.g. a system-wide ngmix install).
mamba run --prefix "${ENV_PREFIX}" \
    env PYTHONNOUSERSITE=1 \
        PYTHONPATH="${SHAPEPIPE_PATH}/src:${PYTHONPATH:-}" \
    python "${SHAPEPIPE_PATH}/${TEST_SCRIPT}" \
    2>&1 | tee "${LOG}"

echo
echo "Done. Results in ${RESULTS_DIR}"
