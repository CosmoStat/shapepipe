#!/bin/bash
# run_all.sh
# Run all centroid bias validation tests.
# One-time setup steps are skipped if already done.

path_to_val_cen=$HOME/astro/repositories/github/shapepipe/scripts/validation/centroid
path_to_envs=$path_to_val_cen/envs_config
script=$path_to_val_cen/run_bias_test.sh

RUNDIR="$(pwd)"
SP_WORKTREE=/automnt/n17data/mkilbing/astro/repositories/github/shapepipe_test_centroid_bug
SP_MAIN=/automnt/n17data/mkilbing/astro/repositories/github/shapepipe

export SP_ENV_DIR="${RUNDIR}/envs"

# ---------------------------------------------------------------------------
# One-time setup (each step is skipped if already done)
# ---------------------------------------------------------------------------

## Worktree for test_centroid_bug branch
if [[ ! -d "${SP_WORKTREE}" ]]; then
    echo "Creating worktree for test_centroid_bug..."
    git -C "${SP_MAIN}" worktree add "${SP_WORKTREE}" test_centroid_bug
else
    echo "Worktree already exists: ${SP_WORKTREE}"
fi

## Conda envs
for cfg in centroid_bug centroid_fix centroid_v2; do
    env_name="sp_${cfg}"
    if [[ ! -d "${SP_ENV_DIR}/${env_name}" ]]; then
        echo "Creating conda env ${env_name}..."
        CONDARC=~/.condarc mamba env create -y \
            -f "${path_to_envs}/${cfg}.yml" \
            --prefix "${SP_ENV_DIR}/${env_name}"
    else
        echo "Conda env already exists: ${env_name}"
    fi
done

## Install shapepipe into each env (editable, no deps)
declare -A ENV_SP_PATH=(
    [sp_centroid_bug]="${SP_WORKTREE}"
    [sp_centroid_fix]="${SP_WORKTREE}"
    [sp_centroid_v2]="${SP_MAIN}"
)
for env_name in sp_centroid_bug sp_centroid_fix sp_centroid_v2; do
    pip_bin="${SP_ENV_DIR}/${env_name}/bin/pip"
    sp_path="${ENV_SP_PATH[$env_name]}"
    if [[ -f "${pip_bin}" ]]; then
        if ! "${pip_bin}" show shapepipe &>/dev/null; then
            echo "Installing shapepipe into ${env_name}..."
            "${pip_bin}" install -e "${sp_path}" --no-deps -q
        else
            echo "shapepipe already installed in ${env_name}"
        fi
    fi
done

# ---------------------------------------------------------------------------
# Run tests
# ---------------------------------------------------------------------------

cd "${RUNDIR}"
"${script}" centroid_bug
"${script}" centroid_fix
"${script}" centroid_v2
