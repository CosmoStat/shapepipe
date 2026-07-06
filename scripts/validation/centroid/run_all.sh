#!/bin/bash
# run_all.sh
# Run centroid bias validation tests.
# One-time setup steps are skipped if already done.
#
# Usage: run_all.sh [case]
#   case  One of: centroid_bug, centroid_fix, centroid_v2 (default: all three)

ALL_CASES=(centroid_bug centroid_fix centroid_v2)

if [[ $# -gt 0 ]]; then
    case "$1" in
        centroid_bug|centroid_fix|centroid_v2)
            RUN_CASES=("$1") ;;
        *)
            echo "Unknown case '$1'. Valid options: ${ALL_CASES[*]}" >&2
            exit 1 ;;
    esac
else
    RUN_CASES=("${ALL_CASES[@]}")
fi

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
for cfg in "${ALL_CASES[@]}"; do
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
for env_name in sp_centroid_bug sp_centroid_fix; do
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

## sp_centroid_v2: shapepipe + pip-only deps
## pip install is unreliable here (user ~/.local shadows env, requires-python mismatch),
## so we use direct --target installs and a .pth file for shapepipe.
python_v2="${SP_ENV_DIR}/sp_centroid_v2/bin/python"
pip_bin="${SP_ENV_DIR}/sp_centroid_v2/bin/pip"
if [[ -f "${python_v2}" ]]; then
    site_pkg=$(PYTHONNOUSERSITE=1 "${python_v2}" -c "import site; print(site.getsitepackages()[0])")

    # shapepipe: .pth file into env site-packages (bypasses requires-python check)
    pth_file="${site_pkg}/shapepipe.pth"
    if [[ ! -f "${pth_file}" ]]; then
        echo "Adding shapepipe src to sp_centroid_v2 via .pth file..."
        echo "${SP_MAIN}/src" > "${pth_file}"
    else
        echo "shapepipe.pth already present in sp_centroid_v2"
    fi

    # pip-only deps: use --target to install into env site-packages regardless
    # of what is visible in ~/.local (pip --prefix sees user packages as satisfied)
    for pkg in cs_util modopt sqlitedict; do
        import_name="${pkg//-/_}"
        if ! PYTHONNOUSERSITE=1 "${python_v2}" -c "import ${import_name}" 2>/dev/null; then
            echo "Installing ${pkg} into sp_centroid_v2..."
            "${pip_bin}" install "${pkg}" -q --target "${site_pkg}"
        fi
    done
fi

# ---------------------------------------------------------------------------
# Run tests
# ---------------------------------------------------------------------------

cd "${RUNDIR}"
for case in "${RUN_CASES[@]}"; do
    "${script}" "${case}"
done
