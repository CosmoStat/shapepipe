#!/bin/bash
#SBATCH --account=def-waerbeke
#SBATCH --job-name=veto-napoli
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=napoli-scan-%j.log
# Submit from the checkout root; OUT_DIR may name a private evidence directory.
# The cluster chooses the partition. Nothing in this script runs on the login node.
set -euo pipefail
module load apptainer/1.4.5
SIF=/project/def-mjhudson/cdaley/containers/shapepipe-develop-240b37e4.sif
OUT_DIR=${OUT_DIR:-$PWD/napoli-evidence}
mkdir -p "$OUT_DIR"
export APPTAINERENV_PYTHONPATH=$PWD/src:$PWD
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
{
    date -Is
    git rev-parse HEAD
    git status --short
    printf 'SLURM_JOB_ID=%s SIF=%s\n' "$SLURM_JOB_ID" "$SIF"
    sha256sum scripts/validation/masked_pixels/bias_vs_distance.py \
        tests/helpers/defect_response.py src/shapepipe/modules/ngmix_package/ngmix.py
    apptainer exec "$SIF" python -c 'import importlib.metadata as m; print({p: m.version(p) for p in ("numpy", "scipy", "galsim", "ngmix")})'
} > "$OUT_DIR/provenance.txt"
apptainer exec -B "$SLURM_TMPDIR" "$SIF" python \
    scripts/validation/masked_pixels/bias_vs_distance.py napoli \
    "$OUT_DIR/napoli.json" "$SLURM_CPUS_PER_TASK"
apptainer exec "$SIF" python scripts/validation/masked_pixels/bias_vs_distance.py \
    napoli_plot "$OUT_DIR/napoli.png" "$OUT_DIR/napoli.json"
