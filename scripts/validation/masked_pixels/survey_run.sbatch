#!/bin/bash
#SBATCH --account=def-waerbeke
#SBATCH --job-name=veto-survey
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/cdaley/veto-survey/reproduce-%j.log
# Submit from /scratch/cdaley/wt-veto. The targeted job lists reproduce
# the measured follow-ups, not an adaptive search silently changing its grid.
set -euo pipefail
module load apptainer/1.4.5
REPO=/scratch/cdaley/wt-veto
OUT=${OUT:-/scratch/cdaley/veto-survey/reproduction}
SIF=/project/def-mjhudson/cdaley/containers/shapepipe-develop-240b37e4.sif
SCRIPTS=$REPO/scripts/validation/masked_pixels
mkdir -p "$OUT/tmp" "$OUT/mpl"
cd "$REPO"
export APPTAINERENV_PYTHONPATH=$REPO/src:$REPO:$SCRIPTS
export APPTAINERENV_PYTHONDONTWRITEBYTECODE=1
export APPTAINERENV_OPENBLAS_NUM_THREADS=1 APPTAINERENV_OMP_NUM_THREADS=1 APPTAINERENV_MKL_NUM_THREADS=1
export APPTAINERENV_TMPDIR=$OUT/tmp APPTAINERENV_MPLCONFIGDIR=$OUT/mpl
run() { apptainer exec "$SIF" "$@"; }
{
    date -Is
    git rev-parse HEAD
    git status --short
    echo "SIF=$SIF allocation=$SLURM_JOB_ID"
    run python -c 'import importlib.metadata as m; print({p:m.version(p) for p in ("numpy","scipy","galsim","ngmix")})'
    sha256sum "$SCRIPTS"/survey_*.py
} > "$OUT/provenance.txt"
run python "$SCRIPTS/survey_prevalence.py" --out "$OUT" --nproc 12
run python "$SCRIPTS/survey_fixtures.py" --out "$OUT/fixtures.json" --nproc 14
for suffix in core interp stability; do
    nproc=14
    [[ $suffix == stability ]] && nproc=2
    run python "$SCRIPTS/survey_fixtures.py" --out "$OUT/fixtures-$suffix.json" \
        --nproc "$nproc" --jobs "$SCRIPTS/survey_results/${suffix}_jobs.json"
done
run python "$SCRIPTS/survey_combine.py" --root "$OUT"
run python "$SCRIPTS/survey_report.py" --root "$OUT"
run python -m pytest tests/unit/test_veto_survey.py tests/unit/test_decisions.py \
    -o addopts='' -q -p no:cacheprovider > "$OUT/tests.txt" 2>&1
run python "$SCRIPTS/survey_mutations.py" > "$OUT/mutations.txt" 2>&1
