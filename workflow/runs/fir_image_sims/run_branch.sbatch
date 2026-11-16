#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --output=/scratch/cdaley/imsim/logs/run2/%x-%j.out
# One shear branch of the Fir image simulations through ShapePipe, in one batch
# job (sp run --in-allocation: snakemake's local executor on this node).
# From the repository root:
#
#   sbatch -c 64 --mem=384G -t 16:00:00 -J im2-sp-1p2z \
#       workflow/runs/fir_image_sims/run_branch.sbatch 1p2z_grid_2
#
# A failure is recorded in /scratch/cdaley/imsim/export-grid2/FAILED.
set -euo pipefail
trap 'echo "$(date) ${SLURM_JOB_NAME} ${SLURM_JOB_ID} failed (sp run $1)" >> /scratch/cdaley/imsim/export-grid2/FAILED' ERR
mkdir -p /scratch/cdaley/imsim/export-grid2
module load apptainer
echo "start $(date) $1"
env SP_SNAKEMAKE_ENV=/scratch/cdaley/imsim/env/snakemake SP_PROFILE=fir \
    workflow/bin/sp run --in-allocation -c workflow/runs/fir_image_sims/$1.yaml
echo "done $(date)"
