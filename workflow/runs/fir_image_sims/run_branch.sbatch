#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --output=/scratch/cdaley/imsim/logs-october/%x-%j.out
# One shear branch of the Fir image simulations through ShapePipe, in one batch
# job: `sp run` under snakemake's local executor, with the run directory, the
# snakemake state and the tile store on the node's NVMe ($SLURM_TMPDIR), so
# /scratch holds only the products. A job that ends early starts its
# unfinished tiles again. From the repository root:
#
#   sbatch -c 32 --mem=192G -t 8:00:00 -J oct-sp7690-sp-1z2z \
#       workflow/runs/fir_image_sims/run_branch.sbatch 1z2z_grid_2-sp7690
#
# A failure is recorded in /scratch/cdaley/imsim/logs-october/FAILED, and the
# run's logs are kept beside this job's output.
set -euo pipefail
LOGS=/scratch/cdaley/imsim/logs-october
RUN=$SLURM_TMPDIR/run
trap 'echo "$(date) ${SLURM_JOB_NAME} ${SLURM_JOB_ID} failed (sp run $1)" >> $LOGS/FAILED
      (cd $RUN && find . -name logs -prune -print0 -o -path "*/output/*.log" -print0 |
         tar czf $LOGS/${SLURM_JOB_NAME}-${SLURM_JOB_ID}-logs.tgz --null -T -) || true' ERR
module load apptainer
echo "start $(date) $1"
CONFIG=$SLURM_TMPDIR/$1.yaml
cat workflow/runs/fir_image_sims/$1.yaml - > $CONFIG <<EOF
outputs:
  run_dir: $RUN
tile_store_root: $SLURM_TMPDIR/tile-store
EOF
MEM_MB=$(( SLURM_MEM_PER_NODE * 95 / 100 ))
env SP_SNAKEMAKE_ENV=/scratch/cdaley/imsim/env/snakemake SP_PROFILE=fir \
    workflow/bin/sp run -c $CONFIG \
    --executor local --cores $SLURM_CPUS_PER_TASK --resources mem_mb=$MEM_MB
echo "done $(date)"
