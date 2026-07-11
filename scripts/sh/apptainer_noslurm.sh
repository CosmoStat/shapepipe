#!/bin/bash
# Run apptainer exec after stripping all SLURM/PMI/PMIX/OMPI environment
# variables that cause OpenMPI to crash when running inside a SLURM job.
# Usage: apptainer_noslurm.sh <sif> [apptainer_args...] -- <cmd> [cmd_args...]

for v in $(env | grep -Eo "^(SLURM|PMI|PMIX|OMPI)[^=]*"); do
  unset "$v"
done

exec apptainer exec "$@"
