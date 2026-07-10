---
title: 'Snakemake crash: RuntimeError can''t start new thread during pipeline_all'
depends-on:
    - shapepipe
created-at: 2026-07-08T09:32:38.198900721+02:00
outcome: 'Root cause: user thread/process limit exhaustion. ulimit -u soft=1200 hard=1280 on candide. Snakemake''s scheduler Timer thread was just the victim; the real consumers are the N parallel apptainer ShapePipe jobs. apptainer_noslurm.sh strips only SLURM/PMI/PMIX/OMPI vars, so host OMP_NUM_THREADS leaks into every container job: -j $OMP_NUM_THREADS gives N jobs x N BLAS/OpenMP threads each ~ N^2 threads. Fix (agreed, NOT yet applied): in scripts/image_sims_pipeline/Snakefile add --env OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,NUMEXPR_NUM_THREADS=1 to APP_CMD and SPV_CMD, and run with a fixed conservative -j 4or5 instead of $OMP_NUM_THREADS. The completed job''s log file is valid, so rerunning snakemake resumes where it stopped.'
---

Run: snakemake -s Snakefile --configfile config.yaml -j $OMP_NUM_THREADS pipeline_all in ~/v2.0/image_sims/grids on candide (node c02). Crash in snakemake 9.19.0 job_scheduler._schedule_reevalutation -> threading.Timer().start() right after tile 288.279 job 16 completed.
