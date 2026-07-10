---
title: Blank input tiles grid_4
depends-on:
    - snakemake-crash-runtimeerror-20c2564f
created-at: 2026-07-08T17:09:48.720667191+02:00
outcome: '12 of 37 input tiles in /n09data/hervas/skills_out/*_grid_4/images/SP_tiles are all-zero (image AND weight), identical set in all 4 shear sims: 276.280 276.281 277.282 280.280 280.282 281.281 282.281 284.279 285.279 285.280 286.279 286.280. Upstream sim-generation problem, needs regeneration by Fabian (hervas). Workaround: new tile_IDs_exclude config key (Snakefile + info.py), active in run config.yaml. Also added: log directive err_job_{tile}_{bit}.txt with tee so errors survive failure; config.yaml no longer tracked input of rule init (edits cascaded into full re-run); standard command now uses --rerun-triggers=mtime and fixed -j 4.'
---

Snakemake pipeline_all failed at job 512 (vignetmaker) for tile 281.281 in two sims: SExtractor (job 256) succeeds but finds 0 objects, writes header-only sexcat, vignetmaker raises ValueError empty/corrupt catalogue. Snakemake deletes outputs of failed jobs, so the shapepipe log vanished; real error was in tiles/<xxx>/<tile>/output/run_sp_tile_ViVi_*/vignetmaker_runner_run_1/logs/process-*.log
