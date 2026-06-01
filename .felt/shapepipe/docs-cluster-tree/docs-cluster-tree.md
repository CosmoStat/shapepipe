---
name: Machine-specific cluster docs tree + freshness pass
tags:
    - shapepipe
    - docs
created-at: 2026-05-31T22:15:20.539028984+02:00
outcome: 'PR #739 (off develop): introduced a single clusters.md (general pattern + candide/CANFAR/ccin2p3 sections) under a ''Running on a cluster'' toctree caption — chosen over a standalone general page (too thin: shared content is ~3 paragraphs, submission models diverge) and over per-machine pages (canfar walkthrough would unbalance). Deleted obsolete canfar.md/pipeline_v2.0.md/work_flow_v2.0.md (planning rescued to [[shapepipe/v2-run-plan]]); kept pipeline_canfar.md as the linked deep CANFAR walkthrough. Rewrote dependencies.md against pyproject (SSOT framing, ngmix fork, dropped CDSclient); stripped removed rho-stats + dead prepare_tiles_for_final from post_processing/random_cat; cosmetics. Audit found install/container/testing/API pages already fresh. Local sphinx build green. Sequencing: candide section assumes #737''s SLURM scripts; publishability needs [[shapepipe/docs-versioning]] (#738).'
---
