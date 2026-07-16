"""Phase A — prepare: tile staging, find_exposures, per-unit star cats.

Split into two static invocations with build_index.py (a plain script) between:
  * prepare_tiles      — tile_get_images -> tile_uncompress -> tile_find_exposures
  * [build_index.py]   — plain script, threshold-gated (not a DAG node)
  * prepare_exposures  — exp_get_images (exposure.smk) + per-unit star cats

Nibi compute nodes have internet, so downloads run in-DAG (no login-node tier).
Rules stay group-compatible: no mid-chain localrules, no run:/pipe outputs.
"""

# get_images reads a one-line tile_numbers.txt the wrapper drops in $SP_RUN.
# NUMBER_LIST is NEVER injected here (download stage; nothing on disk to
# validate against — #746 would hard-fail).
rule tile_get_images:
    output:
        directory(str(RUN_DIR / "tiles/{tile}/output/run_sp_tile_Git"))
    params:
        cmd=lambda wc: sp_rule("tile_get_images", "config_tile_Git_vos.ini",
                               "tile", wc.tile, isolate=False)
    shell:
        "{params.cmd} --threads {threads}"

rule tile_uncompress:
    input:
        rules.tile_get_images.output
    output:
        directory(str(RUN_DIR / "tiles/{tile}/output/run_sp_tile_Uz"))
    params:
        cmd=lambda wc: sp_rule("tile_uncompress", "config_tile_Uz.ini",
                               "tile", wc.tile, isolate=True)
    shell:
        "{params.cmd} --threads {threads}"

# find_exposures reads the tile FITS HISTORY header -> exp_numbers-<ID>.txt, the
# data-derived tile->exposure edge build_index.py aggregates (as a plain script).
rule tile_find_exposures:
    input:
        rules.tile_uncompress.output
    output:
        directory(str(RUN_DIR / "tiles/{tile}/output/run_sp_tile_Fe"))
    params:
        cmd=lambda wc: sp_rule("tile_find_exposures", "config_tile_Fe.ini",
                               "tile", wc.tile, isolate=True)
    shell:
        "{params.cmd} --threads {threads}"

# Star catalogues for masking are pre-generated offline (create_star_cat.py on a
# networked login node) and consumed as DIRECTORIES: the mask configs read
# INPUT_DIR $SP_RUN/star_cat_{exp,tiles} with per-CCD numbering (exp cats are
# star_cat-<exp>-<ccd>.fits, 40/exposure; tile cats star_cat-<IDra>-<IDdec>.fits)
# — the v2.0 mechanism. The wrapper materialises the two dir symlinks in every
# unit work dir (sp_rule.materialise_unit); there is NO per-unit star-cat DAG
# node: the store is pre-run input like the image store, and a missing cat fails
# the mask stage's count-floor loudly. (Earlier per-unit file rules linked names
# that don't exist — one file per exposure vs the store's 40 — and nothing read
# them: the configs only see the dir symlinks.) When star-cat GENERATION moves
# in-workflow (P2), it becomes real per-unit rules producing into the store.
