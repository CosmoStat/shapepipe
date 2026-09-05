"""Invocation 1 — PREPARE: one static chain per tile.

    tile_get_images -> tile_uncompress -> tile_find_exposures

Known from the tile list alone, cheap, wide, idempotent. find_exposures parses
the tile FITS HISTORY header into ``exp_numbers-<IDra>-<IDdec>.txt`` — the
data-derived tile->exposure edge that invocation 2's parse aggregates into the
index. Nibi compute nodes have internet, so downloads run in-DAG (no login-node
tier).

All three rules carry ``group: "tile_prep"``, so one tile's whole chain is ONE
sbatch instead of three (medians 0:41 / ~0:40 / 0:15 — all far under the
15-minute runtime limit Alliance policy asks us to bundle away, and at DR6 scale three
submissions per tile is a scheduler load out of all proportion to the work).
Group membership is per connected DAG component and distinct tiles share no
edge, so this is exactly one group job per tile, never a cross-tile bundle.
Rules stay group-compatible: shell only, no mid-chain localrules, no pipe outputs.

Group resource composition (snakemake 9.23, ``GroupResources.basic_layered`` in
snakemake/resources.py): jobs are laid out per toposort level; within a level
non-additive resources (mem_mb, cpus) SUM — split into layers when a global
constraint is exceeded, the group's width being the widest layer — while the
additive resource ``runtime`` is maxed within a layer and SUMMED across layers.
This chain is strictly linear, one job per level, so the group asks for
max(mem_mb) = 8000*attempt, max(threads) = 4 and sum(runtime) = 150 min.
Attempt scaling survives grouping: ``GroupJob.attempt``'s setter clears the
cached group resources and re-sets ``attempt`` on every member (jobs.py), and
``GroupJob.restart_times`` is the max over members — so tile_get_images'
``retries: 2`` still governs. A retry re-runs the whole group, which is safe
because every rule ``rm -rf``s its own run dir at start.

Star catalogues for masking are NOT a prepare-phase concern and not pre-run
input: the compute DAG fetches the campaign footprint's stars once
(``star_catalogue``) and cuts them per exposure (``exp_star_cat``), both in
exposure.smk, into a run-independent store. The tile side has no star-cat node
because it has no mask rule yet — see tile.smk.
"""

# No NUMBER_LIST for get_images — a download stage has nothing on disk to
# validate against (see exposure.smk's docstring for the convention).
rule tile_get_images:
    group: "tile_prep"
    output:
        manifest = f"{TILE_DIR}/manifests/tile_get_images.json"
    log:
        f"{TILE_DIR}/logs/tile_get_images.json"
    params:
        pre = lambda wc: unit_pre("tile_get_images", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = lambda wc, attempt: 4000 * attempt,
        runtime = 60
    shell:
        sp_shell("tile_get_images", "config_tile_Git.ini")

rule tile_uncompress:
    group: "tile_prep"
    input:
        rules.tile_get_images.output.manifest
    output:
        manifest = f"{TILE_DIR}/manifests/tile_uncompress.json"
    log:
        f"{TILE_DIR}/logs/tile_uncompress.json"
    params:
        pre = lambda wc: unit_pre("tile_uncompress", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 4
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 60
    shell:
        sp_shell("tile_uncompress", "config_tile_Uz.ini")

rule tile_find_exposures:
    group: "tile_prep"
    input:
        rules.tile_uncompress.output.manifest
    output:
        manifest = f"{TILE_DIR}/manifests/tile_find_exposures.json"
    log:
        f"{TILE_DIR}/logs/tile_find_exposures.json"
    params:
        pre = lambda wc: unit_pre("tile_find_exposures", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: 2000 * attempt,
        runtime = 30
    shell:
        sp_shell("tile_find_exposures", "config_tile_Fe.ini")
