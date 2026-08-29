"""Invocation 1 — PREPARE: one static chain per tile.

    tile_get_images -> tile_uncompress -> tile_find_exposures

Known from the tile list alone, cheap, wide, idempotent. find_exposures parses
the tile FITS HISTORY header into ``exp_numbers-<IDra>-<IDdec>.txt`` — the
data-derived tile->exposure edge that invocation 2's parse aggregates into the
index. Nibi compute nodes have internet, so downloads run in-DAG (no login-node
tier).

Rules stay group-compatible: shell only, no mid-chain localrules, no pipe outputs.

Star catalogues for masking are NOT a prepare-phase concern and not pre-run
input: the compute DAG fetches the campaign footprint's stars once
(``star_catalogue``) and cuts them per exposure (``exp_star_cat``), both in
exposure.smk, into a run-independent store. The tile side has no star-cat node
because it has no mask rule yet — see tile.smk.
"""

# NUMBER_LIST is never set for get_images (download stage; nothing on disk to
# validate against — #746 would hard-fail the unit). The committed
# config_tile_Git.ini simply has no NUMBER_LIST, so this is now a property of the
# config, not of an injection step.
rule tile_get_images:
    output:
        manifest = f"{TILE_DIR}/manifests/tile_get_images.json"
    log:
        f"{TILE_DIR}/logs/tile_get_images.json"
    params:
        pre = lambda wc: unit_pre("tile_get_images", "tile", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 1
    retries: 2
    resources:
        mem_mb = lambda wc, attempt: 4000 * attempt,
        runtime = 60
    shell:
        sp_shell("tile_get_images", "config_tile_Git.ini")

rule tile_uncompress:
    input:
        rules.tile_get_images.output.manifest
    output:
        manifest = f"{TILE_DIR}/manifests/tile_uncompress.json"
    log:
        f"{TILE_DIR}/logs/tile_uncompress.json"
    params:
        pre = lambda wc: unit_pre("tile_uncompress", "tile", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 4
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 60
    shell:
        sp_shell("tile_uncompress", "config_tile_Uz.ini")

rule tile_find_exposures:
    input:
        rules.tile_uncompress.output.manifest
    output:
        manifest = f"{TILE_DIR}/manifests/tile_find_exposures.json"
    log:
        f"{TILE_DIR}/logs/tile_find_exposures.json"
    params:
        pre = lambda wc: unit_pre("tile_find_exposures", "tile", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 1
    resources:
        mem_mb = lambda wc, attempt: 2000 * attempt,
        runtime = 30
    shell:
        sp_shell("tile_find_exposures", "config_tile_Fe.ini")
