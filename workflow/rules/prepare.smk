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
        "{params.cmd}"

rule tile_uncompress:
    input:
        rules.tile_get_images.output
    output:
        directory(str(RUN_DIR / "tiles/{tile}/output/run_sp_tile_Uz"))
    params:
        cmd=lambda wc: sp_rule("tile_uncompress", "config_tile_Uz.ini",
                               "tile", wc.tile, isolate=True)
    shell:
        "{params.cmd}"

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
        "{params.cmd}"

# Star catalogues for masking are pre-generated offline (create_star_cat.py on a
# networked login node — the one genuinely networked prepare step). They are
# re-keyed PER UNIT here: each tile/exposure's cat is its own DAG node, so
# appending tile 51 adds nodes rather than mutating a run-level output that would
# dirty every finished mask (the growth-invariant fix, finding 11).
rule exp_star_cat:
    output:
        str(RUN_DIR / "exp/{exp}/star_cat-{exp}.fits")
    params:
        src=lambda wc: str(STAR_CATS / "exp" / f"star_cat-{wc.exp}.fits")
    shell:
        "ln -sf {params.src} {output}"

rule tile_star_cat:
    output:
        str(RUN_DIR / "tiles/{tile}/star_cat-{tile}.fits")
    params:
        src=lambda wc: str(STAR_CATS / "tiles" / f"star_cat-{wc.tile}.fits")
    shell:
        "ln -sf {params.src} {output}"
