"""Exposure chain — per exposure, keyed by exp_id (dedup is structural).

Gie -> Sp -> Ma -> psf, each in the exposure's own work dir, chained by
ShapePipe's last:/run_sp_exp_* INPUT_DIR resolver. Completeness is the count
floor (completeness.py).

NO temp() on these shared exposure directories: exposures overlap tiles by
construction (~7-10 tiles per exposure), so temp() would delete an exposure the
moment its forests are built and cascade destructive reruns across spatial
neighbours whenever a tile is appended (findings 4/5/6/19). The store is
persistent; space is reclaimed by the explicit index-driven `sp clean-exposures`
verb (scripts/clean_exposures.py), which deletes an exposure's intermediates
only once every consuming tile in the index has its final_cat present.

NUMBER_LIST is injected only for exp_split (numbering scheme == exposure id);
NEVER for get_images (download) / exp_mask / exp_psf (per-CCD schemes where the
#746 startup validation would turn tolerated per-CCD attrition into a
whole-exposure hard failure).
"""

EXP_OUT = str(RUN_DIR / "exp/{exp}/output")

# Symlink/download exposure image/weight/flag; reads the fabricated pseudo-Fe
# exp_numbers-000-000.txt the wrapper drops (via last:find_exposures_runner).
rule exp_get_images:
    output:
        directory(f"{EXP_OUT}/run_sp_exp_Gie")
    params:
        cmd=lambda wc: sp_rule("exp_get_images", "config_exp_Gie_vos.ini",
                               "exp", wc.exp, isolate=False,
                               extra=f"--exp-name {EXP[wc.exp]}")
    shell:
        "{params.cmd} --threads {threads}"

# Split multi-HDU exposure into single-CCD files (+ headers-*.npy). Bulk output;
# persistent (see module docstring on temp()).
rule exp_split:
    input:
        rules.exp_get_images.output
    output:
        directory(f"{EXP_OUT}/run_sp_exp_Sp")
    params:
        cmd=lambda wc: sp_rule("exp_split", "config_exp_Sp.ini", "exp", wc.exp)
    shell:
        "{params.cmd} --threads {threads}"

# Mask per-CCD; this exposure's star cat is a declared input (per-unit re-key).
rule exp_mask:
    input:
        rules.exp_split.output,
        star=str(RUN_DIR / "exp/{exp}/star_cat-{exp}.fits"),
    output:
        directory(f"{EXP_OUT}/run_sp_exp_Ma")
    params:
        cmd=lambda wc: sp_rule("exp_mask", "config_exp_Ma_onthefly.ini",
                               "exp", wc.exp, isolate=False)
    shell:
        "{params.cmd} --threads {threads}"

# SExtractor -> setools star selection -> PSFEx model -> psfex_interp, per CCD.
# setools may reject a sparse CCD (~0.2% attrition) — tolerated by the floor.
rule exp_psf:
    input:
        rules.exp_mask.output
    output:
        directory(f"{EXP_OUT}/run_sp_exp_SxSePsfPi")
    params:
        cmd=lambda wc: sp_rule("exp_psf", "config_exp_psfex.ini", "exp", wc.exp,
                               isolate=False)
    shell:
        "{params.cmd} --threads {threads}"
