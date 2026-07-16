"""Tile post-chain — per tile: gather exposures, then detect/PSF/shape/catalogue.

A per-tile exposure "forest" (symlink view of exactly this tile's exposures'
products) gives the tile modules a deterministic $SP_EXP to glob. But the forest
is only a runtime convenience, NEVER the DAG edge: every tile rule that reads
through $SP_EXP declares the underlying exposure directory() outputs as input via
an index dict lookup, so the last real reader is the last DAG consumer (findings
5/19). All rules are group-compatible (shell only, no mid-chain localrules).
"""

TILE_OUT = str(RUN_DIR / "tiles/{tile}/output")
FOREST   = str(RUN_DIR / "tiles/{tile}/exp_forest")


def _exp_out(wc, prefix):
    return [str(exp_dir(e) / "output" / prefix) for e in TILE_EXP.get(wc.tile, [])]

def tile_exp_split(wc):    return _exp_out(wc, "run_sp_exp_Sp")
def tile_exp_mask(wc):     return _exp_out(wc, "run_sp_exp_Ma")
def tile_exp_psf(wc):      return _exp_out(wc, "run_sp_exp_SxSePsfPi")
def tile_exp_products(wc):
    """All of this tile's exposures' product dirs — the forest's DAG inputs."""
    return tile_exp_split(wc) + tile_exp_mask(wc) + tile_exp_psf(wc)


# Build the per-tile symlink forest (exp_forest/<base>/output/run_sp_exp_* -> the
# real exposure products). Declaring the exposure dirs as input makes this the
# tile's wait-on-its-exposures edge; the forest itself is just the SP_EXP view.
rule tile_exp_forest:
    input:
        tile_exp_products
    output:
        directory(FOREST)
    params:
        # --forest {output} lives in the shell string: Snakemake formats shell
        # ONCE — an {output} placeholder inside params.cmd survives literally
        # (same trap as {threads}; all four forest jobs then race one './{output}').
        cmd=lambda wc: (f"python {SCRIPTS}/build_forest.py --tile {wc.tile} "
                        f"--run-dir {RUN_DIR} --index {INDEX_DB}")
    shell:
        # no --threads: build_forest.py is single-threaded symlinking
        "{params.cmd} --forest {output}"

# Merge single-exposure WCS headers into the tile-level sqlite log
# (log_exp_headers-<IDra>-<IDdec>.sqlite — the tile suffix downstream Sx/PiViVi/Ng
# consume via NUMBERING_SCHEME -000-000). Reads headers-*.npy through the forest,
# so the exposure split dirs are explicit inputs.
rule tile_merge_headers:
    input:
        forest=rules.tile_exp_forest.output,
        split=tile_exp_split,
    output:
        directory(f"{TILE_OUT}/run_sp_tile_Mh_exp")
    params:
        cmd=lambda wc: sp_rule("tile_merge_headers", "config_tile_Mh_exp.ini",
                               "tile", wc.tile, isolate=True,
                               exp_forest=FOREST.format(tile=wc.tile))
    shell:
        "{params.cmd} --threads {threads}"

# Mask tiles (star cats: pre-run store, read via the wrapper's dir symlink).
rule tile_mask:
    input:
        git=rules.tile_get_images.output,
    output:
        directory(f"{TILE_OUT}/run_sp_tile_Ma")
    params:
        cmd=lambda wc: sp_rule("tile_mask", "config_tile_Ma_onthefly.ini",
                               "tile", wc.tile, isolate=True)
    shell:
        "{params.cmd} --threads {threads}"

# SExtractor object detection on the tile.
rule tile_detect:
    input:
        mask=rules.tile_mask.output,
        mh=rules.tile_merge_headers.output,
    output:
        directory(f"{TILE_OUT}/run_sp_tile_Sx")
    params:
        cmd=lambda wc: sp_rule("tile_detect", "config_tile_Sx.ini",
                               "tile", wc.tile, isolate=True)
    shell:
        "{params.cmd} --threads {threads}"

# PSFEx interpolation to galaxies + vignet postage stamps (reads exposure PSF
# models + split through the forest — those dirs are explicit inputs).
rule tile_vignets:
    input:
        sx=rules.tile_detect.output,
        forest=rules.tile_exp_forest.output,
        split=tile_exp_split,
        psf=tile_exp_psf,
    output:
        directory(f"{TILE_OUT}/run_sp_tile_PiViVi")
    params:
        cmd=lambda wc: sp_rule("tile_vignets", "config_tile_PiViVi_canfar_sx.ini",
                               "tile", wc.tile, isolate=True,
                               exp_forest=FOREST.format(tile=wc.tile))
    shell:
        "{params.cmd} --threads {threads}"

# ngmix shape measurement — static N chunks. Each chunk computes its own closed
# object-ID range in-job from the tile's own sexcat (a declared input); chunks
# never rewrite the shared log (--no-log-sync) and read literal INPUT_DIRs.
rule tile_ngmix:
    input:
        vignets=rules.tile_vignets.output,
        sx=rules.tile_detect.output,
    output:
        directory(f"{TILE_OUT}/run_sp_tile_ngmix_Ng{{chunk}}u")
    params:
        cmd=lambda wc: sp_rule(
            "tile_ngmix", "config_tile_Ng_template_batch.ini", "tile", wc.tile,
            extra=f"--ngmix-chunk {wc.chunk} --ngmix-nchunks {NGMIX_CHUNKS}"),
    shell:
        "{params.cmd} --threads {threads}"

def ngmix_chunks(wc):
    return [f"{RUN_DIR}/tiles/{wc.tile}/output/run_sp_tile_ngmix_Ng{k}u"
            for k in range(1, NGMIX_CHUNKS + 1)]

# Merge the N chunk catalogues (DAG-serialized: the single per-tile log-sync
# point before make_cat). N_SPLIT_MAX == chunk count.
rule tile_merge_cats:
    input:
        ngmix_chunks
    output:
        directory(f"{TILE_OUT}/run_sp_Ms")
    params:
        cmd=lambda wc: sp_rule("tile_merge_cats", "config_merge_sep_cats_template.ini",
                               "tile", wc.tile, isolate=True,
                               extra=f"--n-split-max {NGMIX_CHUNKS}"),
    shell:
        "{params.cmd} --threads {threads}"

# Build the final catalogue: the run's science product. protected() so a
# code/params-drift rerun cannot silently discard a finished catalogue — the
# recompute must go through `sp rerun` (--forcerun).
rule tile_make_cat:
    input:
        rules.tile_merge_cats.output
    output:
        protected(str(RUN_DIR / "tiles/{tile}/final_cat-{tile}.fits"))
    params:
        cmd=lambda wc: sp_rule("tile_make_cat", "config_make_cat_psfex_nosm.ini",
                               "tile", wc.tile, isolate=True,
                               extra=f"--final-cat {RUN_DIR}/tiles/{wc.tile}/final_cat-{wc.tile}.fits"),
    shell:
        "{params.cmd} --threads {threads}"
