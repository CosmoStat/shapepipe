"""Tile post-chain — per tile: gather exposures, then detect / PSF / shape / catalogue.

    tile_exp_forest
    tile_merge_headers -> tile_detect -> tile_vignets -> tile_ngmix x N
                                                      -> tile_merge_cats -> tile_make_cat

The DAG edge to the exposures is always the exposures' MANIFESTS, looked up
through the index (TILE_EXP). The per-tile exposure "forest" — a symlink view of
exactly this tile's exposures' products — exists only so the ShapePipe configs
have one deterministic ``$SP_EXP`` to glob; it is NEVER the edge. Its 2-char
shard level is not cosmetic: ``exp_utils.get_exp_output_files`` hardwires
``<SP_EXP>/<prefix>/<base>/output/run_sp_*`` into its glob, so a flat forest
makes every tile gather stage fail "No split_exp_runner output found".

All rules are group-compatible (shell only, no mid-chain localrules).

Note there is no `tile_mask` rule: the committed config chain is the
"sx_nomask" tile_detect variant (config_tile_Sx.ini reads Git + Uz + Mh, no mask
run), and no tile-mask config was committed in the S2 sweep. Adding the masked
variant is a config + one rule, at the config selector the PRD describes.
"""

def tile_exp(wc):
    return TILE_EXP.get(wc.tile, [])

def exp_manifests(wc, stage):
    return [exp_manifest(e, stage) for e in tile_exp(wc)]

def tile_exp_split(wc):  return exp_manifests(wc, "exp_split")
def tile_exp_mask(wc):   return exp_manifests(wc, "exp_mask")
def tile_exp_psf(wc):    return exp_manifests(wc, "exp_psf")
def tile_exp_all(wc):    return tile_exp_split(wc) + tile_exp_mask(wc) + tile_exp_psf(wc)


# Build the per-tile symlink forest. Declaring the exposure manifests as input
# makes this wait on its exposures; the forest itself is only the $SP_EXP view.
# Its output stays a directory() (it has no ShapePipe run dir and no manifest —
# it is not a shapepipe_run at all).
rule tile_exp_forest:
    input:
        tile_exp_all
    output:
        forest = directory(f"{TILE_DIR}/exp_forest")
    params:
        cmd = lambda wc: (f"python {SCRIPTS}/build_forest.py --tile {wc.tile} "
                          f"--run-dir {RUN_DIR} --index {INDEX_DB}"),
        # build_forest.py's own content hash rides here and nowhere else.
        script_hash = FOREST_HASH
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 20
    shell:
        # --forest {output} lives in the shell string: snakemake formats shell
        # ONCE, so an {output} placeholder inside params.cmd would survive
        # literally and every forest job would race one './{output}'.
        "{params.cmd} --forest {output.forest}"

# Merge single-exposure WCS headers into the tile-level sqlite
# (log_exp_headers-<IDra>-<IDdec>.sqlite, which Sx / PiViVi / ngmix consume).
# Reads headers-*.npy through the forest -> the split manifests are the edge.
rule tile_merge_headers:
    input:
        forest = rules.tile_exp_forest.output.forest,
        split  = tile_exp_split,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_merge_headers.json"
    params:
        pre = lambda wc: unit_pre("tile_merge_headers", "tile", wc.tile,
                                  forest=forest_dir(wc.tile)),
        script_hash = SCRIPT_HASH
    threads: 4
    resources:
        mem_mb = lambda wc, attempt: 8000 * attempt,
        runtime = 120
    shell:
        sp_shell("tile_merge_headers", "config_tile_Mh_exp.ini")

# SExtractor object detection on the tile.
rule tile_detect:
    input:
        uz = f"{TILE_DIR}/manifests/tile_uncompress.json",
        mh = rules.tile_merge_headers.output.manifest,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_detect.json"
    params:
        pre = lambda wc: unit_pre("tile_detect", "tile", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 180
    shell:
        sp_shell("tile_detect", "config_tile_Sx.ini")

# PSFEx interpolation to galaxies + vignet postage stamps: the last stage that
# reads exposure products, and the bulk intra-tile intermediate.
#
# The vignette store is declared as a SECOND, temp(directory()) output alongside
# the manifest. This is the ONE scoped exception to "no directory() outputs"
# (D5): the store is ~tens of GB per tile and must be reclaimed when its last
# intra-tile reader finishes, but it must not become DAG currency — so the
# manifest stays the edge, and the directory rides along purely so native temp()
# fires at the right moment. Its readers (tile_ngmix, tile_make_cat) declare
# BOTH. `--notemp` keeps it for debugging.
rule tile_vignets:
    input:
        sx     = rules.tile_detect.output.manifest,
        forest = rules.tile_exp_forest.output.forest,
        split  = tile_exp_split,
        psf    = tile_exp_psf,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_vignets.json",
        store    = temp(directory(f"{TILE_DIR}/output/run_sp_tile_PiViVi")),
    params:
        pre = lambda wc: unit_pre("tile_vignets", "tile", wc.tile,
                                  forest=forest_dir(wc.tile)),
        script_hash = SCRIPT_HASH
    threads: 16
    resources:
        mem_mb = lambda wc, attempt: 32000 * attempt,
        runtime = 240
    shell:
        sp_shell("tile_vignets", "config_tile_PiViVi.ini")

# ngmix shape measurement — N chunks per tile (D4). Each chunk computes its own
# CLOSED object-ID range at EXECUTION time from this tile's own sexcat: a params
# function cannot, because params evaluate before the sexcat exists. Closed, not
# open-ended: `ID_OBJ_MAX = -1` on the last chunk was the 13-hour straggler's
# root cause (ngmix treats id_obj_max <= 0 as unbounded).
#
# Chunks write nothing shared: each has its own run_sp_tile_ngmix_Ng<k>u, and
# merge_sep_cats — DAG-serialised after all chunks — is the gather.
rule tile_ngmix:
    input:
        vignets  = rules.tile_vignets.output.manifest,
        store    = rules.tile_vignets.output.store,
        sx       = rules.tile_detect.output.manifest,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.json",
        # temp(directory()) for the same reason as the vignette store above.
        chunkdir = temp(directory(f"{TILE_DIR}/output/run_sp_tile_ngmix_Ng{{chunk}}u")),
    params:
        pre = lambda wc: unit_pre(
            "tile_ngmix", "tile", wc.tile,
            env={"SP_NGMIX_CHUNK": wc.chunk, "NGMIX_N_CHUNKS": NGMIX_CHUNKS},
            pre_run=[f'eval "$(python {SCRIPTS}/ngmix_range.py --run-dir '
                     f'"$SP_RUN" --chunk {wc.chunk} --n-chunks {NGMIX_CHUNKS})"']),
        script_hash = SCRIPT_HASH
    threads: 4
    retries: 2
    benchmark:
        f"{TILE_DIR}/manifests/tile_ngmix_{{chunk}}.benchmark.tsv"
    resources:
        mem_mb = lambda wc, attempt: 14000 * attempt,
        runtime = 720
    shell:
        sp_shell("tile_ngmix", "config_tile_Ng_template.ini")


def ngmix_manifests(wc):
    return [f"{tile_dir(wc.tile)}/manifests/tile_ngmix_{k}.json"
            for k in range(1, NGMIX_CHUNKS + 1)]

def ngmix_chunkdirs(wc):
    return [f"{tile_dir(wc.tile)}/output/run_sp_tile_ngmix_Ng{k}u"
            for k in range(1, NGMIX_CHUNKS + 1)]

# The gather: merge the N chunk catalogues. N_SPLIT_MAX comes from the workflow's
# own chunk count via $NGMIX_N_CHUNKS (env-expanded by the module).
rule tile_merge_cats:
    input:
        manifests = ngmix_manifests,
        chunkdirs = ngmix_chunkdirs,
    output:
        manifest = f"{TILE_DIR}/manifests/tile_merge_cats.json"
    params:
        pre = lambda wc: unit_pre("tile_merge_cats", "tile", wc.tile,
                                  env={"NGMIX_N_CHUNKS": NGMIX_CHUNKS}),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 120
    shell:
        sp_shell("tile_merge_cats", "config_merge_sep_cats.ini")

# The run's science product. make_cat also reads the vignette store's
# psfex_interp output, so it — not ngmix — is the store's last reader.
#
# No protected(): the full default rerun-triggers govern, and protected() only
# ever forced people through a `--forcerun` detour.
rule tile_make_cat:
    input:
        ms    = rules.tile_merge_cats.output.manifest,
        store = rules.tile_vignets.output.store,
    output:
        manifest  = f"{TILE_DIR}/manifests/tile_make_cat.json",
        final_cat = f"{TILE_DIR}/final_cat-{{tile}}.fits",
    params:
        pre = lambda wc: unit_pre("tile_make_cat", "tile", wc.tile),
        script_hash = SCRIPT_HASH
    threads: 8
    resources:
        mem_mb = lambda wc, attempt: 16000 * attempt,
        runtime = 120
    shell:
        # Publish the catalogue next to the manifest: a real file, so it is a
        # real declared output (and it persists — never temp()).
        "{params.pre}\n"
        "rc=0\n"
        'shapepipe_run -c "$SP_CONFIG/config_tile_Mc.ini" -b {threads} || rc=$?\n'
        f"python {SCRIPTS}/completeness.py check tile_make_cat {{output.manifest}} || rc=1\n"
        "if [ $rc -eq 0 ]; then\n"
        '  cp -f "$(ls -1 "$SP_RUN"/output/run_sp_Mc/make_cat_runner/output/final_cat*.fits'
        ' | head -1)" {output.final_cat}\n'
        "fi\n"
        "exit $rc\n"
