# cfis_image_sims — module configs for ShapePipe on simulated tiles

This directory is the single `.ini` tree for image-simulation runs. It is
selected by `run_job_sp_canfar_v2.0.bash -t image_sims` (which sets
`config_dir = example/cfis_image_sims`, `retrieve = symlink`, and forces
`tile_det = sx`) and consumed stage-by-stage by `job_sp_canfar_v2.0.bash`.

The pipeline runs as a bit-coded chain of jobs. `run_job_sp_canfar_v2.0.bash`
loops over the bits set in `-j`, initialises the tile/exposure work directories,
runs each completeness check, and delegates the actual `shapepipe_run` calls to
`job_sp_canfar_v2.0.bash -j <bit>`, which is where the `.ini` file for each bit
is selected. A few bits are special-cased for sims *before* that delegation —
those cases are documented in the last column below and, at more length, under
[Sim special-casing](#sim-special-casing).

## Job-bit dispatch

Every row is derived from the two bash scripts. "Module(s)" is the ShapePipe
runner(s) the selected `.ini` names; "`.ini` selected" is what
`job_sp_canfar_v2.0.bash` picks for that bit under sim settings
(`retrieve=symlink`, `psf=psfex`, `tile_det=sx`, `star_cat_for_mask=onthefly`).

| Bit | Stage | Module(s) | `.ini` selected (sim settings) | Sim special-casing |
|----:|-------|-----------|--------------------------------|--------------------|
| 1 | retrieve tile image + weight | `get_images_runner` | `config_tile_Git_symlink.ini` | symlink retrieval (see below); completeness check expects 2 files (image + weight) vs. 4 for data |
| 2 | uncompress tile weight | *(none — faked)* | *(none)* | **Faked weight-uncompress.** run_job skips the `config_tile_Uz.ini` run entirely; sim weights are already uncompressed, so it fabricates a `run_sp_tile_Uz*/uncompress_fits_runner/output` dir and symlinks `input_tiles/CFIS_simu_weight-<ID-dashed>.fits` into it, satisfying downstream `last:uncompress_fits_runner` lookups |
| 4 | find exposures | `find_exposures_runner` | `config_tile_Fe.ini` | — |
| 8 | retrieve exposure images | `get_images_runner` | `config_exp_Gie_symlink.ini` | symlink retrieval; completeness check expects 3 files vs. 6 for data |
| 16 | split exposures, merge WCS headers | `split_exp_runner` | `config_exp_Sp.ini` | — |
| 32 | mask exposures | `mask_runner` | `config_exp_Ma_onthefly.ini` | — |
| 64 | exposure PSF model | *(none — placeholder)* | *(none)* | **Placeholder.** For data this runs full exposure PSF modelling. For sims run_job writes a placeholder log and does nothing here; the sim PSF (`fake_psf_runner`) actually runs inside bit 512 |
| 128 | merge exposure WCS headers → tile sqlite log | `merge_headers_runner` | `config_tile_Mh_exp.ini` | — |
| 256 | object detection on tiles | `sextractor_runner` | `config_tile_Sx_nomask.ini` | `tile_det` is forced to `sx`, so the SExtractor-no-mask branch is always taken (the `uc` external-catalogue branch is never reached for sims) |
| 512 | fake PSF + postage stamps | `fake_psf_runner`, then `vignetmaker_runner` ×2 | `config_exp_psfex.ini` (fake PSF), then `config_tile_PiViVi_canfar_sx.ini` (vignets) | **Two sub-runs.** run_job first calls the job script with `-j 64` → `config_exp_psfex.ini`, which despite its name runs `fake_psf_runner` (needs the sexcat from bit 256; run dir `run_sp_tile_fpsf`), then `-j 512` → `config_tile_PiViVi_canfar_sx.ini` for the two `vignetmaker_runner` runs. Data instead runs `psfex_interp_runner` + vignets here |
| 1024 | multi-epoch shape measurement | `ngmix_runner` | `config_tile_Ng_batch_psfex_sx.ini` | — |
| 2048 | create final catalogue | `make_cat_runner` | `config_tile_Mc_psfex.ini` | — |

Notes on cross-cutting conventions used above:

- **Dashed tile-id convention.** Data tile IDs use dot format (`233.293`); sim
  input files use dash format (`233-293`). run_job writes `tile_numbers.txt` in
  dash form for `image_sims` (`${ID//./-}`), and the faked bit-2 weight symlink
  and sim input patterns all key off the dashed ID.
- **Symlink retrieval.** `-t image_sims` sets `retrieve=symlink`, so bits 1 and 8
  select the `*_symlink.ini` get-images configs (`RETRIEVE = symlink`), which
  link sim images out of `$SP_DIR/input_tiles` / `$SP_DIR/input_exp` rather than
  downloading from VOSpace.
- **Forced `tile_det=sx`.** The `image_sims` type branch forces `tile_det=sx`, so
  bit 256 always runs SExtractor detection and bit 512 selects the `_sx`-suffixed
  vignet/ngmix/PiViVi configs.

## Sim special-casing

Three bits diverge from the data chain, all handled in
`run_job_sp_canfar_v2.0.bash` before the per-bit delegation:

- **Bit 2 — faked weight-uncompress.** Sim weights ship uncompressed, so instead
  of running `uncompress_fits_runner`, run_job fabricates the expected Uz output
  directory and symlinks `CFIS_simu_weight-<ID-dashed>.fits` in, so downstream
  `last:uncompress_fits_runner` references resolve.
- **Bit 64 — placeholder.** run_job writes a placeholder completeness log and
  runs nothing; the sim PSF is produced by `fake_psf_runner` inside bit 512.
- **Bit 512 — fake PSF + vignets.** run_job runs `fake_psf_runner`
  (`config_exp_psfex.ini`) followed by the two `vignetmaker_runner` runs
  (`config_tile_PiViVi_canfar_sx.ini`). `fake_psf_runner` consumes the sexcat
  from bit 256, so bit 256 must precede it.

## Parallelism model and the `PSF_DICT` requirement

The orchestrator (the sp_validation Snakemake workflow) fans out one job per
tile; within a single job ShapePipe uses its own multiprocessing, controlled by
`-N/--N_SMP` (threaded through to each module's `SMP_BATCH_SIZE`). Every sim
config sets `MODE = SMP` — the sims chain launches no MPI.

`fake_psf_runner` reads the sims' pickled PSF dictionary from the `$PSF_DICT`
environment variable (`config_exp_psfex.ini` sets `PSF_DICT_PATH = $PSF_DICT`,
expanded from the environment). The orchestrating workflow exports it (from its
`psf_dict` config key); on candide it points at the SKiLLS PSF dictionary. Bit
512 will fail if `$PSF_DICT` is unset.
