# Image Simulations Calibration Pipeline

## Overview

The image simulations calibration pipeline derives multiplicative (m) and additive (c) shear bias corrections from synthetic image simulations. It uses ShapePipe to measure galaxy shapes in five pre-sheared grids, then computes bias by comparing measurements across shear directions.

**Key outputs:**
- `m_bias_results.yaml` / `m_bias_results.txt` — final m₁, m₂, c₁, c₂ bias estimates with bootstrap errors
- `mbias_cumulative.yaml` — convergence history (m/c as function of n_tiles)
- `mbias_convergence.png` — m and c vs tile count (2 panels)
- `mbias_errors.png` — error shrinkage vs tile count (2 panels)

---

## Quick Start

### Prerequisites
- container: `/n17data/mkilbing/sp_validation_im_sims.sif`
- Pipeline scripts: `/path/to/shapepipe/scripts/image_sims_pipeline/`

### Setup (one time)

```bash
# Copy template config to run directory and rename to config.yaml
cp /path/to/shapepipe/scripts/image_sims_pipeline/config.yaml.template <run_dir>/config.yaml

# Edit config.yaml with your tile IDs, grid number, etc.
```

### Run (fully automatic, snakemake)

To get help on how to run snakemake, and the different targets and options, type
```bash
python ~/shapepipe/scripts/image_sims_pipeline/info.py -h
```

To monitor the status of the run and completed number of tiles, type
```bash
python ~/shapepipe/scripts/image_sims_pipeline/info.py -m [-v]
```

First, for convenience define the short cut
```bash
alias sn_ims="snakemake -s ~/shapepipe/scripts/image_sims_pipeline/Snakefile --configfile config.yaml -j $OMP_NUM_THREADS'
```

To go step by step (one target at a time), do

```bash
# Initialise subdirs (one for each sheared verions)
sn_ims init_all

# Optional: scan input weight images for tile coverage. Tiles with a low
fraction of non-zero weight pixels (<50%, can be set by the user) has a small
number of detected galaxies. In some cases this number can be zero which would
lead to a failure of the ShapePipe module vignet_maker_runner.

This call creates a tile-by-tile coverage report and adds the below-threshold
tiles to an exclusion list in config.yaml.

# Run shapepipe
sn_ims pipeline_all
```

```bash

# Full pipeline to calibrated catalogues
snakemake -s ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/Snakefile \
  --configfile config.yaml -j 5 calibrate_all

# M-bias computation with convergence tracking
snakemake -s ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/Snakefile \
  --configfile config.yaml -j 1 mbias

# Monitor progress
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/info.py -m -v

# Incremental m-bias as tiles finish
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/monitor_mbias.py -v
```

### Run (non-automatic)

Activate Docker container with
```bash
app_gen --bind /n09data /n17data/mkilbing/shapepipe_im_sims.sif
```

Get config file
```bash
cp /path/to/shapepipe/scripts/image_sims_pipeline/config.yaml.template ./config.yaml
```

Edit config file.

Initialise subdir:
```bash
init_run_v2.0.py [-s <subdir>]
```
with e.g. `subdir=1p2z_grid_1`.

Start pipeline runs

```bash
cd <subdir>
run_job_sp_canfar_v2.0.bash -t image_sims -j 1
```

This will run the first job (`-j 1`) for all tiles found in `tile_numbers.txt`.

To get help on the options, use the flag `-h`.

To test it one a single tile for testing, add the option `-e ID`, e.g. `-e 272.281`.

To run all pipeline jobs, either call `run_job_sp_canfar_v2.0.bash` for each subsequent job,
e.g. `-j 1`, `-j 2`, `-j 4`, or group them by (bit-)adding the job numbers. The entire pipeline
will be run with `-j 4095`.


## Configuration

### Snakemake Config: `config.yaml`

```yaml
# Run directory
base:    /n17data/mkilbing/astro/Runs/shapepipe/CFIS/v2.0/image_sims

# Tile IDs (all tiles in the survey, or subset for testing)
tile_IDs: ["254.286", "254.287"]

# catalogue subversion (mask config: mask_v1.X.{sample}_im_sim.yaml)
sample:  9

# Simulation type (grid or others)
type:    grid

# Grid number (e.g., grid_2)
num:     2

# SMP batch size per ShapePipe job (-1 = use container default)
n_smp:   -1
```

### M-bias Calibration Settings

M-bias computation parameters are hard-coded in the Snakefile `mbias` rule:

```python
shear_amplitude:     0.02      # Input shear amplitude (absolute value)
match_radius_deg:    0.0002    # Position matching radius for cross-matching
w_col:               w_des     # Weight column name
n_bootstrap:         500       # Bootstrap resamples for error estimation
catalog_name:        shape_catalog_cut_ngmix.fits
```

To modify these parameters for a different run, edit the `mbias` rule in the Snakefile.

---

## Pipeline Stages

The pipeline runs in order, with each stage dependent on the previous:

### 1. `init_all` — Initialize run directories

Creates per-sim directories, parameter files (`params.py`), and mask configs.

**Outputs:** `{grids}/{sim}/params.py`, `{grids}/{sim}/config_mask.yaml`

### 2. `pipeline_all` — Run ShapePipe on all tiles

Executes ShapePipe's full job sequence (bits 1→2→4→...→2048) for each tile in each simulation.

**Note:** `-j 5` means 5 concurrent Snakemake jobs (tile×sim pairs run sequentially within each job).

**Outputs:** ShapePipe catalogues, logs at `{grids}/{sim}/logs/log_job_{tile}_{bit}.txt`

**Time:** ~hours per grid depending on tile count and cluster load.

### 3. `merge_all` — Merge ShapePipe tile outputs

Combines per-tile catalogs into a single HDF5 file per simulation using `create_final_cat.py`.

**Outputs:** `{grids}/{sim}/final_cat_{sim}.hdf5` (HDF5 with tile counts as attributes)

### 4. `extract_all` — Extract comprehensive catalogues

Reads merged HDF5 and extracts shape information, PSF quantities, and pre-calibration columns using `extract_info.py`.

**Outputs:** `{grids}/{sim}/shape_catalog_comprehensive_ngmix.hdf5` (with `n_tiles` attribute in header)

### 5. `calibrate_all` — Apply cuts and calibrate

Applies selection masks (flags, magnitude, signal-to-noise) and computes shear calibration (m/c per object) using `calibrate_comprehensive_cat.py`.

**Outputs:** `{grids}/{sim}/shape_catalog_cut_ngmix.fits` (final catalogue, ready for m-bias)

### 6. `mbias` — Compute m-bias with convergence tracking

Computes m₁, m₂, c₁, c₂ from the five shear pairs. With `--cumulative`, tracks convergence as tile counts grow.

**Outputs:**
- `results/m_bias_results.yaml` — final results (YAML)
- `results/m_bias_results.txt` — final results (human-readable text)
- `results/mbias_cumulative.yaml` — convergence history
- `results/mbias_convergence.png` — m/c vs n_tiles (2 panels)
- `results/mbias_errors.png` — error convergence (2 panels)

---

## Monitoring Convergence

### Live Status: `info.py`

Monitor pipeline progress with a status table:

```bash
cd /path/to/run/image_sims
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/info.py -m -v
```

**Status table (example output):**

| Simulation | #Jobs | Job Bits | Merge | Extract | Calibrate |
|-----------|-------|----------|-------|---------|-----------|
| 1m2z_grid_2 | 12/12 | ✓✓✓✓✓✓✓✓✓✓✓✓ | 2 | ✓ | ✓ |
| 1p2z_grid_2 | 12/12 | ✓✓✓✓✓✓✓✓✓✓✓✓ | 2 | ✓ | ✓ |
| 1z2m_grid_2 | 12/12 | ✓✓✓✓✓✓✓✓✓✓✓✓ | 2 | ✓ | ✓ |
| 1z2p_grid_2 | 12/12 | ✓✓✓✓✓✓✓✓✓✓✓✓ | 2 | ✓ | ✓ |
| 1z2z_grid_2 | 12/12 | ✓✓✓✓✓✓✓✓✓✓✓✓ | 2 | ✓ | ✓ |

**Column meanings:**
- **#Jobs** — ShapePipe job completion count (max 12 bits)
- **Job Bits** — per-bit status (✓ = all tiles done, · = incomplete)
- **Merge** — tile count in final HDF5 catalogue
- **Extract** — comprehensive catalog status
- **Calibrate** — calibrated catalogue status

### Incremental M-bias: `monitor_mbias.py`

Recompute m-bias each time new tiles finish all ShapePipe jobs. Automatically detects completed tiles and updates convergence tracking.

```bash
cd /path/to/run/image_sims
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/monitor_mbias.py -v
```

Run this repeatedly as tiles complete to watch m-bias converge.

---

## Results Interpretation

### M-Bias Results Files

**YAML format** (`m_bias_results.yaml`):
```yaml
m1: -0.7909173336321427
m1_err: 0.27407603929732716
c1: -0.0030413689177988708
c1_err: 0.005899515490890933
m2: -1.350016719584532
m2_err: 0.20884008300375373
c2: -0.004950048084875509
c2_err: 0.004208875641566414
```

**Text format** (`m_bias_results.txt`):
```
Multiplicative and additive shear bias from image simulations
============================================================

m1 = -0.790917 ± 0.274076
c1 = -0.003041 ± 0.005900

m2 = -1.350017 ± 0.208840
c2 = -0.004950 ± 0.004209

Errors computed via bootstrap resampling (n=500 resamples)
```

### Convergence History (`mbias_cumulative.yaml`)

Tracks m/c evolution as tile count grows:

```yaml
'2':
  c1: -0.0030413689177988708
  c1_err: 0.005899515490890933
  m1: -0.7909173336321427
  m1_err: 0.27407603929732716
  # ... (c2, m2, errors)
'4':
  c1: -0.0025...
  # ... (more tiles)
```

### Understanding Error Bars

Errors are computed via **bootstrap resampling**:

1. Draw N=500 random resamples (with replacement) from measured galaxy ellipticities
2. Recompute m and c for each resample
3. Error = standard deviation of the bootstrap distribution

This captures:
- **Photometric noise** — measurement uncertainties per galaxy
- **Cosmic variance** — shape correlations across the sky
- **Calibration uncertainties** — from shear responsivity scatter

**Error shrinkage:** Error ∝ 1/√(n_tiles). More tiles → tighter constraints.

### Convergence Plots

**`mbias_convergence.png`** — m and c with error bars vs n_tiles
- Left panel: multiplicative bias (m₁, m₂) with error bars
- Right panel: additive bias (c₁, c₂) with error bars
- Shows systematic trends and statistical uncertainties

**`mbias_errors.png`** — error shrinkage vs n_tiles
- Left panel: m error bars (m₁, m₂) only
- Right panel: c error bars (c₁, c₂) only
- Shows constraint tightening as data accumulates
- Useful for deciding when m-bias is "converged" (errors small enough for science)

---

## File Structure

```
/n17data/mkilbing/astro/Runs/shapepipe/CFIS/v2.0/image_sims/
├── config.yaml                     # Run configuration
├── .snakemake/                     # Snakemake metadata
├── logs/                           # Snakemake logs
├── grids/
│   ├── {sim}_grid_2/
│   │   ├── logs/                   # ShapePipe job logs
│   │   ├── tiles/                  # Per-tile ShapePipe outputs
│   │   ├── final_cat_{sim}.hdf5    # Merged catalogue
│   │   ├── shape_catalog_comprehensive_ngmix.hdf5
│   │   └── shape_catalog_cut_ngmix.fits
│   └── results/
│       ├── m_bias_results.yaml
│       ├── m_bias_results.txt
│       ├── mbias_cumulative.yaml
│       ├── mbias_convergence.png
│       └── mbias_errors.png
└── monitoring/                     # Incremental m-bias workspace
    └── {sim}/
        ├── final_cat_{sim}.hdf5
        ├── shape_catalog_comprehensive_ngmix.hdf5
        └── shape_catalog_cut_ngmix.fits
```

---

## References

- **Pipeline code:** `~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/`
- **ShapePipe:** `~/astro/repositories/github/shapepipe/`
- **SP Validation:** `~/astro/repositories/github/sp_validation/`

---

*Last updated: 2026-06-26*
