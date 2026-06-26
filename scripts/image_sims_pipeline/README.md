# Image Simulations Calibration Pipeline

Snakemake workflow for deriving shear m-bias from image simulations.

## Setup

### 1. Create run directory and copy config

```bash
mkdir -p /path/to/run/image_sims
cd /path/to/run/image_sims

# Copy template config and edit for your run
cp ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/config.yaml.template \
   config.yaml

# Edit config.yaml with your tile IDs, grid number, etc.
```

### 2. Run Snakemake from the run directory

```bash
cd /path/to/run/image_sims

# Full pipeline
snakemake -s ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/Snakefile \
  --configfile config.yaml -j 5 calibrate_all

# M-bias computation
snakemake -s ~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/Snakefile \
  --configfile config.yaml -j 1 mbias
```

## Scripts

**Snakefile** — Main pipeline definition  
**info.py** — Monitor pipeline progress with status table  
**monitor_mbias.py** — Incrementally compute m-bias as tiles finish  
**config.yaml.template** — Configuration template (copy and edit for each run)

## Usage Examples

### Monitor progress (from run directory)

```bash
cd /path/to/run/image_sims
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/info.py -m -v
```

### Incremental m-bias computation

```bash
cd /path/to/run/image_sims
~/astro/repositories/github/shapepipe/scripts/image_sims_pipeline/monitor_mbias.py -v
```

Run this repeatedly while `pipeline_all` is in progress to watch m-bias converge.

## Output Structure

```
/path/to/run/image_sims/
├── config.yaml
├── .snakemake/
├── logs/                          # Snakemake logs
├── grids/
│   ├── 1m2z_grid_2/
│   ├── ...
│   └── results/
│       ├── m_bias_results.yaml
│       ├── m_bias_results.txt
│       ├── mbias_cumulative.yaml
│       ├── mbias_convergence.png
│       └── mbias_errors.png
└── monitoring/                    # Temp workspace for incremental m-bias
```

## Documentation

Full documentation: `~/astro/repositories/github/shapepipe/docs/source/image_sims_calibration.md`
