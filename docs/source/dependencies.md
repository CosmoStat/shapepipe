# Dependencies

ShapePipe's dependencies are installed automatically — the recommended path is
the [container](container.md), which bundles all of them (see also
[Installation](installation.md)). `pyproject.toml` is the authoritative list: it
declares the **abstract minimum versions** ShapePipe is compatible with, while
`uv.lock` pins the **exact** versions that ship in a given build. The tables
below describe the main packages; for the complete, current set always refer to
`pyproject.toml`.

## Python Dependencies

Scientific stack:

| Package | References |
|---------|------------|
| [Astropy](https://www.astropy.org/) | {cite:p}`astropy:2013,astropy:2018` |
| [GalSim](https://github.com/GalSim-developers/GalSim) | {cite:p}`rowe:15` |
| [ngmix](https://github.com/aguinot/ngmix) | {cite:p}`sheldon:15` |
| [MCCD](https://github.com/CosmoStat/mccd) | {cite:p}`liaudat:21` |
| [ModOpt](https://cea-cosmic.github.io/ModOpt/) | {cite:p}`farrens:20` |
| [python-pysap](https://github.com/CEA-COSMIC/pysap) | |
| [Numpy](https://numpy.org/) | {cite:p}`harris:20` |
| [Numba](https://numba.pydata.org/) | |
| [Pandas](https://pandas.pydata.org/) | {cite:p}`pandas:10,pandas:20` |
| [Matplotlib](https://matplotlib.org/) | {cite:p}`hunter:07` |
| [Joblib](https://joblib.readthedocs.io/en/latest/) | {cite:p}`joblib:20` |
| [mpi4py](https://mpi4py.readthedocs.io/en/stable/) | {cite:p}`dalcin:05,dalcin:08,dalcin:11` |
| [reproject](https://reproject.readthedocs.io/) | |
| [h5py](https://www.h5py.org/) | |
| [sf_tools](https://github.com/sfarrens/sf_tools) | |

Data access &amp; infrastructure (CANFAR / UNIONS):

| Package | Purpose |
|---------|---------|
| [vos](https://github.com/opencadc/vostools) | CADC / CANFAR VOSpace access |
| canfar | CANFAR session monitoring |
| [astroquery](https://astroquery.readthedocs.io/) | external catalogue queries |
| [cs_util](https://github.com/CosmoStat/cs_util) | shared CosmoStat utilities |
| [sqlitedict](https://github.com/RaRe-Technologies/sqlitedict) | on-disk pipeline state |

```{note}
`ngmix` is pinned to the
[`aguinot/ngmix@stable_version`](https://github.com/aguinot/ngmix/tree/stable_version)
fork until the fixes land upstream — do not modernise that line (see the note in
`pyproject.toml`).
```

## System Dependencies

The container also provides the non-Python tools ShapePipe calls, all from Debian
packages (no source builds), plus the MPI stack:

| Package | References |
|---------|------------|
| [Source Extractor](https://www.astromatic.net/software/sextractor/) | {cite:p}`bertin:96` |
| [PSFEx](https://www.astromatic.net/software/psfex/) | {cite:p}`bertin:11` |
| [WeightWatcher](https://www.astromatic.net/software/weightwatcher/) | {cite:p}`marmo:08` |
| OpenMPI (5.0.x) | |

Python dependencies themselves are managed with [uv](https://docs.astral.sh/uv/);
see [Container Workflow](container.md) for how `pyproject.toml`, `uv.lock`, and
the `Dockerfile` fit together.
