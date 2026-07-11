# Installation

```{attention}
ShapePipe was not designed to be a stand-alone Python library. Instead users
are expected to install the full ShapePipe environment on the system(s) where
data should be processed.
```

## Container Installation (Recommended)

The easiest way to install ShapePipe is via a container. Docker images are automatically built and pushed to the [Github Container Registry (GHCR)](ghcr.io/cosmostat/shapepipe) on every push. Two image targets are published per branch — the default tag is the rich `dev` image; `<tag>-runtime` is a slim variant for batch jobs. See [Container Workflow](container.md) for the full rationale and the relationship between `pyproject.toml`, `uv.lock`, and the `Dockerfile`.

We recommend running the image with **Apptainer** (formerly Singularity) which is installed on most HPC clusters. To simply run the image, use the following command:

```bash
# build writeable "sandbox" container in the current directory
# ./shapepipe will be a directory that functions like a vm
apptainer build --sandbox shapepipe docker://ghcr.io/cosmostat/shapepipe:develop

# open a shell in the container and run the example
apptainer shell --writable shapepipe
cd /app && shapepipe_run -c /app/example/config.ini
```

You can also run the image with **Docker**:

```bash
docker run --rm -it ghcr.io/cosmostat/shapepipe:develop shapepipe_run -c /app/example/config.ini
```

For canfar batch jobs or downstream images, the slim runtime tag is preferred:

```bash
docker pull ghcr.io/cosmostat/shapepipe:develop-runtime
```

```{attention}
We do not currently build images for Apple Silicon/amr64; however the amd64 images should work on these systems, albeit with reduced performance.
```

The image bundles the astromatic binaries (`source-extractor`, `psfex`,
`weightwatcher`), MPI (`mpi4py` + OpenMPI), and every Python dependency, so
there is nothing else to install or build. To process data on a cluster with
MPI, run the pipeline through Apptainer the same way you would any MPI job.

## Installing the ShapePipe Library Only

The ShapePipe library, i.e. the core Python package *without* the system
executables (`source-extractor`, `psfex`, …) or the bundled examples, can be
installed with `pip`. This is enough to import `shapepipe` or to develop
against the package, but not to run the full pipeline.

After cloning the repository:

```bash
pip install .
```

Without cloning the repository:

```bash
pip install git+https://github.com/CosmoStat/shapepipe.git
```

For a development checkout with the test, lint, and doc tools, install the
`dev` extra (ideally into a fresh virtual environment, e.g. one managed by
[uv](https://docs.astral.sh/uv/)):

```bash
pip install -e ".[dev]"
```

## Troubleshooting

If you encounter any problems installing ShapePipe please
[open an issue](https://github.com/CosmoStat/shapepipe/issues/new/choose) and
we will do our best to help you.
