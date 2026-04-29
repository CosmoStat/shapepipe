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

## Conda Installation (Deprecated)

```{tip}
:class: margin
Check out [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for a
light weight and easy installation of Conda.
```

The standard installation of ShapePipe manages [dependencies](dependencies.md)
and scripts using a [Conda](https://docs.conda.io/en/latest/) environment.
Therefore, to follow the standard installation Conda must be available on the
system.

The ShapePipe package should first be cloned (or downloaded) from the
[GitHub repository](https://github.com/CosmoStat/shapepipe).

```{note}
:class: margin
Developers should simply clone the repository as usual.
```
```bash
git clone -b <VERSION> --depth 1 git@github.com:CosmoStat/shapepipe.git
cd shapepipe
```

where `<VERSION>` is a
[tagged release](https://github.com/CosmoStat/shapepipe/releases) of ShapePipe
(e.g. `v1.0.1`). It is recommend to use the
[latest release](https://github.com/CosmoStat/shapepipe/releases/latest)
unless you want to reproduce an older set of results.

Then, the entire ShapePipe environment, including dependencies, can be built
using the `install_shapepipe` script as follows.

```bash
./install_shapepipe
```

The `install_shapepipe` script will create the recommended Conda environment
along with all of the required core and module dependencies. The script also
provides a checklist for the success or failure of installing each of the
dependencies.

Once the installation is complete the `shapepipe` environment needs to be
activated.

```bash
conda activate shapepipe
```

A list installation options can be seen using the `--help` option.

```bash
./install_shapepipe --help
```

In particular, the `--develop` flag can be used to install additional tools
for developers.

## MPI Installation

The standard installation of ShapePipe will install and enable MPI on a given
node (i.e. a given machine). However, `mpi4py` needs to be built from source
in order to take advantage of a preinstalled MPI distribution on the compute
nodes of a cluster.

This can be done as follows

```bash
./install_shapepipe --mpi-root=<PATH TO MPI>
```

where `<PATH TO MPI>` is the full path to MPI root directory (i.e. where the
`bib`, `include` and `lib` directories can be found). If the installation is
successful, ShapePipe will be able to submit jobs to all of the nodes in the
cluster.

## Uninstalling ShapePipe

The `--uninstall` flag can be passed to `install_shapepipe` to remove the
entire ShapePipe environment.

```bash
./install_shapepipe --uninstall
```

## Installing the ShapePipe Library Only

The ShapePipe library, i.e. the core package not including module
dependencies, can be installed in the following ways.

After cloning the repository.

```{warning}
:class: margin
Note, this method will not include any system executables or examples.
```
```bash
pip install .
```

Without cloning the repository.

```bash
pip install git+https://github.com/CosmoStat/shapepipe.git
```

## Troubleshooting

If you encounter any problems installing ShapePipe please
[open an issue](https://github.com/CosmoStat/shapepipe/issues/new/choose) and
we will do our best to help you.
