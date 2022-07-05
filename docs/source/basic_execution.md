# Basic Execution

ShapePipe pipelines are launched and managed via the `shapepipe_run` script.

A list of command line arguments can be displayed using the `--help`
option:

```{seealso}
:class: margin
The `shapepipe` environment will need to be built and activated in order to run this script (see [Installation](installation.md)).
```
```bash
shapepipe_run --help
```

The options for defining a pipeline are managed via a
[configuration file](configuration.md).

## Running ShapePipe with Joblib

By default ShapePipe is run using
[Joblib](https://joblib.readthedocs.io/en/latest/) to manage parallel
processes. The `shapepipe_run` script can run as follows

```bash
shapepipe_run
```

which will by default look for a file called `config.ini` in the current
directory. To specify the path to a given configuration file (with any name)
you would run

```bash
shapepipe_run -c <PATH TO CONFIG FILE>
```

## Running the Pipeline with MPI

ShapePipe can also use [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
for managing parallel processes on clusters with multiple nodes.
The `shapepipe_run` script can be run with MPI as follows

```bash
mpiexec -n <NUMBER OF CORES> shapepipe_run
```

where `<NUMBER OF CORES>` is the number of cores to allocate to the run.
