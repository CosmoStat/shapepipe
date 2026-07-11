# Basic Execution

ShapePipe pipelines are launched and managed via the `shapepipe_run` script.

A list of command line arguments can be displayed using the `--help`
option:

```{seealso}
:class: margin
ShapePipe runs inside its container — there is no environment to activate. See [Installation](installation.md).
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
to spread work across multiple nodes of a cluster. Set `MODE = mpi` in the
`[EXECUTION]` section of the config and launch with an MPI runner:

```bash
mpiexec -n <NUMBER OF RANKS> shapepipe_run -c <PATH TO CONFIG FILE>
```

where `<NUMBER OF RANKS>` is the number of MPI processes to start.

### Through the container (the supported way on a cluster)

On a cluster you run ShapePipe from the published image as a standard Apptainer
*hybrid* MPI job: the **host** `mpirun`/`mpiexec` launches one container rank per
slot, and the OpenMPI bundled in the image wires the ranks together.

```bash
# one-time: pull the runtime image
apptainer pull shapepipe.sif docker://ghcr.io/cosmostat/shapepipe:develop-runtime

# load a host MPI in the same family as the image's OpenMPI (5.0.x), then launch
module load openmpi
mpirun -n <NUMBER OF RANKS> \
    apptainer exec --bind "$PWD:$PWD" shapepipe.sif \
    shapepipe_run -c <PATH TO CONFIG FILE>
```

The image ships **OpenMPI 5.0.x** so that its PMIx matches modern cluster
launchers. The host and container MPI must be compatible: if you see *N* copies
of `rank 0 of 1` instead of one *N*-rank job, load a host OpenMPI in the 5.0.x
family. See `example/pbs/candide_mpi.sh` for a complete SLURM batch script.
