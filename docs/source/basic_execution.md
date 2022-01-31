# Basic Execution

## Running the ShapePipe with SMP

The pipeline can be run with SMP as follows:

```bash
shapepipe_run
```

A list of command line arguments can be displayed using the `--help`
option:

```bash
shapepipe_run --help
```

## Running the Pipeline with MPI

The pipeline can be run with MPI as follows:

```bash
mpiexec -n <number_of_cores> shapepipe_run
```

where `<number_of_cores>` is the number of cores to allocate to the run.
