[Home](./shapepipe.md) | [Environments](./environment.md)

# CANDIDE Set Up

> Environment Status Notes  
> - Website: https://candideusers.calet.org/
> - No internet access on compute nodes, see [tutorial](https://github.com/CosmoStat/shapepipe/blob/master/docs/wiki/tutorial/pipeline_tutorial.md#mask-images) for how to manage `mask_runner`
> - Current stable OpenMPI version: `4.0.2`

## Contents

1. [Introduction](#Introduction)
1. [Installation](#Installation)
1. [Execution](#Execution)
1. [Troubleshooting](#Troubleshooting)

## Introduction

The [CANDIDE cluster](https://candideusers.calet.org/) is hosted and maintained at the Institut d’Astrophysique de Paris by Stephane Rouberol.

### CANDIDE Account

To request and account on CANDIDE send an email to [Henry Joy McCracken](mailto:hjmcc@iap.fr) and [Stephane Rouberol](mailto:rouberol@iap.fr) at IAP with a short description of what you want to do and with whom you work.

### SSH

Once you have an account on CANDIDE you can connect via SSH as follows:

```bash
$ ping -c 1 -s 999 candide.iap.fr; ssh <mylogin>@candide.iap.fr
```

## Installation

The CANDIDE system uses [Environment Modules](https://modules.readthedocs.io/en/latest/) to manage various software packages. You can view the modules currently available on the system by running:

```bash
$ module avail
```

ShapePipe requires `conda`, which on CANDIDE is provided via `intelpython/3`. To load the newest version of this package simply run:

```bash
$ module load intelpython/3-2020.1
```

For the installation of `astromatic` software (`SExtractor`, `PSFEx`) the `BLAS` library is required. This is made available by loading the Intel MLK module,
```bash
module load intel/19.0
```

> You can add these commands to your `.bash_profile` to ensure that this module is available when you log in.

You can list the modules already loaded by running:

```bash
$ module list
```

### With MPI

To install ShapePipe with MPI enabled on CANDIDE you also need to load the `openmpi` module. To do so run:

```bash
$ module load openmpi
```

You can also specify a specific version of OpenMPI to use.

```bash
$ module load openmpi/<VERSION>
```

Then you need to identify the root directory of the OpenMPI installation. A easy way to get this information is by running:

```bash
$ module show openmpi
```

which should reveal something like `/softs/openmpi/<VERSION>-torque-CentOS7`. Provide this path to the `mpi-root` option of the installation script as follows:

```bash
$ ./install_shapepipe --mpi-root=/softs/openmpi/<VERSION>-torque-CentOS7
```

> Be sure to check the output of the **Installing MPI** section, as the final check only tests if the `mpiexec` command is available on the system.

You can rebuild the MPI component at any time by doing the following:

```bash
$ pip uninstall mpi4py
$ ./install_shapepipe --no-env --no-exe --mpi-root=/softs/openmpi/<VERSION>-torque-CentOS7
```

### Without MPI

To install ShapePipe without MPI enabled simply pass the `no-mpi` option to the installation script as follows:

```bash
$ ./install_shapepipe --no-mpi
```

## Execution

CANDIDE uses [TORQUE](https://en.wikipedia.org/wiki/TORQUE) for handling distributed jobs.

TORQUE uses standard [Portable Batch System (PBS) commands](https://www.cqu.edu.au/eresearch/high-performance-computing/hpc-user-guides-and-faqs/pbs-commands) such as:

- `qsub` - To submit jobs to the queue.
- `qstat` - To check on the status of jobs in the queue.
- `qdel` - To kill jobs in the queue.

Additionally, the availability of compute nodes can be seen using the command

```bash
$ cnodes
```

Jobs should be submitted as bash scripts. *e.g.*:

```bash
$ qsub candide_smp.sh
```

In this script you can specify:

- The number of nodes to use (*e.g.* `#PBS -l nodes=10`)
- A specific machine to use with a given number of cores (*e.g.* `#PBS -l nodes=n04:ppn=10`)
- The maximum computing time for your script (*e.g.* `#PBS -l walltime=10:00:00`)

### Example SMP Script

[`candide_smp.sh`](../../example/pbs/candide_smp.sh)

```bash
#!/bin/bash

##########################
# SMP Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_smp
# Join output and errors in one file
#PBS -j oe
# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00
# Request number of cores
#PBS -l nodes=4

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Activate conda environment
module load intelpython/3
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_smp.ini

# Return exit code
exit 0
```

> Make sure the number of nodes requested matches the `SMP_BATCH_SIZE` in the config file.

### Example MPI Script

[`candide_mpi.sh`](../../example/pbs/candide_mpi.sh)

```bash
#!/bin/bash

##########################
# MPI Script for CANDIDE #
##########################

# Receive email when job finishes or aborts
#PBS -M <name>@cea.fr
#PBS -m ea
# Set a name for the job
#PBS -N shapepipe_mpi
# Join output and errors in one file
#PBS -j oe
# Set maximum computing time (e.g. 5min)
#PBS -l walltime=00:05:00
# Request number of cores (e.g. 4 from 2 different machines)
#PBS -l nodes=2:ppn=2
# Allocate total number of cores to variable NSLOTS
NSLOTS=`cat $PBS_NODEFILE | wc -l`

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Load moudules and activate conda environment
module load intelpython/3
module load openmpi/4.0.2
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/mpiexec -n $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
```

## Troubleshooting

### OpenBLAS

If you get the following error

```
error while loading shared libraries: libopenblas.so.0: cannot open shared object file: No such file or directory
```

simply run

```bash
$ export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
```

> You can add the command to your `.bash_profile`.
