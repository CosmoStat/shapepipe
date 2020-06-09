[Home](./shapepipe.md) | [Environments](./environment.md)

# CCIN2P3 Set Up

## Contents

1. [Introduction](#Introduction)
1. [Installation](#Installation)
1. [Execution](#Execution)
1. [Troubleshooting](#Troubleshooting)

## Introduction

[CC-IN2P3](https://doc.cc.in2p3.fr/en-index.html) is the computing centre (CC) of the French National Institute of Nuclear and Particle Physics (INP3) based in Lyon.


### CCIN2P3 Account

To get an account on CCIN2P3 you need to follow [their instructions](https://doc.cc.in2p3.fr/en/Getting-started/access.html). In summary:
- Fill in the [application form](http://cctools.in2p3.fr/cclogon/) and print the resulting PDF
- Get this document signed by your local "czar". For CEA this should be [Georgette Zoulikha](mailto:zou.georgette@cea.fr).
- [Submit a ticket](https://cc-usersupport.in2p3.fr/otrs/customer.pl) with a scan of the signed document.
- If required indicate `Euclid` as project.

### SSH

Once you have an account on CCIN2P3 you can connect via SSH as follows:

```bash
$ ssh <username>@cca.in2p3.fr
```

## Installation

The CCIN2P3 system uses an internal system called `ccenv` to manage various software packages. You can view the packages currently available on the system by running:

```bash
$ ccenv --list
```

ShapePipe requires `conda`, which on CCIN2P3 is provided via `anaconda`. To load this package simply run:

```bash
$ ccenv anaconda
```

> You can add this command to your `.bash_profile` to ensure that this module is available when you log in.

### With MPI

To install ShapePipe with MPI enabled on CCIN2P3 you also need to load the `openmpi` module. To do so run:

```bash
$ ccenv openmpi
```

You can also specify a specific version of OpenMPI to use.

```bash
$ ccenv openmpi <VERSION>
```

Then you need to identify the root directory of the OpenMPI installation. A easy way to get this information is by running:

```bash
$ which mpiexec
```

which should reveal something like `/pbs/software/centos-7-x86_64/openmpi/<VERSION>`. Provide this path to the `mpi-root` option of the installation script as follows:

```bash
$ ./shapepipe_install --mpi-root=/pbs/software/centos-7-x86_64/openmpi/<VERSION>
```

> Be sure to check the output of the **Installing MPI** section, as the final check only tests if the `mpiexec` command is available on the system.

You can rebuild the MPI component at any time by doing the following:

```bash
$ pip uninstall mpi4py
$ ./install_shapepipe --no-env --no-exe --mpi-root=/pbs/software/centos-7-x86_64/openmpi/<VERSION>
```

### Without MPI

To install ShapePipe without MPI enabled simply pass the `no-mpi` option to the installation script as follows:

```bash
$ ./shapepipe_install --no-mpi
```

## Execution

CCIN2P3 uses [UGE](https://en.wikipedia.org/wiki/Univa_Grid_Engine) for handling distributed jobs.

UGE uses standard [Portable Batch System (PBS) commands](https://www.cqu.edu.au/eresearch/high-performance-computing/hpc-user-guides-and-faqs/pbs-commands) such as:

- `qsub` - To submit jobs to the queue.
- `qstat` - To check on the status of jobs in the queue.
- `qdel` - To kill jobs in the queue.

Jobs should be submitted as bash scripts. *e.g.*:

```bash
$ qsub cc_smp.sh
```

In this script you can specify:

- The number of nodes to use for SMP (*e.g.* `#$ -pe multicores 4`)
- The number of nodes to use for MPI (*e.g.* `#$ -pe openmpi 8`)
- The maximum computing time for your script (*e.g.* `#$ -l h_cpu=10:00:00`)

### Example SMP Script

[`cc_smp.sh`](../../example/pbs/cc_smp.sh)

```bash
#!/bin/bash

##########################
# SMP Script for ccin2p3 #
##########################

# Receive email when job finishes or aborts
#$ -M <name>@cea.fr
#$ -m bea
# Set a name for the job
#$ -N shapepipe_smp
# Set a group for the job
#$ -P P_euclid_sci
# Join output and errors in one file
#$ -j y
# Set maximum computing time (e.g. 5min)
#$ -l h_cpu=00:05:00
# Request muliprocessing resources
#$ -l os=cl7
# Request number of cores
#$ -pe multicores 4

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Activate conda environment
ccenv anaconda
source activate $SPENV

# Run ShapePipe using full paths to executables
$SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_smp.ini

# Return exit code
exit 0
```

> Make sure the number of nodes requested matches the `SMP_BATCH_SIZE` in the config file.

### Example MPI Script

[`cc_mpi.sh`](../../example/pbs/cc_mpi.sh)

```bash
#!/bin/bash

##########################
# MPI Script for ccin2p3 #
##########################

# Receive email when job finishes or aborts
#$ -M <name>@cea.fr
#$ -m bea
# Set a name for the job
#$ -N shapepipe_mpi
# Set a group for the job
#$ -P P_euclid_sci
# Join output and errors in one file
#$ -j y
# Set maximum computing time (e.g. 5min)
#$ -l h_cpu=00:05:00
# Request muliprocessing resources
#$ -l os=cl7
# Request number of cores
#$ -pe openmpi 4

# Full path to environment
export SPENV="$HOME/.conda/envs/shapepipe"
export SPDIR="$HOME/shapepipe"

# Activate conda environment
ccenv anaconda
ccenv openmpi
source activate $SPENV

# Run ShapePipe
$SPENV/bin/mpiexec -n $NSLOTS $SPENV/bin/shapepipe_run -c $SPDIR/example/pbs/config_mpi.ini

# Return exit code
exit 0
```

## Troubleshooting

The Euclid-SDC slack [channel](#euclid-sdc-france.slack.com) is very helpful in case of problems of job submission, environment set up, library and program versions etc. You can write a PM directly at some of the Euclid sys ads, Quentin Le Boulc'h and Gabriele Mainetti (in English or French). They are quick to help with any issue.
