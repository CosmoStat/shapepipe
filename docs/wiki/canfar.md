[Home](./shapepipe.md) | [Environments](./environment.md)

Authors:

Sam Farrens <samuel.farrens@cea.fr>

Martin Kilbinger <martin.kilbinger@cea.fr>

Axel Guinot <axel.guinot@cea.fr>

# CANFAR Set Up

Here are some instructions on how to set up a VM on CANFAR to run the pipeline.

## Contents

1. [Current Set Up](#Current-Set-Up)
1. [Virtual Machine](#Virtual-Machine)
1. [Batch System](#Batch-System)
1. [Troubleshooting](#Troubleshooting)
1. [Interactive Mode](#Interactive-Mode)
1. [Running a CFIS job](#Running-A-CFIS-Job)

## Current Set Up

### Instances, snapshots, and VMs

An `instance` contains the base setting, which is used to create VMs, where
in turn jobs are submitted and run.

A `snapshot` or `image` is created from an instance, and contains a frozen version of
the instance. An instance can have a number of snapshots, e.g. to use
different version of some code or library.

A  `Virtual Machine` (VM) is a copy of a snapshot on which a job is run.


### Available Instances

Once logged in to computecanada, the user can see all instances
for the project(s) they are registered on the page
https://arbutus-canfar.cloud.computecanada.ca/project/instances.

For `ShapePipe` we use the following instance:
- `ShapePipe2`:
- 90 GB RAM
- 20 GB Disk
- 8 CPUs
- Flavour: c8-90gb-186
- IP Address: 206.12.92.159

### Available Snapshots

All snapshots are listed here:
https://arbutus-canfar.cloud.computecanada.ca/project/images
Snapshots created from the `ShapePipe2` instance are typically called
`ShapePipe2-mk-<date>`.

Note that if `vos` down- or up-loads to the `canfar` storage are performed by a job,
this requires a `cadc` certificate, which is valid for 10 days. In that case, a snapshot
will not function anymore later than 10 days after its creation.

## Virtual Machine

The virtual machine (VM) is a space where we can install software under a given Linux distribution with given CPU, RAM and storage limits. Once we are happy with a given set-up we can freeze these conditions (*i.e.* all the software versions *etc.* currently installed) by creating a *snapshot* that acts like a container for the VM. Jobs can then be submitted through the batch system using a given snapshot.

In general the processing should be done through the batch system and not run directly on the VM.
However, tests of the VM and the code to be submitted can be run in (interactive mode)[interactive(-mode] on the VM.

The following steps show how to set up a VM.

1. Create a VM:

    Follow the instructions on [CANFAR quick start](https://www.canfar.net/en/docs/quick_start/).

    VMs can be managed on [OpenStack](https://arbutus-canfar.cloud.computecanada.ca/).

    > Note: An IP address has to be assigned to the VM in order to be able to log in and there are a limited number of IPs per workspace.

2. SSH to VM:

    Run the following command to connect to a given VM:

    ```bash
    ssh ubuntu@IP_ADDRESS
    ```

    This will connect you to a generic *ubuntu* user space, shared between all users. Once connected software *etc.* can be installed and tested.

    > Note: You should only really be connecting to the VM with the intention of creating a new snapshot or running tests with the current set-up. Avoid making any software changes not intended for a new snapshot.
    > Note: The person who creates the VM will have to manually added the SSH keys of any other potential user.

3. For `ShapePipe`, install the following tools:

    ```bash
    sudo apt update
    sudo apt install git
    sudo apt install make
    sudo apt install autoconf
    sudo apt install libtool
    ```

4. For `ShapePipe`, install miniconda:

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    bash
    ```

5. Install `ShapePipe`:
    * Add VM SSH key to the `ShapePipe` GitHub repository

    ```bash
    ssh-keygen -t rsa -b 4096 -C “EMAIL_ADDRESS”
    cat .ssh/id_rsa.pub
    ```
    On github, click on your thumbnail in the upper right corner, choose the Settings menu entry, and click on "SSH Keys" on the left. Create a new ssh key by copying the contents of `ssh/id_rsa.pub` to the corresponding field, and provide a name.

    * Clone repository

    ```bash
    git clone git@github.com:CosmoStat/shapepipe.git
    ```

    * Run the install script

    ```bash
    cd ShapePipe
    ./install_shapepipe --vos
    ```

    * Activate the `SshapePipe` environment

    ```bash
    conda activate shapepipe
    ```

6. Generate certificate to access VOSPACE:

    ```bash
    cadc-get-cert -u USERNAME
    ```
    When asked, enter your CADC password.
    A CADC certificate is also needed to transfer data to/from the VOSPACE.

7. Test that pipeline is working and can access multiple VM CPUs:

    ```bash
    ./shapepipe_run.py -c example/config.ini
    ```

    If the test log returns the expected results the VM should be ready for a snapshot.

8. Create a snapshot of the VM status:

    On [OpenStack](https://arbutus-canfar.cloud.computecanada.ca/) under "Instances" click the "Create Snapshot" button for the corresponding VM. Be sure to follow the snapshot naming scheme defined for the VM above.

It will take a few to a few ten minutes until a snapshot is ready to be used.

9. Update VM and create a new snapshot:

    The VM set-up only needs to be done once, afterwards the VM can simply be modified for new snapshots. *e.g.* pull the latest changes to the pipeline repository and install any new dependencies then repeat step 9.

## Batch System

The batch system is a server where jobs can be submitted to the CANFAR cluster using a previously defined snapshot.

> Note: You will have to request access from [Sebastian Fabbro](sebfabbro@gmail.com) to the batch system before you can connect.

1. SSH to batch system:

    You can connect to the batch system as follows:

    ```bash
    ssh USERNAME@batch.canfar.net
    ```

    You will be connected to a personal user space.

2. Source OpenStack environment variables:

    ```bash
    source lensing-openrc.sh
    ```
    When asked, enter your CADC password.
    This is a necessary step before submitting jobs.

3. Create a bash script, for example:

    The bash script defines the command lines to be run on the snapshot. The following example script demonstrates how to:
     - activate the ShapePipe environment,
     - create an output directory,
     - copy a configuration file to the snapshot from the VOSPACE,
     - run ShapePipe,
     - and copy the output back to the VOSPACE

    ```bash
    #!/bin/bash
    export VM_HOME=/home/ubuntu
    export SP_ROOT=$VM_HOME/ShapePipe
    source $VM_HOME/miniconda3/bin/activate shapepipe
    mkdir output
    vcp vos:cfis/cosmostat/sfarrens/config.ini .
    shapepipe_run -c config.ini
    vcp --certfile=$VM_HOME/.ssl/cadcproxy.pem output vos:cfis/cosmostat/sfarrens
    ```

    > Note: The default path for a snapshot is not the `/home/ubuntu` directory, hence the definition of the `VM_HOME` environment variable.

4. Create a job file, for example:

    The job file defines the script to be run (*i.e.* the bash script previously defined), the corresponding outputs and the computational requirements for the job.

    ```txt
    executable     = shapepipe.bash

    output         = shapepipe.out
    error          = shapepipe.err
    log            = shapepipe.log

    # Make sure the requested resources do not exceed what was
    # specified for the VM
    request_cpus   = 1
    request_memory = 8G
    request_disk   = 10G

    queue
    ```

5. Submit a job:

    Jobs are submitted using the `canfar_submit` command followed by the previously defined job file, the name of the desired snapshot and the *flavour* of the corresponding VM.

    ```bash
    > canfar_submit JOB_FILE SNAP_SHOT FLAVOUR
    ```
 Make sure that the snapshot is ready to be used, see below to check its status.

6. Check queue:

    ```bash
    condor_status -submitter
    ```
  This command tells you running, idle, and held jobs for you and other users.

  Information for your own jobs only:
  ```bash
  condor_q [-nobatch]
  ```
  From there you can get the job ID, which lets you examine your job more closely:
  ```bash
  condor_q -better-analyse <ID>
  ```

  You can do an `ssh` to the VM that is (or will be) running your job for checking:
  ```bash
  condor_ssh_to_job -auto-retry <ID>
  ```
  For multi-job submissions, the JOB_IDS has subnumbers, e.g. `1883.0-9`.
  You can `ssh` to each of those VMs, with e.g.
  ```bash
  condor_ssh_to_job -auto-retry 1883.6
  ```

## Troubleshooting

If the above `condor` commands do not help, try:
```bash
cloud_status -m
```
to check the status of all VMs.

Sometime an snap shot image is not (yet) active and shared, since its creation can take a lot of time. Check the status with:
```bash
openstack image show -c visibility -c status <SnapShotName>
```
When `status = active`, the job can be started. The field `visibiltiy` has the value `private` before first use, which afterwards changes to `shared`.

In general, a job should be started within 5 - 10 minutes. This time will increase if the queue is full. If the job is launched before the snap shot status is `active`, it might be stuck in the queue for a long time (for ever?).

Contact Seb on the canfar slack channel, he usually replies quickly, sometimes there are issues that only he can fix.


## Interactive Mode

The VM can be used interactively for testing, e.g. to make sure an
executable bash script is set up correctly and runs without errors.

1. Copy executable to VM

   Copy via `scp` the bash script to the virtual machine. Do this from your
   machine whose ssh keys are stored, from where you can ssh to the batch
   system and the VM. Copy the script to `/mnt/scratch`, which is a file
   system that is not stored and updated with the VM, so the VM status is
   not changed:
   ```bash
   scp USERNAME@batch.canfar.net:script.bash .
   scp script.bash ubuntu@IP_ADDRESS:/mnt/scratch
   ```

2.  SSH to VM:

    As above, run the following to connect to a given VM:

    ```bash
    ssh ubuntu@IP_ADDRESS
    ```

3. Run the executable

   Change to the directory outside the VM and run the script:
   ```bash
   cd /mnt/scratch
   bash script.bash
   ```
   In case the script needs to write to VM-directories in $VM_ROOT,
   you can create symbolic links to /mnt/scratch.
   
## Running a CFIS job

### Job preparation

1. Make sure the virtual machine is active, the correct version/branch of `ShapePipe` is installed, and the `cadc` certificate is valid.

2. Make sure the desired configuration files are uploaded to `vos`. The corresponding files are in the directory `example/cfis`, and can be copied to `vos` with:
```bash
cd /path/to/shapepipe/example
vcp cfis vos:cfis/cosmostat/kilbinger
```

3. Make sure the `results` directory on `vos` exists:
```bash
vls vos:cfis/cosmostat/kilbinger
```
Results (log and catalogue FITS files) will be uploaded there for each tiles.

It is recommended that this directory is empty, and does not have files from previous runs. The simplest way to clean up is
```bash
vrmdir vos:cfis/cosmostat/kilbinger/results
vmkdir vos:cfis/cosmostat/kilbinger/results
```

### Set up job file and submit job

The two following alternative ways exist to set up and submit a job. For both it is recommended to `cd` into a new subdirectory on the batch system.

1. Enter tile IDs by hand

   CFIS tiles can be run with the bash script [canfar_sp.bash](https://github.com/CosmoStat/shapepipe/blob/master/scripts/sh/canfar_sp.bash). On the VM via a    job script, for example:
   ```bash
   executable     = ./canfar_sp.bash

   output         = log_sp_tile_$(arguments).out
   error          = log_sp_tile_$(arguments).err
   log            = log_sp_tile_$(arguments).log

   request_cpus   = 8
   request_memory = 19G
   request_disk   = 100G


   queue arguments from (
   277.282
   )
   ```
   To launch a job with more tiles, simply add the corresponding tile IDs to the `queue arguments from (` list.
   See (#batch-system) point 5. how to submit a job.

   In interactive mode, type
   ```bash
   bash canfar_sp.bash 277.282
   ```
   To run the script on more than one tile, add the IDs as command line arguments.

2. Using tile IDs from file

   Use the bash script `canfar_submit_selection.sh` to automatically read tile IDs from an ASCII file, create the job file,
   and submit the job.
   ```bash
   bash ~/shapepipe/scripts/sh/canfar_submit_selection.sh ~/shapepipe/aux/CFIS/tiles_202007/tiles_all_order.txt ShapePipe2-mk-20200701b [-n]
   ```
   Add the option `-n` (dry run) to test without submitting the job.

### Analysis

For a summary of the status of submitted jobs from the current directory, type
```bash
python3 ~/shapepipe/scripts/python/canfar_run_analyse.py
```
Results of jobs are uploaded to `vos`. These can be complete results from jobs that finished with success, or partial
results (e.g. log files) from jobs that were stopped due to errors. Download and post-processing analysis of those results should
to be performed on a different machine, e.g. [candide](https://github.com/CosmoStat/shapepipe/blob/master/docs/wiki/candide.md).
