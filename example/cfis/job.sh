#!/usr/bin/bash

if [ "$1" == "-h" ]
then
	echo "Usage [on candide]: qsub -d . -N <name> job.sh -F <config_file>"
	exit 0
fi

#PBS -S /usr/bin/bash

#PBS -N ShapePipe

#PBS -j oe

#PBS -l nodes=1:ppn=4,walltime=10:00:00

#PBS -d /home/mkilbing/area_W3

echo -n "pwd = "
pwd

cmd="$HOME/ShapePipe/shapepipe_run.py -c $1"
echo "Running command $cmd"
$cmd
ex=$?
exit $ex


