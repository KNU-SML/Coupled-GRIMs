#!/bin/sh
#PBS -l walltime=0:30:00
#PBS -l rmsnodes=1:12
#PBS -j oe
#
set -ax
for envvar in `env | grep PBS | cut -d'=' -f1`
do
        export $envvar
done
export RSM_NODES
export RSM_PROCS
