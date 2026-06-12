#!/bin/sh
#PBS -q dque
#PBS -A 
#PBS -N rsm
#PBS -l nodes=1:ppn=12,
#PBS -l walltime=0:30:0
#PBS -o file.out
#PBS -e file.err
#PBS -V
#
for envvar in `env | grep PBS | cut -d'=' -f1`
do
        export $envvar
done
