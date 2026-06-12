#!/bin/ksh
#$ -cwd -S /bin/sh -j y
#$ -q parallel -pe mpi 12
#$ -N out
#$ -l h_vmem=1024M
#
limit cputime unlimited
#
. /usr/local/codine/codine_settings.sh
. /opt/modules/modules/init/sh
module load MIPSpro mpt
#
