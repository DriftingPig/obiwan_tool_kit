#!/bin/bash -l

#SBATCH -N 1         #Use 2 nodes
#SBATCH -t 00:30:00  #Set 30 minute time limit
#SBATCH -q regular   #Submit to the regular QOS
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes

./ObiwanCorr.sh
