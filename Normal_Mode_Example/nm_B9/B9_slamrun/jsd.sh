#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N PbF2_f5
#$ -pe smp 2

module load intel/2018.3
mpirun -np $NSLOTS /home/uccawkj/bin/slam.240122.mpi.x slam.xyz 1000 > slam.out 
