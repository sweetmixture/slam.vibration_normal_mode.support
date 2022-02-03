#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N PbF2_f5
#$ -pe smp 2

module load intel/2018.3
mpirun -np $NSLOTS /home/uccawkj/src_tool/bin/slam_v_2.2_opti.mpi.x.bipb out.xyz 4500 > opti.out
