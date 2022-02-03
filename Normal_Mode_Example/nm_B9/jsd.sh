#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N PbF2_f5
#$ -pe smp 2

module load intel/2018.3
bash /home/uccawkj/bin/slam_vib_package_v_2.2/slam_vib.auto_geo_gen.sh 0.05 > vibres.out
