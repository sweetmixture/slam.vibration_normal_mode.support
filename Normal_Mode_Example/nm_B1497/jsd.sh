#!/bin/bash -l
#$ -cwd
#$ -m e
#$ -M uccawkj@master
#$ -q all.q
#$ -N PbF2_f5
#$ -pe smp 1

module load intel/2018.3
bash /home/uccawkj/src_tool/SLAM/slam_vib_package_v_2.2/slam_vib.auto_geo_gen.sh 0.015 > res.out
