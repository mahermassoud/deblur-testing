#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=1:ppn=3

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_analysis_only

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss1/ss1_10"
o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss1/ss1_10_a_out"

source activate q2env
mkdir -p $o_dir
time ranalysis -i $i_dir -o $o_dir -tl 150 -tl 100 -tl 90
