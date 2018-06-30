#!/bin/bash

#PBS -l nodes=1:ppn=2

#PBS -l walltime=200:00:00

#PBS -l mem=220gb

#PBS -N dbt_emp_analysis_all

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/analysis_out"
cdf="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/collapse_all.csv"

source activate q2env
mkdir -p $o_dir
cd $o_dir

time ranalysis -ip ~/dbt_analysis_pathfile \
                -o $o_dir \
                -cfp $cdf \
                -tl 150 \
                -tl 100 \
                -tl 90
