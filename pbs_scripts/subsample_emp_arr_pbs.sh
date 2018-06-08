#!/bin/bash

#PBS -l nodes=1:ppn=5

#PBS -l walltime=05:00:00

#PBS -l mem=37gb

#PBS -N subsample_emp_arr_dbt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe
#PBS -t 1-5


fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
fp_100nt="/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
fp_90nt="/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom"
o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss$PBS_ARRAYID"

source activate q2env
mkdir -p $o_dir
time ssbiom -mi $fp_150nt -oi $fp_100nt -oi $fp_90nt -s 10 -e 100 -c 5 -o $o_dir
