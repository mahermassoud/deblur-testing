#!/bin/bash

#PBS -l nodes=1:ppn=3

#PBS -l walltime=05:00:00

#PBS -l mem=24gb

#PBS -N deblur_test_ss1_emp

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
fp_100nt="/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
fp_90nt="/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom"
i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss1"

cd $i_dir
printf "identifier\truntime\n" > ss1_10/time_elapsed.tsv
source activate q2env

time sbiom -i emp_deblur_150bp_ss_10.biom \
           -i emp_deblur_100bp_ss_10.biom \
           -i emp_deblur_90bp_ss_10.biom  \
           -o ss1_10/ \
           -to ss1_10/time_elapsed.tsv \
           -toa ss_10_arr_1
