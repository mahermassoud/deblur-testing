#!/bin/bash

#PBS -l nodes=1:ppn=5

#PBS -l walltime=05:00:00

#PBS -l mem=37gb

#PBS -N deblur_test_ss_emp

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss1"

source activate q2env
cd $i_dir
printf "identifier\truntime\n" > time_elapsed.tsv

time sbiom -i emp_deblur_150bp_ss_10.biom \
           -i emp_deblur_100bp_ss_10.biom \
           -i emp_deblur_90bp_ss_10.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_10_arr-1 \
           -tl 150 -tl 100 -tl 90
