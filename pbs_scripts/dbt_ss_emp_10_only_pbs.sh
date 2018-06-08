#!/bin/bash

#PBS -l nodes=1:ppn=3

#PBS -l walltime=03:00:00

#PBS -l mem=30gb

#PBS -N deblur_test_ss_emp_arr_10_only

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss2"

source activate q2env
cd $i_dir
mkdir -p ss_10_only_out
printf "identifier\truntime\n" > ss_10_only_out/time_elapsed.tsv

time sbiom -i emp_deblur_150bp_ss_10.biom \
           -i emp_deblur_100bp_ss_10.biom \
           -i emp_deblur_90bp_ss_10.biom  \
           -o ss_10_only_out/ \
           -to ss_10_only_out/time_elapsed.tsv \
           -toa ss_10_only
