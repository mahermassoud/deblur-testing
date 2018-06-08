#!/bin/bash

#PBS -l nodes=1:ppn=5

#PBS -l walltime=05:00:00

#PBS -l mem=37gb

#PBS -N deblur_test_ss_emp

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss5/"

source activate q2env
cd $i_dir
printf "identifier\truntime\n" > time_elapsed_single.tsv

mkdir -p ss_10
time sbiom -i emp_deblur_150bp_ss_10.biom \
           -i emp_deblur_100bp_ss_10.biom \
           -i emp_deblur_90bp_ss_10.biom  \
           -o ss_10/ \
           -to time_elapsed_single.tsv \
           -toa ss_10

mkdir -p ss_32
time sbiom -i emp_deblur_150bp_ss_32.biom \
           -i emp_deblur_100bp_ss_32.biom \
           -i emp_deblur_90bp_ss_32.biom  \
           -o . \
           -o ss_32/ \
           -to time_elapsed_single.tsv \
           -toa ss_32

mkdir -p ss_55
time sbiom -i emp_deblur_150bp_ss_55.biom \
           -i emp_deblur_100bp_ss_55.biom \
           -i emp_deblur_90bp_ss_55.biom  \
           -o ss_55/ \
           -to time_elapsed_single.tsv \
           -toa ss_55

mkdir -p ss_77
time sbiom -i emp_deblur_150bp_ss_77.biom \
           -i emp_deblur_100bp_ss_77.biom \
           -i emp_deblur_90bp_ss_77.biom  \
           -o ss_77/ \
           -to time_elapsed_single.tsv \
           -toa ss_77

mkdir -p ss_100
time sbiom -i emp_deblur_150bp_ss_100.biom \
           -i emp_deblur_100bp_ss_100.biom \
           -i emp_deblur_90bp_ss_100.biom  \
           -o ss_100/ \
           -to time_elapsed_single.tsv \
           -toa ss_100
