#!/bin/bash

#PBS -l nodes=1:ppn=5

#PBS -l walltime=05:00:00

#PBS -l mem=37gb

#PBS -N deblur_test_ss_emp

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe
#PBS -t 1-5


fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
fp_100nt="/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
fp_90nt="/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom"
i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_arr/ss$PBS_ARRAYID"

source activate q2env
cd $i_dir
printf "identifier\truntime\n" > time_elapsed.tsv

time sbiom -i emp_deblur_150bp_ss_10.biom \
           -i emp_deblur_100bp_ss_10.biom \
           -i emp_deblur_90bp_ss_10.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_10_arr-$PBS_ARRAY_ID

time sbiom -i emp_deblur_150bp_ss_32.biom \
           -i emp_deblur_100bp_ss_32.biom \
           -i emp_deblur_90bp_ss_32.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_32_arr-$PBS_ARRAY_ID

time sbiom -i emp_deblur_150bp_ss_55.biom \
           -i emp_deblur_100bp_ss_55.biom \
           -i emp_deblur_90bp_ss_55.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_55_arr-$PBS_ARRAY_ID

time sbiom -i emp_deblur_150bp_ss_77.biom \
           -i emp_deblur_100bp_ss_77.biom \
           -i emp_deblur_90bp_ss_77.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_77_arr-$PBS_ARRAY_ID

time sbiom -i emp_deblur_150bp_ss_100.biom \
           -i emp_deblur_100bp_ss_100.biom \
           -i emp_deblur_90bp_ss_100.biom  \
           -o . \
           -to time_elapsed.tsv \
           -toa ss_100_arr-$PBS_ARRAY_ID
