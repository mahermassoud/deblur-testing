#!/bin/bash

#PBS -l nodes=1:ppn=10

#PBS -l walltime=00:30:00

#PBS -l mem=10gb

#PBS -N check pr_parallel_pt_ppn10

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe
#PBS -t 1-3

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/check_parallel"
o_name="rand_1000x1000_arr${PBS_ARRRAY_ID}_ppn_10"

source activate q2env
cd $i_dir

time rpost_trim -i  rand_1000x1000.qza \
                -o . \
                -tl 100 \
                -on $o_name \
                -to times.tsv \
                -toa $o_name
