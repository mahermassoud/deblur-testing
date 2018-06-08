#!/bin/bash

#PBS -l nodes=1:ppn=3

#PBS -l walltime=00:30:00

#PBS -l mem=10gb

#PBS -N check_ppn3_cmp_par_s_pt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt"
o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/ppn3_par_s"
o_name="rand_1000x1000_ppn_3_par_s"

source activate q2env
mkdir -p $o_dir
cd $i_dir
printf "identifier\truntime\n" > $o_dir"/"time_elapsed.tsv

time rpost_trim -i  rand_1000x1000.qza \
                -o $o_dir\
                -tl 100 \
                -on $o_name \
                -to elapsed.tsv \
                -toa $o_name \
                -pc 3
