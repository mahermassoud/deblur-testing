#!/bin/bash

#PBS -l nodes=1:ppn=32

#PBS -l walltime=20:00:00

#PBS -l mem=30gb

#PBS -N post_trim_emp_100nt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

fp_100nt="/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out"

source activate q2env
mkdir -p $o_dir
cd $i_dir

time rpost_trim -ib $fp_100nt \
                -o $o_dir\
                -tl 100 \
                -on $o_name \
                -to elapsed.tsv \
                -toa $o_name \
                -pc 1
