#!/bin/bash

#PBS -l nodes=1:ppn=8

#PBS -l walltime=150:00:00

#PBS -l mem=240gb

#PBS -N post_trim_emp_90nt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/90nt"
o_name="emp_deblur_150bp_release1_pt_90"

source activate q2env
mkdir -p $i_dir
cd $i_dir

time rpost_trim -ib $fp_150nt \
                -o $i_dir \
                -tl 90 \
                -on $o_name \
                -to elapsed.tsv \
                -toa $o_name \
                -pc 8 \
                -sb
