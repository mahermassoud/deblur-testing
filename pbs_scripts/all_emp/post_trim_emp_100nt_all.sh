#!/bin/bash

#PBS -l nodes=1:ppn=8

#PBS -l walltime=150:00:00

#PBS -l mem=220gb

#PBS -N post_trim_emp_100nt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/100nt"
o_name="emp_deblur_150bp_release1_pt_100nt"

source activate q2env
mkdir -p $o_dir
cd $o_dir

time rpost_trim -ib $fp_150nt \
                -o $o_dir \
                -tl 100 \
                -on $o_name \
                -to elapsed.tsv \
                -toa $o_name \
                -pc 8 \
                -sb
