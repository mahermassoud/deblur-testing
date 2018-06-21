#!/bin/bash

#PBS -l nodes=1:ppn=2

#PBS -l walltime=150:00:00

#PBS -l mem=100gb

#PBS -N clps_emp_100nt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
i_path="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/100nt/emp_deblur_150bp_release1_pt_100nt100.biom"
o_path="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/100nt/collapse.csv"

source activate q2env

time clps_count -i $i_path \
                -o $o_path
