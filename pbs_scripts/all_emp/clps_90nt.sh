#!/bin/bash

#PBS -l nodes=1:ppn=2

#PBS -l walltime=50:00:00

#PBS -l mem=100gb

#PBS -N clps_emp_90nt

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

i_path="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/90nt/emp_deblur_150bp_release1_pt_9090.biom"
o_path="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out/90nt/collapse.csv"

source activate q2env

time clps_count -i $i_path \
                -o $o_path
