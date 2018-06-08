#!/bin/bash

#PBS -l nodes=1:ppn=5

#PBS -l walltime=05:00:00

#PBS -l mem=20gb

#PBS -N emp_pt_150_to_100

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe
#PBS -t 1-2


fp_150nt="/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
fp_100nt="/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
fp_90nt="/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom"
i_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_all_out"

source activate q2env
mkdir -p $i_diir
cd $i_dir

time rpost_trim -i $fp_150nt \
                -o . \
                -tl 100 \
                -on emp_deblur_150_pt_