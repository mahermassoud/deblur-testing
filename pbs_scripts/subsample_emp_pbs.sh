#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=3:ppn=1

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N subsample_emp_deblur_testing

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

source activate q2env
cd home/mamaher/deblur-testing
pip install -e .

echo $(date +"%T")
cd /panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_2
time ssbiom -mi /projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom -oi /projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom -oi /projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom -s 10 -e 100 -c 5 -o .
