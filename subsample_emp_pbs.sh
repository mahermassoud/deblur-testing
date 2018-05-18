#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=3:ppn=1

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_unittest_job1

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

time source activate q2env
echo $(date +"%T")
cd /panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss
time ssbiom -mi /projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom -si /projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom -si barnacle:/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom -s 10 -e 100 -c -o .
