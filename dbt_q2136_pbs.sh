#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=1:ppn=3

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_unittest_job1

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

time source activate q2env
echo $(date +"%T")
time sbiom -i /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136/reference-hit-100nt_q2136.biom -i /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136/reference-hit-150nt_q2136.biom -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136/


# comment PBS -o /home/mamaher/dbt_job_out
# comment PBS -e /home/mamaher/dbt_job_err
# comment PBS -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_out
# comment PBS -e /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_err
