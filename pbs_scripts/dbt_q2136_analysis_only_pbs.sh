#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=1:ppn=3

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_q2136_analysis_only

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

time source activate q2env
echo $(date +"%T")
touch /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136_2/touch_file
time ranalysis -i /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136_2 -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136_a_out/ -n 2 --trim-incr 33
rplot -i /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136_a_out/ -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/q2136_a_out/

# comment PBS -o /home/mamaher/dbt_job_out
# comment PBS -e /home/mamaher/dbt_job_err
# comment PBS -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_out
# comment PBS -e /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_err
