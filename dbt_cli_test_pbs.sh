#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=1:ppn=3

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_unittest_job1

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

source activate q2env
rall -i /panfs/panfs1.ucsd.edu/panscratch/mamaher/deblur-testing/test/data/mock-3/emp-single-end-sequences -m /panfs/panfs1.ucsd.edu/panscratch/mamaher/deblur-testing/test/data/mock-3/sample-metadata.tsv -mbc BarcodeSequence -l 150 -nc 3 -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/rall_out/


# comment PBS -o /home/mamaher/dbt_job_out
# comment PBS -e /home/mamaher/dbt_job_err
# comment PBS -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_out
# comment PBS -e /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_err
