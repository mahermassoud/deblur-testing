#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=3:ppn=1

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N dbt_unittest_job1

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

cd $PBS_O_WORKDIR
source activate q2env
echo "---------------------------------------------------------------------"
echo "-------------------Running test_methods.py---------------------------"
echo "---------------------------------------------------------------------"
touch /panfs/panfs1.ucsd.edu/panscratch/mamaher/a_foo
#python /panfs/panfs1.ucsd.edu/panscratch/mamaher/deblur-testing/test/test_methods.py

echo "---------------------------------------------------------------------"
echo "-------------------Running test_script.py---------------------------"
echo "---------------------------------------------------------------------"
#python /panfs/panfs1.ucsd.edu/panscratch/mamaher/deblur-testing/test/test_methods.py

# comment PBS -o /home/mamaher/dbt_job_out
# comment PBS -e /home/mamaher/dbt_job_err
# comment PBS -o /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_out
# comment PBS -e /panfs/panfs1.ucsd.edu/panscratch/mamaher/dbt_job_err
