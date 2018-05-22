#!/bin/bash

# 3 nodes, 1 proccess per node
#PBS -l nodes=3:ppn=1

#PBS -l walltime=02:00:00

#PBS -l mem=24gb

#PBS -N emp_deblur_ss_get_shapes

#PBS -m bea
#PBS -M mahermassoud@gmail.com

#PBS -k oe

time source activate q2env
echo $(date +"%T")
echo "sample observation filename" > /panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_2/shapes
for file in /panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_2/*.biom; do
  echo $file >> /home/mamaher/emp_deblur_ss2/shapes
  biom summarize-table -i $file | awk -F ":" 'FNR == 1 {print $2}; FNR == 2 {print $2};' >> /panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_ss_2/shapes
done
