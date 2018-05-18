"""
Quick and dirty python script to find samples that are shared accross
150,100,90nt tables of emp. Should be ran on barnacle
"""
import biom
import time
import numpy as np

emp_150nt_fp = "/projects/emp/03-otus/04-deblur/emp_deblur_150bp.release1.biom"
emp_100nt_fp = "/projects/emp/03-otus/04-deblur/emp_deblur_100bp.release1.biom"
emp_90nt_fp = "/projects/emp/03-otus/04-deblur/emp_deblur_90bp.release1.biom"
OUTPUT_FP = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/emp_shared_ids"

print("Loading 150nt")
start = time.clock()
tbl_150nt = biom.load_table(emp_150nt_fp)
print("Took {}s".format(str(time.clock() - start)))

print("Loading 100nt")
start = time.clock()
tbl_100nt = biom.load_table(emp_100nt_fp)
print("Took {}s".format(str(time.clock() - start)))

print("Loading 90nt")
start = time.clock()
tbl_90nt = biom.load_table(emp_90nt_fp)
print("Took {}s".format(str(time.clock() - start)))

ids = set(tbl_150nt.ids()) & set(tbl_100nt.ids()) & set(tbl_90nt.ids())

with open(OUTPUT_FP, "w") as file:
    for id in ids:
        file.write("{}\n".format(str(id)))