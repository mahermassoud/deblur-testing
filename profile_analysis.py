from click.testing import CliRunner
from scripts import post_trims

pre_fp="/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/rand_1000x1000.qza"
o_dir="/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/ppn1"
o_name="rand_1000x1000_ppn_1_parallel_arr_$PBS_ARRAYID"



runner = CliRunner()
pre_fp = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/rand_1000x1000.qza"
o_fp = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/ppn1"
o_root = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt"
with runner.isolated_filesystem():
    result = runner.invoke(post_trims,    ["-i", pre_fp,
                                             "-o", o_fp,
                                             "-to", o_root + "/prof_elapsed.tsv",
                                             "-toa", "rand500_ppn1_np",
                                             "-tl", 100,
                                             "-pc", 1])

    print(result.exc_info)
    print(result.output)
