from click.testing import CliRunner
from scripts import post_trims




runner = CliRunner()
pre_fp = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/rand_1000x1000.qza"
o_fp = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt/ppn1"
o_root = "/panfs/panfs1.ucsd.edu/panscratch/mamaher/cmp_pt"
with runner.isolated_filesystem():
    result = runner.invoke(post_trims,    ["-i", pre_fp,
                                             "-o", o_fp,
                                             "-to", o_root + "/prof_elapsed.tsv",
                                             "-toa", "rand500_ppn1_np_no_pc",
                                             "-tl", 100])

    print(result.exc_info)
    print(result.output)
