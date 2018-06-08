from click.testing import CliRunner
from scripts import post_trims

runner = CliRunner()
pre_fp = "/Users/massoudmaher/Documents/Code/deblur-testing/rand_bioms/pre_rand_500x500.biom"
o_fp = "/Users/massoudmaher/Documents/Code/deblur-testing/rand_bioms/sbiom_out"
pre_fp = "/Users/massoudmaher/Documents/Code/deblur-testing/rand_bioms/sbiom_out/pre_rand_500x500.biom.qza"
with runner.isolated_filesystem():
    result = runner.invoke(post_trims,    ["-i", pre_fp,
                                             "-o", o_fp,
                                             "-to", o_fp + "/time.tsv",
                                             "-toa", "rand500",
                                             "-tl", 100,
                                             "-tl", 90])

    print(result.exc_info)
    print(result.output)
