from unittest import TestCase, main
import unittest
import traceback
from os import listdir
import numpy as np
from pandas.util.testing import assert_frame_equal

from scripts import *
from click.testing import CliRunner


NUM_CORES = 4 # Set for your system
dir_path = os.path.dirname(os.path.realpath(__file__))

"""
Integration tests
"""
class TestDoDemux(TestCase):
    def setUp(self):
        self.seq_fp = dir_path + "/data/small/"
        self.md_fp = dir_path + "/data/small-metadata.tsv"
        self.exp_out = ["demux_small.qza"]

    def test_do_demux(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(do_demux, ["-i", self.seq_fp,
                                              "-o", "demux_small",
                                              "-m", self.md_fp,
                                              "-metadata_bc_col",
                                              "BarcodeSequence"])
            out_files = listdir(os.getcwd())

            self.assertEqual(0, result.exit_code,
                             msg="output:\n{}\nexc.info:\n{}"\
                             .format(result.output, result.exc_info))
            self.assertCountEqual(self.exp_out, out_files)

@unittest.skip("pre trim takes a long time")
class TestPreTrims(TestCase):
    def setUp(self):
        self.demux_fp = dir_path + "/data/mock-3-demo/demux.qza"
        self.exp_out = ["deblurred_pre_150.qza", "deblurred_pre_135.qza",
                        "deblur.log"]


    def test_pre_trims(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(pre_trims, ["-i", self.demux_fp, "-l", 150,
                                               "-o", ".", "-nc", NUM_CORES,
                                               "-n", 2])
            out_files = listdir(os.getcwd())

            self.assertEqual(0, result.exit_code, msg=result.exc_info)
            self.assertCountEqual(self.exp_out, out_files)


class TestPostTrims(TestCase):
    def setUp(self):
        self.db_path = dir_path + "/data/mock-3-demo/deblurred_pre_150.qza"
        self.biom_path = dir_path + "/data/pre_rand_500x500.biom"
        self.exp_out = ["deblurred_pt_100.qza",
                        "deblurred_pt_90.qza", "collapse.csv"]
        self.exp_out_biom = ["rand_500x500_pt_100.qza", "rand_500x500_pt_100.biom",
                             "collapse.csv", "elapsed.tsv"]

    def test_post_trims(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(post_trims, ["-i", self.db_path,
                                                "-o", "."])
            out_files = listdir(os.getcwd())

            self.assertEqual(0, result.exit_code, msg=result.exc_info)
            self.assertCountEqual(self.exp_out, out_files)


    def test_post_trims_biom_in(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(post_trims, ["-ib", self.biom_path,
                                                "-o", ".",
                                                "-tl", 100,
                                                "-on", "rand_500x500_pt_",
                                                "-to", "elapsed.tsv",
                                                "-toa", "rand_500x_500",
                                                "-pc", 2,
                                                "-sb"])
            out_files = listdir(os.getcwd())

            self.assertEqual(0, result.exit_code, msg=result.exc_info)
            self.assertCountEqual(self.exp_out_biom, out_files)

class TestAnalysis(TestCase):
    def setUp(self):
        self.db_path = dir_path + "/data/mock-3-done/"
        self.exp_out = ["pairwise_mantel.csv", "pre_post.csv", "pre_post_sample.csv",
                        "counts.csv", "read_changes.csv"]

    def test_analysis(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(analysis, ["-i", self.db_path, "-o", "."])
            out_files = listdir(os.getcwd())

            self.assertEqual(0, result.exit_code, msg=result.exc_info)
            self.assertCountEqual(self.exp_out, out_files)

class TestSubSampleBiom(TestCase):
    def setUp(self):
        self.tbl_100_1 = biom.Table(np.array([[1,1],[0,0],[0,0],[0,0]]).transpose(),
                                         ['X','Y'],["S1","S2","S3","S4"])
        self.tbl_100_2 = biom.Table(np.array([[2,2,2],[2,2,2],[2,2,2],[2,2,2]]).transpose(),
                                               ['X','Y',"Z"],["S1","S2","S3","S4"])
        self.tbl_100_3 = biom.Table(np.array([[3,3,3],[3,3,3],[3,3,3],[3,3,3]]).transpose(),
                                    ['X','Y',"Z"],["S1","S2","S3","S4"])

    def test_subsample_biom(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open("tbl_100_1.biom", "w") as f1:
                self.tbl_100_1.to_json(generated_by="test_script", direct_io=f1)
            with open("tbl_100_2.biom", "w") as f2:
                self.tbl_100_2.to_json(generated_by="test_script", direct_io=f2)
            with open("tbl_100_3.biom", "w") as f3:
                self.tbl_100_3.to_json(generated_by="test_script", direct_io=f3)

            result = runner.invoke(subsample_biom, ["-mi", "tbl_100_1.biom",
                                                    "-oi", "tbl_100_2.biom",
                                                    "-oi", "tbl_100_3.biom",
                                                    "-s", 1,
                                                    "-e", 3,
                                                    "-c", 3,
                                                    "-o", "."])
            out_files = listdir(os.getcwd())
            print(out_files)
            print("output:")
            print(result.output)
            self.assertEqual(0, result.exit_code,
                             msg="output:\n{}\nexc.info:\n{}\n traceback:{}" \
                             .format(result.output, result.exc_info,
                                     traceback.extract_tb(result.exc_info[2])))


#@unittest.skip("Full integration test takes a long time")
class TestAll(TestCase):
    """Integration test"""
    def setUp(self):
        self.seq_path = dir_path + "/data/mock-3/emp-single-end-sequences"
        self.md_path = dir_path + "/data/mock-3/sample-metadata.tsv"
        self.exp_out = ["demux.qza", "deblurred_pre_150.qza",
                        "deblurred_pre_135.qza", "deblur.log",
                        "deblurred_pt_150.qza", "deblurred_pt_135.qza",
                        "collapse.csv", "pairwise_mantel.csv", "pre_post.csv",
                        "pre_post_sample.csv",
                        "counts.csv", "read_changes.csv", 'collapse.png',
                        'counts.png', 'pairwise_mantel.png', 'pre_post.png',
                        'read_changes.png', 'pre_post_sample.png']

    def test_pre_post(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(pre_post, ["-i", self.seq_path,
                                              "-m", self.md_path,
                                              "-mbc", "BarcodeSequence",
                                              "-l", 150,
                                              "-n", 2,
                                              "-nc", NUM_CORES,
                                              "-o", "."])

            out_files = listdir(os.getcwd())
            self.assertEqual(0, result.exit_code,
                             msg="output:\n{}\nexc.info:\n{}"\
                                 .format(result.output, result.exc_info))
            self.assertCountEqual(self.exp_out, out_files)


if __name__ == '__main__':
    main()
