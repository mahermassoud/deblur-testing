from unittest import TestCase, main
import unittest
from methods import *
from qiime2 import Artifact
from qiime2 import MetadataColumn, CategoricalMetadataColumn
from qiime2.plugins.demux.methods import emp_single
import biom
import numpy as np
from numpy.testing import assert_array_almost_equal
from pandas.util.testing import assert_frame_equal
import sys

NUM_CORES = 4
dir_path = os.path.dirname(os.path.realpath(__file__))

@unittest.skip("Skipping test import because it takes long time")
class TestImport(TestCase):

    def setUp(self):
        print("Importing metadata for expected")
        self.exp_barcode_metadata = \
            Metadata.load(dir_path + "/data/mock-3/sample-metadata.tsv")

        self.exp_demux = Artifact.load(dir_path + "/data/mock-3/exp_demux.qza")
        self.exp_out = [self.exp_demux, self.exp_barcode_metadata]
        self.working_dir_fp = dir_path + "/data/mock-3"

    def test_import_dataset(self):
        obs = import_dataset(self.working_dir_fp, "BarcodeSequence")
        for i, entry in enumerate(obs):
            try:
                self.assertEqual(entry, self.exp_out[i])
            except AssertionError as e:
                # TODO fix this
                if("uuid" in str(e)): # uuid is always unique so we expect this
                    continue
                else:
                    self.fail(str(e))

@unittest.skip("Skipping test deblur because it takes long time")
class TestDeblur(TestCase):

    def setUp(self):
        self.exp_demux = Artifact.load(dir_path + "/data/mock-3/exp_demux.qza")
        self.exp_deblurred = Artifact.load(dir_path + "/data/mock-3/deblurred_150nt.qza")
        self.exp_deblur_biom = self.exp_deblurred.view(biom.Table)
        self.exp_deblurred_pt = Artifact.load(dir_path + "/data/mock-3/deblurred_100nt_pt.qza")
        self.num_parallel = NUM_CORES

    def test_establish_dataset(self):
        obs = do_deblur(self.exp_demux, 150)
        obs_biom = obs.view(biom.Table)
        obs_biom = obs_biom.sort_order(self.exp_deblur_biom.ids(axis='observation'), axis='observation')
        obs_biom = obs_biom.sort_order(self.exp_deblur_biom.ids(axis='sample'), axis='sample')

        self.assertEqual(self.exp_deblur_biom, obs_biom)

    def test_establish_dataset_ncores(self):
        obs = do_deblur(self.exp_demux, 150, num_cores=self.num_parallel)
        obs_biom = obs.view(biom.Table)
        obs_biom = obs_biom.sort_order(self.exp_deblur_biom.ids(axis='observation'), axis='observation')
        obs_biom = obs_biom.sort_order(self.exp_deblur_biom.ids(axis='sample'), axis='sample')

        self.assertEqual(self.exp_deblur_biom, obs_biom)

class TestPairwiseDist(TestCase):
    def setUp(self):
        self.jaccard_in = biom.Table(np.array([[1,1,0],[0,1,1]]),
                                     ['A', 'B'], ['S1', 'S2','S3'])
        self.jaccard_exp = np.array([[ 0. ,  0.5,  1. ],
                                     [ 0.5,  0. ,  0.5],
                                     [ 1. ,  0.5,  0. ]])
        self.bc_in = biom.Table(np.array([[1,1,0],[0,2,3]]),
                                ['A', 'B'], ['S1', 'S2','S3'])
        self.bc_exp = np.array([[ 0.,  0.5,  1.],
                                [ 0.5,  0., 0.33333333],
                                [ 1.,  0.33333333,  0.]])

        self.pre = biom.Table(np.array([[1,1,0],[0,1,1]]),
                        ['A', 'B'], ['S1', 'S2','S3'])
        self.post = biom.Table(np.array([[10,50,0],[0,0,30]]),
                                ['A', 'B'], ['S1', 'S2','S3'])
        self.tl = [100]
        self.exp_df = pd.DataFrame(columns = ["trim_length", "dist_type", "r",
                                              "nsamples"],
                                   data= np.array([
                                                  [100,"jaccard", 0.49999999, 3],
                                                  [100,"braycurtis",0.5, 3]]))
        self.exp_df["r"] = pd.to_numeric(self.exp_df["r"])
        self.exp_df["trim_length"] = pd.to_numeric(self.exp_df["trim_length"])
        self.exp_df["nsamples"] = pd.to_numeric(self.exp_df["nsamples"])
        self.exp_df["r_sq"] = self.exp_df["r"] ** 2

    def test_jaccard(self):
        obs = get_pairwise_dist_mat(self.jaccard_in,"jaccard")
        assert_array_almost_equal(self.jaccard_exp, obs.data)

    def test_bc(self):
        obs = get_pairwise_dist_mat(self.bc_in,"braycurtis")
        assert_array_almost_equal(self.bc_exp, obs.data)

    def test_get_pairwise_diversity(self):
        obs = get_pairwise_diversity_data([self.pre], [self.post], self.tl)
        obs = obs.drop(["pval"], axis=1)
        assert_frame_equal(self.exp_df, obs, check_less_precise=True,
                           check_dtype=False)


class TestPrePostDist(TestCase):
    def setUp(self):
            self.pre_100 = biom.Table(np.array([[0,0],[1,1],[1,1],[1,1]]),
                                  ['A', 'B', 'C','D'], ['S1', 'S2'])
            self.post_100 = biom.Table(np.array([[1,1],[0,1],[1,1],[1,2]]),
                                  ['A', 'B', 'C','D'], ['S1', 'S2'])
            self.pre_100_sample = biom.Table(np.array([[0,0],[1,1],[1,1],[1,1]]).transpose(),
                                             ['X','Y'],["S1","S2","S3","S4"])
            self.pre_100_sample_extra = biom.Table(np.array([[0,0,6],[1,1,6],[1,1,6],[1,1,6]]).transpose(),
                                             ['X','Y',"Z"],["S1","S2","S3","S4"])
            self.post_100_sample = biom.Table(np.array([[1,1],[0,1],[1,1],[1,2]]).transpose(),
                                             ['X','Y'],["S1","S2","S3","S4"])

            self.post1_100_oo = biom.Table(np.array([[1,1],[2,1],[1,1],[1,0],[99,99]]),
                                       ['A', 'D', 'C','B',"E"], ['S2', 'S1'])
            self.pre_50 = self.pre_100
            self.post_50 = self.post_100

            self.exp_pre_100_ov = biom.Table(np.array([[0,0],[1,1],[1,1],[1,1]]),
                                       ['A', 'D', 'C','B'], ['S2', 'S1'])
            self.exp_post1_100_oo_ov = biom.Table(np.array([[1,1],[2,1],[1,1],[1,0]]),
                                           ['A', 'D', 'C','B'], ['S2', 'S1'])

            self.exp_100 = pd.DataFrame(columns=["seq","dist_type","dist"],
                                    data= np.array([
                                        ["A","jaccard", 1.0],
                                        ["A","braycurtis",1.0],
                                        ["B","jaccard",0.5],
                                        ["B","braycurtis",1/3],
                                        ["C","jaccard",0.0],
                                        ["C","braycurtis",0.0],
                                        ["D","jaccard",0.0],
                                        ["D","braycurtis",0.2],
                                    ]))
            self.exp_100_sample = pd.DataFrame(columns=["sample","dist_type","dist"],
                                        data= np.array([
                                            ["S1","jaccard", 1.0],
                                            ["S1","braycurtis",1.0],
                                            ["S2","jaccard",0.5],
                                            ["S2","braycurtis",1/3],
                                            ["S3","jaccard",0.0],
                                            ["S3","braycurtis",0.0],
                                            ["S4","jaccard",0.0],
                                            ["S4","braycurtis",0.2],
                                        ]))
            self.exp_both = pd.DataFrame(columns=["seq","dist_type","dist",
                                                  "length"],
                               data= np.array([
                                   ["A", "jaccard",    1.0, 100],
                                   ["A", "braycurtis", 1.0, 100],
                                   ["B", "jaccard",    0.5, 100],
                                   ["B", "braycurtis", 1/3, 100],
                                   ["C", "jaccard",    0.0, 100],
                                   ["C", "braycurtis", 0.0, 100],
                                   ["D", "jaccard",    0.0, 100],
                                   ["D", "braycurtis", 0.2, 100],
                                   # next
                                   ["A", "jaccard",    1.0, 50],
                                   ["A", "braycurtis", 1.0, 50],
                                   ["B", "jaccard",    0.5, 50],
                                   ["B", "braycurtis", 1/3, 50],
                                   ["C", "jaccard",    0.0, 50],
                                   ["C", "braycurtis", 0.0, 50],
                                   ["D", "jaccard",    0.0, 50],
                                   ["D", "braycurtis", 0.2, 50]
                               ]))
            self.exp_both_sample = pd.DataFrame(columns=["sample","dist_type","dist",
                                                  "length"],
                                         data= np.array([
                                             ["S1", "jaccard",    0.5,  100],
                                             ["S1", "braycurtis", 1/3,  100],
                                             ["S2", "jaccard",    0.25, 100],
                                             ["S2", "braycurtis", 0.25, 100],
                                             ["S1", "jaccard",    0.5,  50],
                                             ["S1", "braycurtis", 1/3,  50],
                                             ["S2", "jaccard",    0.25, 50],
                                             ["S2", "braycurtis", 0.25, 50],
                                         ]))

            self.exp_100["dist"] = pd.to_numeric(self.exp_100["dist"])
            self.exp_100_sample["dist"] = pd.to_numeric(self.exp_100_sample["dist"])
            self.exp_both["dist"] = pd.to_numeric(self.exp_both["dist"])
            self.exp_both["length"] = pd.to_numeric(self.exp_both["length"])
            self.exp_both_sample["dist"] = pd.to_numeric(self.exp_both_sample["dist"])
            self.exp_both_sample["length"] = pd.to_numeric(self.exp_both_sample["length"])

    def test_get_distance_distribution(self):
        obs = get_distance_distribution(self.pre_100,self.post_100)
        assert_frame_equal(self.exp_100, obs)

    def test_get_distance_distribution_samplewise(self):
        obs = get_distance_distribution(self.pre_100_sample, self.post_100_sample, by_sample=True)
        assert_frame_equal(self.exp_100_sample, obs)

    def test_get_pre_post_distances(self):
        obs, obs_sample, po, pto = get_pre_post_distance_data([self.pre_100, self.pre_50],
                                         [self.post_100, self.post_50],
                                         [100, 50])
        obs.index = list(range(obs.shape[0]))
        obs_sample.index = list(range(obs_sample.shape[0]))
        assert_frame_equal(self.exp_both, obs)
        assert_frame_equal(self.exp_both_sample, obs_sample)

    def test_get_pre_post_distances_fail(self):
        with self.assertRaises(ValueError):
              get_pre_post_distance_data([self.pre_100, self.pre_50],
                                         [self.post_100, self.post_50], [100, 50, 1])

    def test_get_overlap_tables_oo(self):
        obs_pre_ov, obs_post_ov = get_overlap_tables(self.pre_100,
                                                     self.post1_100_oo)

        self.assertEqual(self.exp_pre_100_ov, obs_pre_ov,
                         msg=self.exp_pre_100_ov.descriptive_equality(obs_pre_ov))
        self.assertEqual(self.exp_post1_100_oo_ov, obs_post_ov,
                         msg=self.exp_post1_100_oo_ov.descriptive_equality(obs_post_ov))

class TestPostTrim(TestCase):
    def setUp(self):
        self.t = biom.Table(np.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[3,1,2,7]]),
                            ['AATT', 'AATG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.exp_pt = biom.Table(np.array([[7,7,10,17],[8,9,10,11]]),
                            ['AAT', 'ATG'], ['S1', 'S2', 'S3', 'S4'])
        self.exp_md = ["AATT", "AATG", "AATC"]
        self.exp_clps = pd.DataFrame({"seq": ["AAT", "ATG"], "length": [3, 3],
                                      "num_collapses": [3, 1]})

    def test_post_trim(self):
        obs = post_trim(self.t, 3)

        # Sort them
        obs = obs.sort_order(self.exp_pt.ids(axis="observation"),
                             axis="observation")
        self.assertEqual(str(self.exp_pt), str(obs))

    def test_post_trim_parallel_even(self):
        obs = post_trim(self.t, 3, 2)

        # Sort them
        obs = obs.sort_order(self.exp_pt.ids(axis="observation"),
                             axis="observation")
        self.assertEqual(str(self.exp_pt), str(obs))

    def test_post_trim_parallel_1(self):
        obs = post_trim(self.t, 3, 1)

        # Sort them
        obs = obs.sort_order(self.exp_pt.ids(axis="observation"),
                             axis="observation")
        self.assertEqual(str(self.exp_pt), str(obs))

    def test_post_trim_parallel_odd(self):
        obs = post_trim(self.t, 3, 3)

        # Sort them
        obs = obs.sort_order(self.exp_pt.ids(axis="observation"),
                             axis="observation")
        self.assertEqual(str(self.exp_pt), str(obs))

    def test_pt_metadata(self):
        obs = post_trim(self.t, 3)
        obs_md = obs.metadata(id="AAT", axis="observation")["collapsed_ids"]

        self.assertCountEqual(self.exp_md, obs_md)

    def test_get_collapse_count(self):
        obs = get_collapse_counts([post_trim(self.t, 3)])
        obs = obs.sort_values(by=["seq"])
        obs.index = [0, 1]
        print(obs)
        exp_srt = self.exp_clps.sort_values(by=["seq"])
        exp_srt.index = [0, 1]
        print(exp_srt)

        assert_frame_equal(exp_srt, obs)

class TestShortSeq(TestCase):
    def setUp(self):
        barcode_map = pd.Series(['GTCA', 'TCAG', 'GGGG'],
                   index=['sample1', 'sample2', 'sample3'], name="aname")
        barcode_map.index.name = "sample_name"
        barcode_map = CategoricalMetadataColumn(barcode_map)

        seqs_fp = dir_path + "/data/small/"

        seqs = Artifact.import_data("EMPSingleEndSequences",
                            seqs_fp)

        self.demuxed, = emp_single(seqs, barcode_map)
        self.exp = 1

    def test_get_shortest_seq(self):
        obs = get_shortest_seq(self.demuxed)
        self.assertEqual(self.exp, obs)

class TestGetOverlap(TestCase):
    def setUp(self):
        self.some1 = biom.Table(np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]]),
                            ['AATT', 'AATG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.some2 = biom.Table(np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]),
                             ['AATG', 'AATT', 'ATGC','ACTC'], ['S1', 'S2', 'S3', 'S4'])
        self.some1_exp = biom.Table(np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3]]),
                                ['AATG', 'AATT', 'ATGC'], ['S1', 'S2', 'S3', 'S4'])
        self.some2_exp = biom.Table(np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0]]),
                                ['AATG', 'AATT', 'ATGC'], ['S1', 'S2', 'S3', 'S4'])
        self.none = biom.Table(np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]]),
                                ['AA', 'GG', 'CC','DD'], ['S1', 'S2', 'S3', 'S4'])

    def test_some_overlap(self):
        obs1, obs2 = get_overlap_tables(self.some1, self.some2)
        self.assertEqual(self.some1_exp, obs1)
        self.assertEqual(self.some2_exp, obs2)

    def test_no_overlap(self):
        with self.assertRaises(ValueError):
            get_overlap_tables(self.none, self.some1)

class TestPartitionTable(TestCase):
    def setUp(self):
        self.tbl = biom.Table(np.array([[0,1,2,3,4],[5,6,7,8,9]]),["A","B"],
                              ["S1","S2","S3","S4","S5"])
        self.p1 = biom.Table(np.array([[0,1],[5,6]]),["A","B"],
                              ["S1","S2"])
        self.p2 = biom.Table(np.array([[2,3],[7,8]]),["A","B"],
                              ["S3","S4"])
        self.p3 = biom.Table(np.array([[4],[9]]),["A","B"],
                              ["S5"])
        self.exp = [self.p1, self.p2, self.p3]

    def test_partition_drop(self):
        obs = partition_table(self.tbl, 3)
        self.assertEqual(self.exp, obs)

class TestGetCountData(TestCase):
    def setUp(self):
        self.pre1 = biom.Table(np.array([[0,10,20,30],[40,20,10,30],[10,15,12,11],[7,8,8,9]]),
                                ['AATT', 'AATG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.pre2 = biom.Table(np.array([[0,11,21,31],[41,21,11,31],[11,16,13,12],[8,9,9,10]]),
                               ['AATT', 'AATG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.post1 = biom.Table(np.array([[0,1,2,3],[4,2,1,3],[1,1,1,1],[8,9,9,10]]),
                               ['AATT', 'ATTG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.post2 = biom.Table(np.array([[0,1,2,3],[4,2,1,3],[1,1,1,1],[8,9,9,10]]),
                                ['ACCT', 'TATG', 'TCCC','ACCC'], ['S1', 'S2', 'S3', 'S4'])

        self.pre_ov1, self.post_ov1 = get_overlap_tables(self.pre1, self.post1)
        self.pre_ov2, self.post_ov2 = get_overlap_tables(self.pre2, self.post2)

        self.pre_bioms = [self.pre1, self.pre2]
        self.post_bioms = [self.post1, self.post2]
        self.pre_ov = [self.pre_ov1, self.pre_ov2]
        self.post_ov = [self.post_ov1, self.post_ov2]
        self.lengths = [1,2]

        self.exp_cdata = pd.DataFrame(columns = ["trim_length",
                                                "sOTU_overlap_count",
                                                "sOTU_unique_pre",
                                                "sOTU_unique_post",
                                                "diff_otu"],
                                      data = np.array([[1,3,1,1,0],
                                                       [2,0,4,4,0]]))
        self.exp_rps = pd.DataFrame(columns = ['trim_length','S1', 'S2', 'S3', 'S4'],
                                    data = np.array([[1,44,40,37,63],
                                                     [2,47,44,41,67]]))

    def test_get_count_data(self):
        obs_cdata, obs_rps = get_count_data(self.pre_bioms, self.pre_ov,
                                            self.post_bioms, self.post_ov,
                                            self.lengths)

        assert_frame_equal(self.exp_cdata, obs_cdata)
        assert_frame_equal(self.exp_rps, obs_rps, check_dtype=False)


if __name__ == '__main__':
    main()
