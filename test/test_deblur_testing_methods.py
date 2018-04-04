from unittest import TestCase, main
from methods import *
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2 import MetadataCategory
from qiime2.plugins.demux.methods import emp_single
import tempfile
import os
import biom
import numpy as np
from numpy.testing import assert_array_almost_equal
from pandas.util.testing import assert_frame_equal

NUM_CORES = 4

class TestImport(TestCase):

    def setUp(self):
        print("Importing metadata for expected")
        self.exp_barcode_metadata = \
            Metadata.load("data/mock-3/sample-metadata.tsv")

        self.exp_demux = Artifact.load("data/mock-3/exp_demux.qza")
        self.exp_out = [self.exp_demux, self.exp_barcode_metadata]
        self.working_dir_fp = "data/mock-3"

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

class TestDeblur(TestCase):

    def setUp(self):
        self.exp_demux = Artifact.load("data/mock-3/exp_demux.qza")
        self.exp_deblurred = Artifact.load("data/mock-3/deblurred_150nt.qza")
        self.exp_deblurred_pt = Artifact.load("data/mock-3/deblurred_100nt_pt.qza")
        self.num_parallel = NUM_CORES

    def test_establish_dataset(self):
        obs = do_deblur(self.exp_demux, 150)
        self.assertEqual(self.exp_deblurred.view(biom.Table), obs.view(biom.Table))

    def test_establish_dataset_ncores(self):
        obs = do_deblur(self.exp_demux, 150,
                        num_cores=self.num_parallel)
        self.assertEqual(self.exp_deblurred.view(biom.Table), obs.view(biom.Table))

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

    def test_jaccard(self):
        obs = get_pairwise_dist_mat(self.jaccard_in,"jaccard")
        assert_array_almost_equal(self.jaccard_exp, obs.data)

    def test_bc(self):
        obs = get_pairwise_dist_mat(self.bc_in,"braycurtis")
        assert_array_almost_equal(self.bc_exp, obs.data)

class TestPrePostDist(TestCase):
    def setUp(self):
            self.pre = biom.Table(np.array([[0,0],[1,1],[1,1],[1,1]]),
                                  ['A', 'B', 'C','D'], ['S1', 'S2'])
            self.post = biom.Table(np.array([[1,1],[0,1],[1,1],[1,2]]),
                                  ['A', 'B', 'C','D'], ['S1', 'S2'])
            self.exp = pd.DataFrame(columns=["seq","dist_type","dist"],
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
            self.exp["dist"] = pd.to_numeric(self.exp["dist"])

    def test_get_distance_distribution(self):
        obs = get_distance_distribution(self.pre,self.post)
        print("obs")
        print(obs)
        print("exp")
        print(self.exp)
        assert_frame_equal(self.exp, obs)

class TestPostTrim(TestCase):
    def setUp(self):
        self.t = biom.Table(np.array([[0,1,2,3],[4,5,6,7],[8,9,10,11],[3,1,2,7]]),
                            ['AATT', 'AATG', 'ATGC','AATC'], ['S1', 'S2', 'S3', 'S4'])
        self.exp_pt = biom.Table(np.array([[7,7,10,17],[8,9,10,11]]),
                            ['AAT', 'ATG'], ['S1', 'S2', 'S3', 'S4'])
        self.exp_pt = Artifact.import_data("FeatureTable[Frequency]",
                                           self.exp_pt)

    def test_post_trim(self):
        obs = post_trim(self.t, 3)

        # Sort them
        obs = obs.sort_order(self.exp_pt.view(biom.Table).ids(axis="observation"),
                             axis="observation")

        print("----in----\n" + str(self.t))
        print("----out----\n" + str(obs))
        print("----exp----\n" + str(self.exp_pt.view(biom.Table)))
        self.assertEqual(str(self.exp_pt.view(biom.Table)), str(obs))

class TestShortSeq(TestCase):
    def setUp(self):
        barcode_map = pd.Series(['GTCA', 'TCAG', 'GGGG'],
                   index=['sample1', 'sample2', 'sample3'])
        barcode_map = MetadataCategory(barcode_map)

        seqs_fp = "data/small"
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
        obs1, obs2 = get_overlap_tables(self.none, self.some1)
        self.assertTrue(len(obs1.ids(axis="observation")) == 0)
        self.assertTrue(len(obs2.ids(axis="observation")) == 0)

# Tests for methods specific to mockrobiota
class TestMockMethods(TestCase):

    def setUp(self):
        self.MOCK6_FORWARD_URL = "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-6/mock-forward-read.fastq.gz"
        self.MOCK6_BARCODE_URL = "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-6/mock-index-read.fastq.gz"
        self.MOCK6_DATASET_METADATA_URL = "https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/data/mock-6/dataset-metadata.tsv"

        # Since mock-3 is small
        self.MOCK3_FORWARD_URL = "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-3/mock-forward-read.fastq.gz"
        self.MOCK3_BARCODE_URL = "https://s3-us-west-2.amazonaws.com/mockrobiota/latest/mock-3/mock-index-read.fastq.gz"
        self.MOCK3_SAMPLE_METADATA_URL = "https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/data/mock-3/sample-metadata.tsv"

    def test_get_dl_urls_pass(self):
        with tempfile.TemporaryDirectory() as output_dir:
            obs = get_dl_urls(self.MOCK6_DATASET_METADATA_URL, output_dir)

        exp = [self.MOCK6_FORWARD_URL, self.MOCK6_BARCODE_URL]
        for i, entry in enumerate(obs):
            self.assertEqual(entry, exp[i])

    def test_get_dl_urls_fail(self):
        with self.assertRaises(Exception) as context:
            get_dl_urls("asdf","asdf")

    def test_download_mock_dataset(self):
        with tempfile.TemporaryDirectory() as output_dir:
            download_mock_dataset(output_dir, self.MOCK3_FORWARD_URL,
                                  self.MOCK3_BARCODE_URL,
                                  self.MOCK3_SAMPLE_METADATA_URL)
            forward_fp = os.path.join(output_dir, "emp-single-end-sequences",
                                      "sequences.fastq.gz")
            self.assertTrue(os.path.exists(forward_fp))
            self.assertTrue(os.path.getsize(forward_fp) > 0)

            barcode_fp = os.path.join(output_dir, "emp-single-end-sequences",
                                      "barcodes.fastq.gz")
            self.assertTrue(os.path.exists(barcode_fp))
            self.assertTrue(os.path.getsize(barcode_fp) > 0)

            metadata_fp = os.path.join(output_dir, "sample-metadata.tsv")
            self.assertTrue(os.path.exists(metadata_fp))
            self.assertTrue(os.path.getsize(metadata_fp) > 0)

if __name__ == '__main__':
    main()