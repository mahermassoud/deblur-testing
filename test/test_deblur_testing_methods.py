import unittest
from unittest import TestCase, main
from deblur_testing_methods import get_dl_urls, download_mock_dataset, \
    import_mock_dataset, establish_dataset, post_trim
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.demux.methods import emp_single
from qiime2.plugins.quality_filter.methods import q_score
from qiime2.plugins.deblur.methods import denoise_16S
import tempfile
import os
import biom
import numpy as np

NUM_CORES = 4

class TestDeblurTestMethods(TestCase):

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


class TestImport(TestCase):

    def setUp(self):
        print("Importing metadata for expected")
        self.exp_barcode_metadata = \
            Metadata.load("data/mock-3/sample-metadata.tsv")

        self.exp_demux = Artifact.load("data/mock-3/exp_demux.qza")
        self.exp_out = [self.exp_demux, self.exp_barcode_metadata]
        self.working_dir_fp = "data/mock-3"

    def test_import_mock_dataset(self):
        obs = import_mock_dataset(self.working_dir_fp, "BarcodeSequence")
        for i, entry in enumerate(obs):
            try:
                self.assertEqual(entry, self.exp_out[i])
            except AssertionError as e:
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
        obs = establish_dataset(self.exp_demux, 150)
        self.assertEqual(self.exp_deblurred.view(biom.Table), obs.view(biom.Table))

    def test_establish_dataset_ncores(self):
        obs = establish_dataset(self.exp_demux, 150,
                                  num_cores=self.num_parallel)
        self.assertEqual(self.exp_deblurred.view(biom.Table), obs.view(biom.Table))

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
        print("----in----\n" + str(self.t))
        print("----out----\n" + str(obs.view(biom.Table)))
        print("----exp----\n" + str(self.exp_pt.view(biom.Table)))
        self.assertEqual(self.exp_pt.view(biom.Table), obs.view(biom.Table))


if __name__ == '__main__':
    main()