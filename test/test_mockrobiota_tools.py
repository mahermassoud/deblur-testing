from unittest import TestCase, main
from mockrobiota_tools import *
import tempfile
import os

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

if __name__ == "__main":
    main()
