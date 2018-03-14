from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.demux.methods import emp_single
from qiime2.plugins.quality_filter.methods import q_score
from qiime2.plugins.deblur.methods import denoise_16S
import pandas as pd
import wget
import biom
import os

"""
This module holds methods used to test deblur post-trim vs pre-trim results
"""
def download_mock_dataset(working_dir_fp, sequences_url, barcodes_url,
                          metadata_url):
    """
    Downloads all the files we need from a mock dataset setup so qiime2 can
    work with them. Puts all files inside
    working_dir_fp + "/emp-single-end-sequences"

    :param working_dir_fp: filepath where sequences_url + barcodes_url file are downloaded to
                           and put into a directory "emp-single-end-sequences". Ideally, this should be a
                           mock-<n> directory from when you clone the mockrobiota github repo
                           should not end with "/"
    :param sequences_url: URL from which sequences are downloaded from using wget
    :param barcodes_url: URL from which barcodes are downloaded from using wget
    :param metadata_url: URL from which we get sample-metadata.tsv file
    """
    # setup files so qiime can work with them
    output_dir = working_dir_fp + "/emp-single-end-sequences"
    try:
        os.mkdir(output_dir)
    except FileExistsError:
       pass

    print("Downloading " + sequences_url)
    wget.download(sequences_url, output_dir + "/sequences.fastq.gz")
    print("Downloading " + barcodes_url)
    wget.download(barcodes_url, output_dir + "/barcodes.fastq.gz")
    print("Downloading " + metadata_url)
    wget.download(metadata_url, working_dir_fp)


def import_mock_dataset(working_dir_fp, metadata_barcode_column):
    """
    Imports seqs as qiime artifact, demuxes them.
    Requires that fastq.gz files already be in
    working_dir_fp/emp-single-end-seqs and sample-metadata.tsv be in
    working_dir_fp

    :param working_dir_fp: filepath where sequences_url + barcodes_url file are downloaded to
                           and put into a directory "emp-single-end-sequences".
                           Should also contain sample-metadata.tsv.
                           Ideally, this should be a
                           mock-<n> directory from when you clone the mockrobiota github repo
                           should not end with "/"
    :param metadata_barcode_column: column header in sample-metadata.tsv that holds barcode data

    :return tuple of demuxed seqs, loaded metadata, None if fails
    """

    # TODO different if we are doing paired reads
    print("Importing seq data")
    seqs = Artifact.import_data("EMPSingleEndSequences", working_dir_fp +
                                "/emp-single-end-sequences")

    print("Loading metadata")
    barcode_metadata = Metadata.load(working_dir_fp + "/sample-metadata.tsv")

    print("Demuxing")
    demux, = emp_single(seqs,
                        barcode_metadata.get_category(metadata_barcode_column))
    return demux, barcode_metadata


def establish_dataset(demuxed_seqs, pre_trim_length, num_cores = 1):
    """
    Given demuxed sequences, deblurs them and returns the result

    :param demuxed_data: qiime artifact of demuxed data
    :param pre_trim_length: length that we want to trim sequences to before
        deblur is ran

    :return FeatureTable[Frequency] of deblured data. Post-trimmed to length
    post_trim_length if post_trim_length is not None
    """
    print("Quality filtering (with default params)")
    demuxed_qfiltered, demuxed_qf_stats = q_score(demuxed_seqs)

    print("Deblur-ing with trim length {:d}".format(pre_trim_length))
    deblurred, repseq, deblur_stats = \
        denoise_16S(demuxed_qfiltered, pre_trim_length,
                    hashed_feature_ids = False, jobs_to_start = num_cores)

    return deblurred

def get_dl_urls(dataset_metadata_url, working_dir_fp):
    """
    :param dataset_metadata_url is url where we can download dataset metadata like this
    https://github.com/caporaso-lab/mockrobiota/blob/master/data/mock-6/dataset-metadata.tsv
    :param working_dir_fp: filepath where sequences_url + barcodes_url file are downloaded to
                           and put into a directory "emp-single-end-sequences".
                           Should not end with "/"
    :return tuple of 2 strings, URLs of sequence and index files respectively
    """
    wget.download(dataset_metadata_url, working_dir_fp)
    metadata_fp = working_dir_fp  + "/" + os.path.basename(dataset_metadata_url)
    metadata = pd.read_table(metadata_fp, index_col = 0, dtype=str)

    seq_url = metadata.loc["raw-data-url-forward-read", "value"]
    barcode_url = metadata.loc["raw-data-url-index-read", "value"]

    return seq_url, barcode_url


def post_trim(deblurred_biom, post_trim_length):
    """
    Trims a deblurred set of seqs-- utility fn so we dont have to re-run establish_dataset()
    :param deblurred_biom deblurred seqs as biom table
    :param post_trim_length length to trim to
    :return trimmed deblurred seqs
    """
    print("Trimming post-demuxed seqs to {:d}".format(post_trim_length))
    post_trimmed_biom = \
        deblurred_biom.collapse(lambda i, m: i[:post_trim_length],
                                axis="observation", norm=False)

    post_trimmed = Artifact.import_data("FeatureTable[Frequency]",
                                        post_trimmed_biom);
    return post_trimmed