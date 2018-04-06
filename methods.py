from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.demux.methods import emp_single
from qiime2.plugins.quality_filter.methods import q_score
from qiime2.plugins.deblur.methods import denoise_16S
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
from skbio.diversity import beta_diversity
from skbio.stats.distance import mantel
import pandas as pd
import wget
import biom
import os
import scipy
import skbio.io
import numpy as np

def import_dataset(working_dir_fp, metadata_barcode_column,
                   rev_comp_barcodes_in=False,
                   rev_comp_mapping_barcodes_in=False):
    """
    Imports seqs as qiime artifact, demuxes them.
    Requires that fastq.gz files already be in
    working_dir_fp/emp-single-end-seqs and sample-metadata.tsv be in
    working_dir_fp

    :param working_dir_fp: filepath where sequences_url + barcodes_url file are
                           downloaded to
                           and put into a directory "emp-single-end-sequences".
                           Should also contain sample-metadata.tsv.
                           Ideally, this should be a
                           mock-<n> directory from when you clone the
                           mockrobiota github repo
                           should not end with "/"
    :param metadata_barcode_column: column header in sample-metadata.tsv that
                                    holds barcode data
    :param rev_comp_barcodes_i: param to emp_single for reversing barcode seqs
    :param rev_comp_mapping_barcodes_i: param to emp_single for reversing
                                        barcode seqs in metadata

    :return tuple of demuxed seqs, loaded metadata, None if fails
    """
    print("Importing seq data")
    seqs = Artifact.import_data("EMPSingleEndSequences", working_dir_fp +
                                "/emp-single-end-sequences")

    print("Loading metadata")
    barcode_metadata = Metadata.load(working_dir_fp + "/sample-metadata.tsv")

    print("Demuxing")
    demux, = emp_single(seqs,
                        barcode_metadata.get_category(metadata_barcode_column),
                        rev_comp_barcodes = rev_comp_barcodes_in,
                        rev_comp_mapping_barcodes =
                        rev_comp_mapping_barcodes_in)
    return demux, barcode_metadata


def do_deblur(demuxed_seqs, pre_trim_length, num_cores = 1):
    """
    Given demuxed sequences, quality filter,pre_trims them and returns the result

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



def post_trim(deblurred_biom, post_trim_length):
    """
    Trims a deblurred set of seqs
    :param deblurred_biom deblurred seqs as biom table
    :param post_trim_length length to trim to
    :return trimmed deblurred seqs as BIOM
    """
    print("Trimming post-demuxed seqs to {:d}".format(post_trim_length))
    post_trimmed_biom = \
        deblurred_biom.collapse(lambda i, m: i[:post_trim_length],
                                axis="observation", norm=False)

    return post_trimmed_biom

def get_distance_distribution(pre_table_overlap, post_table_overlap):
    """
    Given biom tables of overlapping reads, returns jaccard and bray curtis
    distances between matching reads.
    Two params should have exact same reads
    :param pre_table_overlap biom table for pre trimmed reads
    :param post_table_overlap biom table for post trimmed reads
    :return pandas dataframe with columns ["seq","dist_type","dist"]
    that hold sequence, type of distance and distance value
    """
    distance_functions = [('jaccard', scipy.spatial.distance.jaccard),
                          ('braycurtis', scipy.spatial.distance.braycurtis)]

    results = []
    pre_pa = pre_table_overlap.pa(inplace=False)
    post_pa = post_table_overlap.pa(inplace=False)
    for obs in pre_table_overlap.ids(axis='observation'):
        for fname, f in distance_functions:
            if(fname == "jaccard"):
                a = pre_pa.data(obs, axis='observation', dense=True)
                b = post_pa.data(obs, axis='observation', dense=True)
                print(obs)
                print("jaccard between " + str(a) + " and " + str(b) + " f=" + str(f(a,b)))
            else:
                a = pre_table_overlap.data(obs, axis='observation', dense=True)
                b = post_table_overlap.data(obs, axis='observation',
                                            dense=True)
            results.append((obs, fname, f(a, b))) # TODO append sample too?

    results = pd.DataFrame(results, columns = ["seq", "dist_type", "dist"])
    return results

def get_pairwise_dist_mat(deblur_biom, dist_type):
    if(dist_type == "jaccard"):
        deblur_biom = deblur_biom.pa(inplace=False)

    dist_mat = beta_diversity(dist_type,
                              deblur_biom.transpose().to_dataframe().as_matrix().astype("int64"),
                              ids = deblur_biom.ids(axis="sample"))
    return dist_mat


def get_overlap_tables(pre, post):
    """
    Takes in biom tables and returns the part of them that overlap in
    same order

    :param pre biom table of pre-trimmed deblurred seqs
    :param post biom table of post-trimmed deblurred seqs
    :return tuple of biom tables (pre_o,post_o) where each table only holds
            reads found in pre and post arguments, in same order as post
    """
    pre_ids = pre.ids(axis='observation')
    post_ids = post.ids(axis='observation')

    features_in_common = set(pre_ids) & set(post_ids)
    print(str(len(features_in_common)) + " reads in common when intersecting")
    pre_table_overlap = pre.filter(features_in_common, axis='observation',
                                   inplace=False)
    post_table_overlap = post.filter(features_in_common, axis='observation',
                                     inplace=False)

    # put the tables into the same order on both axes
    try:
        pre_table_overlap = pre_table_overlap.sort_order(post_table_overlap.ids(axis='observation'), axis='observation')
    except Error:
        print("Failed to sort tables, possibly because number of samples changed between pre and post")
        return

    return (pre_table_overlap, post_table_overlap)

def get_pre_post_distances(pre_bioms, post_bioms, trim_lengths):
    """For each otu, get distance between the otu in pre and post. Returns
    all distances in a pandas dataframe. Does jaccard and bray curtis

    Parameters
    ----------
    pre_bioms: array_like of qiime2.Artifacts type FeatureTable[Frequency]
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_bioms
    post_bioms: array_like of qiime2 artifacts
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_bioms
    trim_lengths: array_like
        Trim lengths in descending order, should correspond to other arguments

    Returns
    -------
    Pandas dataframe that holds results for each pre-post mantel test
    """
    if(not (len(pre_bioms) == len(post_bioms) == len(trim_lengths))):
        raise ValueError("Length of 3 arguments lists should be same\n"
                         "pre: {}, post: {}, lengths: {}".format(len(pre_bioms),
                                                                 len(post_bioms),
                                                                 len(trim_lengths)))

    pre_overlaps = []
    post_overlaps = []
    all_dists = pd.DataFrame()
    for i in range(len(pre_bioms)):
        # pre-post distances
        pre_overlap_biom, post_overlap_biom = \
            get_overlap_tables(pre_bioms[i], post_bioms[i])

        pre_overlaps.append(pre_overlap_biom)
        post_overlaps.append(post_overlap_biom)

        dists = get_distance_distribution(pre_overlap_biom,
                                                  post_overlap_biom)
        dists["length"] = trim_lengths[i]
        all_dists = all_dists.append(dists)

    return all_dists

def get_pairwise_diversity(pre_bioms, post_bioms, trim_lengths):
    """For each pre-post pair, gets the pairwise distance matrix of each
    sequence set and does a mantel test between pre and post pariwise distance
    matrices using both jaccard and bray-curtis metrics

    Parameters
    ----------
    pre_bioms: array_like of qiime2.Artifacts type FeatureTable[Frequency]
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_bioms
    post_bioms: array_like of qiime2 artifacts
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_bioms
    trim_lengths: array_like
        Trim lengths in descending order, should correspond to other arguments

    Returns
    -------
    Pandas dataframe that holds results for each pre-post mantel test
    """
    if(not (len(pre_bioms) == len(post_bioms) == len(trim_lengths))):
        raise ValueError("Length of 3 arguments lists should be same\n"
                         "pre: {}, post: {}, lengths: {}".format(len(pre_bioms),
                                                                 len(post_bioms),
                                                                 len(trim_lengths)))

    cols = ["trim_length", "dist_type", "r", "pval", "nsamples"]
    pairwise_diversity = pd.DataFrame(columns=cols)
    for i in range(len(pre_bioms)):
        # pairwise distance matrices
        pre_biom = pre_bioms[i]
        post_biom = post_bioms[i]

        pre_d_j = get_pairwise_dist_mat(pre_biom, "jaccard")
        post_d_j = get_pairwise_dist_mat(post_biom, "jaccard")
        r, p, nsamp = mantel(pre_d_j, post_d_j)
        pairwise_diversity = pairwise_diversity.append(dict(zip(cols, [trim_lengths[i], "jaccard", r, p, nsamp])),
                                                       ignore_index=True)
        pre_d_bc = get_pairwise_dist_mat(pre_biom, "braycurtis")
        post_d_bc = get_pairwise_dist_mat(post_biom, "braycurtis")
        r, p, nsamp = mantel(pre_d_bc, post_d_bc)
        pairwise_diversity = pairwise_diversity.append(dict(zip(cols, [trim_lengths[i], "braycurtis", r, p, nsamp])),
                                                       ignore_index=True)

    pairwise_diversity["r_sq"] = pairwise_diversity["r"]**2
    return pairwise_diversity

def get_shortest_seq(demuxed):
    """
    Given a qiime artifact of demuxed reads, returns the length of the read
    :param qiime2 artifact of demuxed reads
    :return length of shortest read
    """
    directory = demuxed.view(SingleLanePerSampleSingleEndFastqDirFmt)

    lengths = {}
    for file_name, format_fp in directory.sequences.iter_views(FastqGzFormat):
        seqs = skbio.io.read(str(format_fp), format='fastq', phred_offset=33)
        lengths[str(format_fp)] = min([len(s) for s in seqs])

    return min(lengths.values())

# The following methods are specific to mockrobiota dataset
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