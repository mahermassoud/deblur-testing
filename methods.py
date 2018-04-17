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
from collections import Counter

def import_dataset(working_dir_fp, metadata_barcode_column,
                   rev_comp_barcodes_in=False,
                   rev_comp_mapping_barcodes_in=False):
    """Imports seqs as qiime artifact, demuxes them.
    Requires that fastq.gz files already be in
    working_dir_fp/emp-single-end-seqs and sample-metadata.tsv be in
    working_dir_fp

    Parameters
    ----------
    working_dir_fp: str
        filepath where sequences_url + barcodes_url file are
        downloaded to and put into a directory "emp-single-end-sequences".
        Should also contain sample-metadata.tsv. Ideally, this should be a
        mock-<n> directory from when you clone the mockrobiota github repo
        should not end with "/"
    metadata_barcode_column: str
        column header in sample-metadata.tsv that holds barcode data
    rev_comp_barcodes_i: bool
        param to emp_single for reversing barcode seqs
    rev_comp_mapping_barcodes_i: bool
        param to emp_single for reversing barcode seqs in metadata

    Returns
    -------
    demuxed seqs,
    loaded metadata
    OR
    None if fails
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
    """Given demuxed sequences, quality filter,pre_trims them and returns the result

    Parameters
    ----------
    demuxed_data: qiime2.Artifact
        demuxed data of type EmpSingleEndSequences
    pre_trim_length: int
    length that we want to trim sequences to before deblur is ran

    Returns
    -------
    qiime2.Artifact of tpye FeatureTable[Frequency] of deblured data.
    """
    print("Quality filtering (with default params)")
    demuxed_qfiltered, demuxed_qf_stats = q_score(demuxed_seqs)

    print("Deblur-ing with trim length {}".format(str(pre_trim_length)))
    deblurred, repseq, deblur_stats = \
        denoise_16S(demuxed_qfiltered, pre_trim_length,
                    hashed_feature_ids = False, jobs_to_start = num_cores)

    return deblurred



def post_trim(db_biom, length):
    """Trims a deblurred set of seqs

    Parameters
    ----------
    db_biom: biom.Table
        deblurred seqs as biom table
    length: int
        length to trim to

    Returns
    -------
    biom.Table of trimmed,deblurred seqs
    """
    print("Trimming post-demuxed seqs to {:d}".format(length))
    pt_biom = db_biom.collapse(lambda i, m: i[:length], axis="observation",
                               norm=False, include_collapsed_metadata=True)

    return pt_biom

def get_collapse_counts(pt_bioms):
    """Fore each trim length and OTU , says how many times otu was collapsed to

    Parameters
    ----------
    bioms: array_like of biom.Table
        list of post trimmed bioms we are getting counts for. Must have the
        metadata column "collapsed_ids" on observation axis

    Returns
    -------
    pandas.DataFrame with columns ["otu", "length", "num_collapses"]
    """
    otu_col = []
    len_col = []
    counts_col = []
    for bt in pt_bioms:
        otus = bt.ids(axis="observation")
        length = len(otus[0])
        counts = [len(bt.metadata(id=otu, axis="observation")) for otu in otus]

        otu_col = otu_col + otus
        len_col.append([length] * len(otus))
        counts_col.append(counts)


    return pd.DataFrame({"otu" : otu_col, "length" : len_col,
                         "num_collapses": counts_col})

def get_distance_distribution(pre_table_overlap, post_table_overlap):
    """Given biom tables of overlapping reads, returns jaccard and bray curtis
    distances between matching reads.
    Two params should have exact same reads

    Parameters
    ----------
    pre_table_overlap: biom.Table
        pre trimmed reads
    post_table_overlap: biom.Table
        post trimmed reads

    Returns
    -------
    pandas.DataFrame with columns ["seq","dist_type","dist"]
    that hold sequence, type of distance and distance value for otu pair
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
            else:
                a = pre_table_overlap.data(obs, axis='observation', dense=True)
                b = post_table_overlap.data(obs, axis='observation',
                                            dense=True)
            results.append((obs, fname, f(a, b))) # TODO append sample too?

    results = pd.DataFrame(results, columns = ["seq", "dist_type", "dist"])
    return results

def get_pairwise_dist_mat(deblur_biom, dist_type):
    """Returns pairwise distance matrix for deblurred seqs

    Parameters
    ----------
    deblur_biom: biom.Table
        Sequences we want pairwise distances (by sample) for
    dist_type: str
        Distance metric we want. Usually "jaccard" or "braycurtis"

    Returns
    -------
    numpy matrix of pairwise distances
    """
    if(dist_type == "jaccard"):
        deblur_biom = deblur_biom.pa(inplace=False)

    dist_mat = beta_diversity(dist_type,
                              deblur_biom.transpose().to_dataframe().as_matrix().astype("int64"),
                              ids = deblur_biom.ids(axis="sample"))
    return dist_mat


def get_overlap_tables(pre, post):
    """Takes in biom tables and returns the part of them that overlap in
    same order

    Parameters
    ----------
    pre: biom.Table
        pre-trimmed deblurred seqs
    post: biom.Table
        post-trimmed deblurred seqs

    Returns
    -------
    biom.Table of pre reads that overlap with post
    biom.Table of post reads that overlap with pre
    """
    pre_ids = pre.ids(axis='observation')
    post_ids = post.ids(axis='observation')

    features_in_common = set(pre_ids) & set(post_ids)
    #print(str(len(features_in_common)) + " reads in common when intersecting")
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

def get_pre_post_distance_data(pre_bioms, post_bioms, trim_lengths):
    """For each otu, get distance between the otu in pre and post. Returns
    all distances in a pandas dataframe. Does jaccard and bray curtis

    Parameters
    ----------
    pre_bioms: array_like of biom.Table
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_bioms
    post_bioms: array_like of biom.Table
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_bioms
    trim_lengths: array_like
        Trim lengths in descending order, should correspond to other arguments

    Returns
    -------
    Pandas dataframe that holds results for each pre-post mantel test
    array_like of biom.Table of overlapping otu's found in pre
    array_like of biom.Table of overlapping otu's found in post
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

    return all_dists, pre_overlaps, post_overlaps

def get_pairwise_diversity_data(pre_bioms, post_bioms, trim_lengths):
    """For each pre-post pair, gets the pairwise distance matrix of each
    sequence set and does a mantel test between pre and post pariwise distance
    matrices using both jaccard and bray-curtis metrics

    Parameters
    ----------
    pre_bioms: array_like of biom.Table
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_bioms
    post_bioms: array_like of biom.Table
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
    p_div = pd.DataFrame(index=range(2*len(pre_bioms)), columns=cols)
    j = 0
    for i in range(len(pre_bioms)):
        # pairwise distance matrices
        pre_biom = pre_bioms[i]
        post_biom = post_bioms[i]

        pre_d_j = get_pairwise_dist_mat(pre_biom, "jaccard")
        post_d_j = get_pairwise_dist_mat(post_biom, "jaccard")
        r, p, nsamp = mantel(pre_d_j, post_d_j)
        p_div.iloc[j] = [trim_lengths[i], "jaccard", r, p, nsamp]
        j += 1

        pre_d_bc = get_pairwise_dist_mat(pre_biom, "braycurtis")
        post_d_bc = get_pairwise_dist_mat(post_biom, "braycurtis")
        r, p, nsamp = mantel(pre_d_bc, post_d_bc)
        p_div.iloc[j] = [trim_lengths[i], "braycurtis", r, p, nsamp]

    p_div["r_sq"] = p_div["r"]**2
    return p_div

def get_count_data(pre_bioms, pre_overlaps, post_bioms, post_overlaps,
                  trim_lengths):
    """
    Parameters
    ----------
    pre_bioms: array_like of biom.Table
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_bioms
    pre_overlaps: array_like of biom.Table
        biom.Table of pre data after pre/post intersection. Same length as
        post_overlaps
    post_bioms: array_like of biom.Table
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_bioms
    post_overlaps: array_like of biom.Table
        biom.Table of post data after pre/post intersection. Same length as
        pre_overlaps
    trim_lengths: array_like
        Trim lengths in descending order, should correspond to other arguments

    Returns
    -------
    pandas DataFrame with columns:
    "trim_length" : trim_length,
    "sotu_overlap_count" : num overlapping sOTUs,
    "sotu_unique_pre" : number of sOTUs unique to pre,
    "sotu_unique_post" : number of sOTUs unique to post,
    "num_samples_pre" : num samples in pre_trim,
    "num_samples_post" : num samples in post_trim,
    "change_num_samples" : change in number of samples,

    pandas DataFrame with columns:
    "trim_length",
    <a column for each sample with same sample name as input">,
    values are change in reads per sample for each trim lenght, pre-post
    """
    if(not (len(pre_overlaps) == len(post_overlaps) == len(trim_lengths))):
        raise ValueError("Length of 3 arguments lists should be same\n"
                         "pre: {}, post: {}, lengths: {}"
                         .format(len(pre_overlaps),
                                 len(post_overlaps),
                                 len(trim_lengths)))

    count_data = pd.DataFrame()
    count_data["trim_length"] = trim_lengths
    count_data["sOTU_overlap_count"] = \
        [num_sOTUs(table) for table in pre_overlaps]
    count_data["sOTU_unique_pre"] = \
        [num_sOTUs(pre_bioms[i]) - num_sOTUs(pre_overlaps[i])
         for i in range(len(pre_overlaps))]
    count_data["sOTU_unique_post"] = \
        [num_sOTUs(post_bioms[i]) - num_sOTUs(post_overlaps[i])
         for i in range(len(pre_overlaps))]

    # TODO get rid of this because its always 0?
    #count_data["num_samples_pre"] = [num_samples(table) for table in pre_bioms]
    #count_data["num_samples_post"] = [num_samples(table) for table in post_bioms]
    #count_data["change_num_samples"] = count_data["num_samples_pre"] - \
    #                                   count_data["num_samples_post"]

    change_reads_per_sample = pd.DataFrame()
    change_reads_per_sample["trim_length"] = trim_lengths
    pre_sums = pd.DataFrame([total_read_counts(tbl) for tbl in pre_bioms],
                            columns = pre_bioms[0].ids(axis="sample"))
    post_sums = pd.DataFrame([total_read_counts(tbl) for tbl in post_bioms],
                             columns = post_bioms[0].ids(axis="sample"))
    delta = pre_sums - post_sums
    change_reads_per_sample = \
        pd.concat([change_reads_per_sample, delta], axis=1)

    return count_data, change_reads_per_sample

def num_sOTUs(biom_table):
    """Returns number of sOTUs in a biom.Table"""
    return biom_table.length(axis="observation")

def num_samples(biom_table):
    return biom_table.length(axis="sample")

def total_read_counts(biom_table):
    """Returns a list read counts per sample for each sample in biom_table"""
    return biom_table.sum(axis="sample")

def get_shortest_seq(demuxed):
    """Given a qiime artifact of demuxed reads, returns the length of the read

    Parameters
    ----------
    demuxed: qiime2.Artifact
        demuxed reads

    Returns
    -------
    int length of shortest read
    """
    directory = demuxed.view(SingleLanePerSampleSingleEndFastqDirFmt)

    lengths = {}
    for file_name, format_fp in directory.sequences.iter_views(FastqGzFormat):
        seqs = skbio.io.read(str(format_fp), format='fastq', phred_offset=33)
        lengths[str(format_fp)] = min([len(s) for s in seqs])

    return min(lengths.values())

###########################################################
# The following methods are specific to mockrobiota dataset
###########################################################
def get_dl_urls(dataset_metadata_url, working_dir_fp):
    """
    Parameters
    ----------
    dataset_metadata_url: str
        url where we can download dataset metadata like this
        https://github.com/caporaso-lab/mockrobiota/blob/master/data/mock-6/dataset-metadata.tsv
    working_dir_fp: str
        filepath where sequences_url + barcodes_url file are downloaded to
        and put into a directory "emp-single-end-sequences".
        Should not end with "/"

    Returns
    -------
    str url to download sequences from
    str url to download barcodes from
    """
    wget.download(dataset_metadata_url, working_dir_fp)
    metadata_fp = working_dir_fp  + "/" + os.path.basename(dataset_metadata_url)
    metadata = pd.read_table(metadata_fp, index_col = 0, dtype=str)

    seq_url = metadata.loc["raw-data-url-forward-read", "value"]
    barcode_url = metadata.loc["raw-data-url-index-read", "value"]

    return seq_url, barcode_url

def download_mock_dataset(working_dir_fp, sequences_url, barcodes_url,
                          metadata_url):
    """Downloads all the files we need from a mock dataset setup so qiime2 can
    work with them. Puts all files inside
    working_dir_fp + "/emp-single-end-sequences"

    Parameters
    ----------
    working_dir_fp: str
        filepath where sequences_url + barcodes_url file are downloaded to
        and put into a directory "emp-single-end-sequences". Ideally, a
        mock-<n> directory from when you clone the mockrobiota github repo
        should not end with "/"
    sequences_url: str
        URL from which sequences are downloaded from using wget
    barcodes_url: str
        URL from which barcodes are downloaded from using wget
    metadata_url: str
        URL from which we get sample-metadata.tsv file
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