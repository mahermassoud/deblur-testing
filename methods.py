"""
Contains methods used in scripts.py for testing the difference between deblur
pre-trim and post-trim
"""
from qiime2 import Artifact
import time
from qiime2 import Metadata
from qiime2.plugins.demux.methods import emp_single
from qiime2.plugins.quality_filter.methods import q_score
from qiime2.plugins.deblur.methods import denoise_16S
from q2_types.per_sample_sequences import SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat
from skbio.diversity import beta_diversity
from skbio.stats.distance import mantel
import pandas as pd
import pathos.multiprocessing as mp
import biom
import os
import scipy
import skbio.io
import numpy as np
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
                        barcode_metadata.get_column(metadata_barcode_column),
                        rev_comp_barcodes = rev_comp_barcodes_in,
                        rev_comp_mapping_barcodes =
                        rev_comp_mapping_barcodes_in)
    return demux, barcode_metadata


def do_deblur(demuxed_seqs, pre_trim_length, num_cores = 1):
    """Given demuxed sequences, quality filters , deblurs them and returns the result

    Parameters
    ----------
    demuxed_data: qiime2.Artifact
        demuxed data of type EmpSingleEndSequences
    pre_trim_length: int
        length that we want to trim sequences to before deblur is ran
    num_cores: int
        number of jobs to start when deblurring

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



def post_trim(db_biom, length, partition_count=None):
    """Trims a deblurred set of seqs

    Parameters
    ----------
    db_biom: biom.Table
        deblurred seqs as biom table
    length: int
        length to trim to
    partition_count: int
        if not None, partitions table into partition_count tables and
        does post_trimming in parallel

    Returns
    -------
    biom.Table of trimmed, deblurred seqs with metadata for collapsed ids
    """
    print("Trimming post-demuxed seqs to {:d}, partition_count: {}".format(length, partition_count), flush=True)
    if partition_count is None:
        pt_biom = db_biom.collapse(lambda i, m: i[:length], axis="observation",
                                   norm=False, include_collapsed_metadata=True)
    else:
        print("Doing parallel post-trim, mp find {} cpu's".format(mp.cpu_count()), flush=True)
        sub_bioms, pool = partition_table(db_biom, partition_count)
        print("partition_Table() end at " + time.strftime("[%H:%M:%S]"), flush=True)

        args = [(sb, length) for sb in sub_bioms]
        pt_bioms = pool.map(single_post_trim, args)

        args = list(divide_chunks(pt_bioms, 2))
        while len(args) >= 1:
            print("len args {}\n\n".format(len(args)))

            pt_bioms = pool.map(intersect_bioms, args)

            args = list(divide_chunks(pt_bioms, 2))
            if args is None:
              break

        pt_biom = pt_bioms[0]

    return pt_biom

def partition_table(tbl, partition_count, drop=True):
    """
    Partitions a biom table into n parts by sample

    Parameters
    ----------
    tbl: biom.Table
        table we are partitioning
    partition_count: int
        number of partitions to partition into
    drop: bool
        whether to drop columns as we partition, set to True to be
        less memory-expensive

    Returns
    -------
    list of biom.Table of length partition_count
    processor pool we should re-use
    """
    print("partition_Table() starting at " + time.strftime("[%H:%M:%S]"), flush=True)
    sids = tbl.ids()
    id_parts = np.array_split(sids, partition_count)

    pool = mp.ProcessPool(nodes=len(id_parts), maxtasksperchild=1)
    args = [(tbl, x, drop) for x in id_parts]
    results = pool.map(index_tbl, args)
    return results, pool

def index_tbl(tbl_sids, drop=True):
    """
    Wrapper function to filter samples from a table

    Parameters
    ----------
    tbl_sids: 2-tuple
        biom table we are filtering, list of sample ids we want to keep
    drop: bool
        whether to drop from table as we filter

    Returns
    -------
    biom.Table with only the samples in tbl_sids[0]
    """
    tbl = tbl_sids[0]
    sids = tbl_sids[1]

    indexed = tbl.filter(sids, inplace=False)
    if drop:
        tbl.filter(sids, invert=True, inplace=True)
    return indexed

def make_biom(dat_obs_sample):
    """
    Wrapper function of biom.Table()

    Parameters
    ----------
    dat_obs_sample: 3-tuple
        data, sample_ids, observation_ids

    Returns
    -------
    biom.Table
    """
    return biom.Table(dat_obs_sample[0], dat_obs_sample[1], dat_obs_sample[2])

def single_post_trim(db_biom_length):
    """
    Post trims a biom table, called by post_trims() in parallel

    Parameters
    ----------
    db_biom_length: 2-tuple
        biom.Table, length we are post-trimming to

    Returns
    -------
    biom.Table that has been post-trimmed, meaning each observation was
    trimmed to argument length and matching observations were summed
    """
    db_biom = db_biom_length[0]
    length = db_biom_length[1]
    print("Trimming post-demuxed seqs to {:d}".format(length))
    pt_biom = db_biom.collapse(lambda i, m: i[:length], axis="observation",
                         norm=False, include_collapsed_metadata=True)

    return pt_biom

def divide_chunks(l, n):
    """
    Divides a list into sub lists of length n

    Parameters
    ----------
    l: array_like
        list we are dividing
    n: int
        size of chunks

    Returns
    -------
    Generator where each item is n elements of l. If there are not enough
    elements to make last chunk exactly n, it is as big as possible
    """
    if len(l) == 1:
      return None
    for i in range(0, len(l), n):
        yield l[i:i + n]

def intersect_bioms(bioms):
    """
    unions 2 bioms,

    Parameters
    ----------
    bioms: 2-tuple
        holds 2 bioms we want to union

    Returns
    -------
    biom.Table of union
    """
    if len(bioms) == 1:
        return bioms[0]
    else:
        return bioms[0].merge(bioms[1])


def get_collapse_counts(pt_bioms):
    """Fore each trim length and OTU , says how many times otu was collapsed to

    Parameters
    ----------
    bioms: array_like of biom.Table
        list of post trimmed bioms we are getting counts for. Must have the
        metadata column "collapsed_ids" on observation axis

    Returns
    -------
    pandas.DataFrame with columns ["seq", "length", "num_collapses"]
    """
    otu_col = []
    len_col = []
    counts_col = []
    for bt in pt_bioms:
        otus = bt.ids(axis="observation")
        length = len(otus[0])
        counts = [len(bt.metadata(id=otu, axis="observation")["collapsed_ids"])
                  for otu in otus]

        otu_col.extend(otus)
        len_col.extend([length] * len(otus))
        counts_col.extend(counts)


    return pd.DataFrame({"seq" : otu_col, "length" : len_col,
                         "num_collapses": counts_col})

def get_distance_distribution(pre_table_overlap, post_table_overlap,
                              by_sample=False):
    """Given biom tables of overlapping reads, returns jaccard and bray curtis
    distances between matching reads OR samples.
    Two params should have exact same reads and samples

    Parameters
    ----------
    pre_table_overlap: biom.Table
        pre trimmed reads
    post_table_overlap: biom.Table
        post trimmed reads
    by_sample: bool
        True if we want to take each sample as a vector and compare distances

    Returns
    -------
    pandas.DataFrame with columns ["seq","dist_type","dist"]
    that hold sequence, type of distance and distance value for otu pair
    """
    distance_functions = [('jaccard', scipy.spatial.distance.jaccard),
                          ('braycurtis', scipy.spatial.distance.braycurtis)]

    if(by_sample):
        axis="sample"
        columns = ["sample", "dist_type", "dist"]
    else:
        axis="observation"
        columns = ["seq", "dist_type", "dist"]

    results = []
    pre_pa = pre_table_overlap.pa(inplace=False)
    post_pa = post_table_overlap.pa(inplace=False)
    for obs in pre_table_overlap.ids(axis=axis):
        for fname, f in distance_functions:
            if(fname == "jaccard"):
                a = pre_pa.data(obs, axis=axis, dense=True)
                b = post_pa.data(obs, axis=axis, dense=True)
            else:
                a = pre_table_overlap.data(obs, axis=axis, dense=True)
                b = post_table_overlap.data(obs, axis=axis, dense=True)

            d_val = f(a, b)
            if np.isnan(d_val):
                print("d_val was NaN!")
                d_val = 0 # because f(0,0) returns NaN
            results.append((obs, fname, d_val))

    results = pd.DataFrame(results, columns=columns)
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

    print("starting beta_diversity")
    dist_mat = beta_diversity(dist_type,
                              deblur_biom.transpose().to_dataframe().as_matrix().astype("int64"), 
                              ids = deblur_biom.ids(axis="sample"))
    print("end beta_diversity")
    return dist_mat


def get_overlap_tables(pre, post):
    """Takes in biom tables and returns the part of them that overlap in
    same order. Tables must have same samples

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
    if len(features_in_common) == 0:
        raise ValueError("Pre and post do not have any features in common")

    pre_table_overlap = pre.filter(features_in_common, axis='observation',
                                   inplace=False)
    post_table_overlap = post.filter(features_in_common, axis='observation',
                                     inplace=False)

    # put the tables into the same order on both axes
    pre_table_overlap = pre_table_overlap.sort_order(post_table_overlap.ids(axis='observation'), axis='observation')
    pre_table_overlap = pre_table_overlap.sort_order(post_table_overlap.ids(axis='sample'), axis='sample')

    return (pre_table_overlap, post_table_overlap)

def get_pre_post_distance_data(pre_bioms, post_bioms, trim_lengths):
    """For each otu, get distance between the otu and samples in pre and post. Returns
    all distances in two pandas dataframes. Does jaccard and bray curtis.

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
    Pandas dataframe that holds results for each pre-post distance by otu
    Pandas dataframe that holds results for each pre-post distance by sample
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
    all_dists_sample = pd.DataFrame()
    print("len(pre_bioms): {}".format(len(pre_bioms)))
    print("len(post_bioms): {}".format(len(post_bioms)))
    for i in range(len(pre_bioms)):
        # pre-post distances
        pre_overlap_biom, post_overlap_biom = \
            get_overlap_tables(pre_bioms[i], post_bioms[i])

        pre_overlaps.append(pre_overlap_biom)
        post_overlaps.append(post_overlap_biom)

        print("len(pre_overlaps): {}".format(len(pre_overlaps)))
        print("len(post_overlaps): {}".format(len(post_overlaps)))

        dists = get_distance_distribution(pre_overlap_biom,
                                          post_overlap_biom)
        dists_sample = get_distance_distribution(pre_overlap_biom,
                                                 post_overlap_biom, by_sample=True)

        #print("i: {}, dists: {}".format(str(i), str(dists)))

        dists["length"] = trim_lengths[i]
        dists_sample["length"] = trim_lengths[i]
        all_dists = all_dists.append(dists)
        all_dists_sample = all_dists_sample.append(dists_sample)
        #print("all_dists:\n{}".format(str(all_dists)))

    print("final all_dists:\n{}".format(str(all_dists)))
    return all_dists, all_dists_sample, pre_overlaps, post_overlaps

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
    print("enter get_pairwise_diversity")
    np.seterr(all="raise")
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
    print("exit get_pairwise_diversity")
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
    "diff_otu" : num of otus pre-post

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
    count_data["diff_otu"] = \
        [num_sOTUs(pre_bioms[i]) - num_sOTUs(post_bioms[i])
         for i in range(len(pre_bioms))]

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
