import click
from qiime2 import Artifact, Metadata
from qiime2.plugins.demux.methods import emp_single
import methods
import seaborn as sns
import matplotlib.pyplot as plt
import biom
import glob
import os
import re
import pandas as pd

@click.command()
@click.option('-i','--input-fp', required=True,
              type=click.Path(exists=True),
              help="Path to folder that contains sequences and barcodes")
@click.option('-m','-metadata', required=True,
              type=click.Path(exists=True),
              help="Path to sample metadata file")
@click.option('-metadata_bc_col', required=True,
              help="Name of column in metadata file that holds barcodes")
@click.option('-rev_bc', is_flag=True, help="Will reverse barcodes")
@click.option('-rev_map_bc', is_flag=True, help="Will reverse mapping barcodes")
@click.option('-o', '--output-fp',  type=click.Path(), default = None,
              help="Path to where demuxed qza is saved. Must specify filename. qza extension optional. Does not save if not supplied")
# TEST with rdemux -seqs=test/analysis_testing_wd/mock-3/emp-single-end-sequences/
# -metadata=test/analysis_testing_wd/mock-3/sample-metadata.tsv -out=rdemux_out
# -metadata_bc_col=BarcodeSequence
def demux(input_fp, metadata, metadata_bc_col, rev_bc, rev_map_bc, output_fp):
    """Imports data and runs demux on it

    Parameters
    ----------
    input_fp: str
        Path to folder that contains sequences and barcodes. MUST be
        named sequences.fastq.gz and barcodes.fastq.gz
    metadata: str
        Path to metadata file that contains barcode sequences. See
        "sample-metadata.tsv" in mockrobiota datasets for examples.
    metadata_bc_col: str
        Name of column in metadata file that holds the barcode sequences
    rev_bc: bool
        Whether to reverse barcodes
    rev_map_bc: bool
        Whether to reverse mapping barcodes
    output_fp: str, optional or None
        Path to where we output demuxed qza file. Does not save if None

    Returns
    -------
    Artifact
        demuxed sequences
    Metadata
        Associated metadata
    """
    if(input_fp is None or metadata is None or metadata_bc_col is None):
        click.echo("Run \'rdemux --help\' flag to see correct usage")
        return

    click.echo("Importing seq data from " + input_fp)
    art = Artifact.import_data("EMPSingleEndSequences", input_fp)

    print("Loading metadata from " + metadata)
    barcode_metadata = Metadata.load(metadata)

    print("Demuxing")
    demux, = emp_single(art,
                        barcode_metadata.get_category(metadata_bc_col),
                        rev_comp_barcodes = rev_bc,
                        rev_comp_mapping_barcodes = rev_map_bc)

    if(output_fp is None):
        click.echo("Not saving demux output")
    else:
        demux.save(output_fp)

    return demux, barcode_metadata

@click.command()
@click.option('-i', '--input-fp',
              type=click.Path(exists=True), required = True,
              help="Path to demuxed qza")
@click.option('-l', '--trim-length', type=click.INT, default=100,
              help='Trim length, default set to max possible length. '
                   'Supply if possible since calcuation is expensive')
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-o', '--output-fp', type=click.Path(), default = None,
              help='Path to output deblurred qza files, optional')
@click.option('-nc', '--num-cores', type=click.INT, default=1,
              help="Number of cores to parallelize deblur")
def pre_trims(input_fp, trim_length, trim_incr,
              num_trims, output_fp, num_cores):
    """Quality filters and then pre_trims sequences to various pre-trim lengths
    Saves qza's if specified. With naming format "deblurred_pre_<length>nt.qza

    Parameters
    ----------
    input_fp: path
        Path to qza of demuxed sequences
    trim_length: int, optional
        Length to trim to. If not supplied, longest possible length is used.
        This takes a while so do supply if possible
    trim_incr: int, optional
        Percent amount to decrement by.
    num_trims: int, optional
        Number of different lengths to trim to. Each trim_incr % less.
    output_fp: path, optional
        Path to output deblurred qza files
    num_cores: int, optional
        Number of cores to parallelize deblur

    Returns
    -------
    list of length trim_lengths of deblurred seq artifacts
    """
    click.echo("Importing seq data from " + input_fp)
    input_artifact = Artifact.load(input_fp)
    return pre_trims_art(input_artifact, trim_length, trim_incr, num_trims,
                         output_fp, num_cores)

def pre_trims_art(input_artifact, trim_length= 100, trim_incr = 10,
                  num_trims = 5, output_fp = None, num_cores = 1):
    """Quality filters and then pre_trims sequences to various pre-trim lengths
    Saves qza's if specified. With naming format "deblurred_pre_<length>.qza

    Parameters
    ----------
    input_artifact: qiime2.Artifact
        Qiime2 artifact of the demuxed seqs if we are using the artifact api.
        Either this or input_fp needs to be supplied.
    trim_length: int, optional
        Length to trim to. If not supplied, longest possible length is used.
        This takes a while so do supply if possible
    trim_incr: int, optional
        Percent amount to decrement by.
    num_trims: int, optional
        Number of different lengths to trim to. Each trim_incr % less.
    output_fp: path, optional
        Path to output deblurred qza files
    num_cores: int, optional
        Number of cores to parallelize deblur

    Returns
    -------
    list of length trim_lengths of deblurred seq artifacts
    """
    if(trim_length is None):
        click.echo("Determining max possible trim length")
        trim_length = methods.get_shortest_seq(input_artifact)

    trim_lengths, percents = calculate_trim_lengths(trim_length, trim_incr, num_trims)

    deblurreds = []
    for l in trim_lengths:
        db_out = methods.do_deblur(input_artifact, l, num_cores=num_cores)
        deblurreds.append(db_out)
        if(output_fp is not None):
            db_out.save("deblurred_pre_" + str(l) + ".qza")

    return deblurreds

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, dir_okay=False),
              required=True,
              help="Path to deblurred qza that we are post-trimming")
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-cf', '--clps-fp', type=click.Path(), default="collapse.csv",
              help="Output path for collapse data tsv file")
@click.option('-o', '--output-fp', type=click.Path(dir_okay=False), default = None,
              help='Path to output post-trimmed qza files, optional')
def post_trims(input_fp, trim_incr, num_trims, clps_fp, output_fp):
    click.echo("Importing seq data from " + input_fp)
    input_artifact = Artifact.load(input_fp)

    print(input_fp)
    print(input_artifact)
    return post_trims_art(clps_fp, input_artifact, trim_incr, num_trims,
                          output_fp)


def post_trims_art(clps_fp, input_artifact = None, trim_incr = 10,
               num_trims = 5, output_fp = None):
    """Post trims to various specified lengths.
    Saves qza's if specified. With naming format "deblurred_pt_<length>.qza
    eg. If input is length 100, trim_incr=10 and num_trims=5, post trims to
    lengths 100, 90, 80, 70, 60

    Parameters
    ----------
    clps_fp: str
        Path to where we write out collapse data file
    input_fp: str
        Path to qza of demuxed sequences. Either this or input_artifact
        must be supplied.
    input_artifact: qiime2.Artifact
        Qiime2 artifact of the deblurred seqs if we are using the artifact api.
    trim_incr: int, optional
        Percent amount to decrement by.
    num_trims: int, optional
        Number of different lengths to trim to. Each trim_incr % less.
    output_fp: path, optional
        Path to output deblurred qza files

    Returns
    -------
    list of length trim_lengths of post_trimmed seqs as biom tables
    """
    print(input_artifact)
    input_biom = input_artifact.view(biom.Table)
    otus = input_biom.ids(axis="observation")
    trim_length = len(otus[0])
    for otu in otus:
        if(len(otu) != trim_length):
            raise ValueError("Input table reads are not all same length. Invalid")
            return

    trim_lengths, percent = calculate_trim_lengths(trim_length, trim_incr, num_trims)

    pt_bioms = []
    for l in trim_lengths:
        pt_biom = methods.post_trim(input_biom, l)
        pt_bioms.append(pt_biom)
        if(output_fp is not None):
            pt_artifact = Artifact.import_data("FeatureTable[Frequency]",
                                               pt_biom)
            pt_artifact.save("deblurred_pt_" + str(l) + ".qza")

    clps = methods.get_collapse_counts(pt_bioms)
    clps.to_csv(clps_fp, index=False)

    return pt_bioms

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              required=True,
              help="Path to directory holding pre and post qza's. In naming "
              "format output by post_trims and pre_trims. Required that each "
              "post-trim has corresponding pre-trim")
@click.option('-o', '--output-fp',type=click.Path(file_okay=False,exists=True),
              default = None, required=False, help='Path to output csv files')
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
def analysis(input_fp, output_fp, trim_incr, num_trims):

    pres = dict()
    res = [f for f in os.listdir(input_fp)
           if re.match('deblurred_pre_\d+.qza', f)]
    for f in res:
        fm = f.replace("_",".")
        fm = fm.split(".")
        length = int(fm[len(fm)-2])
        artifact = Artifact.load(input_fp + "/" + f)
        pres[length] = artifact
    pre_artifacts = [pres[x] for x in sorted(pres.keys())]
    pre_artifacts.reverse()

    posts = dict()
    res = [f for f in os.listdir(input_fp)
           if re.match('deblurred_pt_\d+.qza', f)]
    for f in res:
        fm = f.replace("_",".")
        fm = fm.split(".")
        length = int(fm[len(fm)-2])
        artifact = Artifact.load(input_fp + "/" + f)
        posts[length] = artifact
    post_artifacts = [posts[x] for x in sorted(posts.keys())]
    post_artifacts.reverse()

    return analysis_art(pre_artifacts, post_artifacts, trim_incr, num_trims,
                        output_fp)


def analysis_art(pre_artifacts, post_artifacts, trim_incr = 10, num_trims = 5,
                 output_fp = None):
    """Returns analysis data on pre/post artifacts

    Parameters
    ----------
    pre_artifacts: array_like of qiime2.Artifacts type FeatureTable[Frequency]
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_artifacts
    post_artifacts: array_like of qiime2 artifacts
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_artifacts
    trim_incr: int, optional
        Percent amount to decrement by. Should correspond to input artifacts
    num_trims: int, optional
        Number of different lengths to trim to. Each trim_incr % less.
        Should correspond to input artifacts
    output_fp: str
        Path to directory which output tsv's are saved, should not end with /

    Returns
    -------
    pandas DataFrame of pariwise mantel test data
    pandas DataFrame of distance between pre and post matching otus
    pandas DataFrame of sOTU and sample counts
    pandas DataFrame of changes in reads per sample from pre to post (pre-post)
    """
    pre_bioms = [art.view(biom.Table) for art in pre_artifacts]
    post_bioms = [art.view(biom.Table) for art in post_artifacts]

    otus = pre_bioms[0].ids(axis="observation")
    length = len(otus[0])
    trim_lengths, prc = calculate_trim_lengths(length, trim_incr,num_trims)

    pairwise_mantel = methods.get_pairwise_diversity_data(pre_bioms, post_bioms,
                                                          trim_lengths)
    pre_post, pre_overlaps, post_overlaps = \
        methods.get_pre_post_distance_data(pre_bioms, post_bioms, trim_lengths)

    counts, read_changes = methods.get_count_data(pre_bioms, pre_overlaps,
                                                  post_bioms, post_overlaps,
                                                  trim_lengths)

    if(output_fp is None):
        pairwise_mantel.to_csv(output_fp + "/pairwise_mantel.csv", index=False)
        pre_post.to_csv(output_fp + "/pre_post.csv", index=False)
        counts.to_csv(output_fp + "/counts.csv", index=False)
        read_changes.to_csv(output_fp + "/read_changes.csv", index=False)

    return pairwise_mantel, pre_post, counts, read_changes

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              required=True,
              help="Path to directory holding csv files made by analysis")
@click.option('-o', '--output-fp',type=click.Path(file_okay=False,exists=True),
              default = None, required=False, help='Path to output csv files')
def do_plots(input_fp, output_fp):
    pairwise_mantel = pd.read_csv(input_fp + "/pairwise_mantel.csv")
    pre_post = pd.read_csv(input_fp + "/pre_post.csv")
    counts = pd.read_csv(input_fp + "/counts.csv")
    read_changes = pd.read_csv(input_fp + "/read_changes.csv")

    return plot_pd(pairwise_mantel, pre_post, counts, read_changes, output_fp)

def plot_pd(pairwise_mantel, pre_post, counts, read_changes, output_fp = None):
    """Plots data generated by analysis()

    Parameters
    ----------
    pairwise_mantel: pandas DataFrame
    pre_post: pandas DataFrame
    counts: pandas DataFrame
    read_changes: pandas DataFrame
        all in formart output by analysis()
    output_fp: str
        directory where to deposit plot png files
    Returns
    -------
    pyplot figure of each respective plot
    """
    mantel_plot = sns.lmplot(x="trim_length", y="r_sq", row="dist_type", hue="dist_type",
                             data=pairwise_mantel, ci=None, fit_reg=False)
    mantel_plot.set(xlabel="Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/pairwise_mantel.png")
    plt.figure()

    pp_plot = sns.boxplot(x="length", y="dist", data=pre_post, hue="dist_type")
    pp_plot.set(xlabel="trim_length", title="Distribution of distances between matching pre/post otus")

    if output_fp is not None:
        plt.savefig(output_fp + "/pre_post.png")
    plt.figure()

    plt.plot(counts["trim_length"], counts["sOTU_overlap_count"])
    plt.plot(counts["trim_length"], counts["sOTU_unique_pre"])
    plt.plot(counts["trim_length"], counts["sOTU_unique_post"])
    plt.legend(counts.columns[1:], loc='best')
    plt.title("sOTU counts")
    plt.xlabel("Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/counts.png")
    plt.figure()

    cols = list(read_changes.columns)
    for colname in read_changes.columns[1:]:
        plt.plot(read_changes["trim_length"], read_changes[colname])
    plt.legend(read_changes.columns[1:], loc='best')
    plt.title("pre - post: reads-per-sample count")
    plt.xlabel("Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/read_changes.png")

    #return mantel_plot, pp_plot, c_plot, rc_plot

def calculate_trim_lengths(length, trim_incr, num_trims):
    """Returns list of lengths we will trim to. Each trim_incr percent less
    eg. if length=100, trim_incr=10, num_trims=3, outputs [100,90,80]

    Parameters
    ----------
    length: int
        Max length
    trim_incr: int
        Percent amount to decrement by.
    num_trims: int
        Number of different lengths to trim to. Each trim_incr % less.

    Returns
    -------
    array_like
        actual trim lengths
    array_like
        trim lengths as % of max length
    """
    trim_lengths = []
    percents = []
    percent = 100
    for i in range(num_trims):
        trim_lengths.append(int(length * (percent/100)))
        percents.append(percent)
        percent = percent - trim_incr

    return trim_lengths, percents


# High priority
# TODO integration tests
# TODO demux unit tests
# TODO how to plot collapse count
# TODO try with normalized data
# TODO test against different environments eg. fecal, skin, soil
# TODO more metrics!! eg. seq depth
    # eg.
    # Look @ distribution of difference-per-feature
    # plot taxa that got dropped out
    # look @ top collapsed features
    # sequencing depth

# Low priority
# TODO input validity checks
# TODO saves to wd regardless of -o in analysis
# TODO fp with / still works
# TODO do eg.s in python console format
# TODO get rid of percents business

# Ideas
# TODO append sample name in pre_post

