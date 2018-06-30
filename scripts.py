import click
from qiime2 import Artifact, Metadata
from qiime2.plugins.demux.methods import emp_single
from random import choice
import methods
import seaborn as sns
import matplotlib.pyplot as plt
import biom
import time
import os
import re
import pandas as pd
import numpy as np

@click.command()
@click.option('-i','--input-fp', required=True,
              type=click.Path(exists=True, file_okay=False),
              help="Path to folder that contains sequences and barcodes")
@click.option('-m','-metadata', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Path to sample metadata file")
@click.option('-mbc', '-metadata_bc_col', required=True,
              help="Name of column in metadata file that holds barcodes")
@click.option('-rev_bc', is_flag=True, help="Will reverse barcodes")
@click.option('-rev_map_bc', is_flag=True, help="Will reverse mapping barcodes")
@click.option('-o', '--output-fp',  type=click.Path(dir_okay=False), default = None,
              help="Path to where demuxed qza is saved. Must specify filename. qza extension optional. Does not save if not supplied")
def do_demux(input_fp, metadata, metadata_bc_col, rev_bc, rev_map_bc, output_fp):
    """Imports raw seqsand runs demux on it
    """
    return do_demux_art(input_fp, metadata, metadata_bc_col, rev_bc,
                        rev_map_bc, output_fp)

def do_demux_art(input_fp, metadata, metadata_bc_col, rev_bc, rev_map_bc, output_fp):
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
    start = time.clock()
    if(input_fp is None or metadata is None or metadata_bc_col is None):
        click.echo("Run \'rdemux --help\' flag to see correct usage")
        return

    click.echo("Importing seq data from " + input_fp)
    art = Artifact.import_data("EMPSingleEndSequences", input_fp)

    click.echo("Loading metadata from " + metadata)
    barcode_metadata = Metadata.load(metadata)

    click.echo("Demuxing")
    demux, = emp_single(art,
                        barcode_metadata.get_column(metadata_bc_col),
                        rev_comp_barcodes=rev_bc,
                        rev_comp_mapping_barcodes=rev_map_bc)

    if(output_fp is None):
        click.echo("Not saving demux output")
    else:
        demux.save(output_fp)

    click.echo("{}s for do_demux".format(str(time.clock() - start)))
    return demux, barcode_metadata

@click.command()
@click.option('-i', '--input-fp',
              type=click.Path(exists=True), required = True,
              help="Path to demuxed qza")
@click.option('-l', '--trim-length', type=click.INT, default=100,
              help='Trim length, default set to max possible length.')
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
    start = time.clock()
    click.echo("Importing seq data from " + input_fp)
    input_artifact = Artifact.load(input_fp)

    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    click.echo("{}s for importing for pre".format(str(time.clock() - start)))
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
    start = time.clock()
    if(trim_length is None):
        click.echo("Determining max possible trim length")
        trim_length = methods.get_shortest_seq(input_artifact)

    trim_lengths, percents = calculate_trim_lengths(trim_length, trim_incr, num_trims)

    deblurreds = []
    for l in trim_lengths:
        click.echo("Pre-trimming to length {}".format(str(l)))
        db_out = methods.do_deblur(input_artifact, l, num_cores=num_cores)
        deblurreds.append(db_out)
        if(output_fp is not None):
            db_out.save(output_fp + "/deblurred_pre_" + str(l) + ".qza")

    click.echo("{}s for pre_trims".format(str(time.clock() - start)))
    return deblurreds

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, dir_okay=False),
              help="Path to deblurred qza that we are post-trimming")
@click.option('-o', '--output-fp', type=click.Path(), default=None,
              required=True,
              help='Path to output post-trimmed qza files, and collapse.csv')
@click.option("-tl", "--trim-lengths", type=click.INT, multiple=True, required=False,
              default=[100,90], help="Trim lengths")
@click.option('-on', '--output-name', default="deblurred_pt_", required=False,
              help="Basename for output post trim qza, length appended")
@click.option("-to", "--time-out",
              help="Path to file we are appending time info to")
@click.option("-toa", "--time-out-append",
              help="Identification string to append to time out file")
@click.option("-pc", "--partition-count", type=click.INT,
              help="How many parallel processes to run when"
                   "post trimming")
@click.option("-ib", "--input-biom-fp", type=click.Path(dir_okay=False, exists=True),
              help="Path to biom table we are demuxing")
@click.option("-sb", "--save-biom", is_flag=True)
def post_trims(input_fp, output_fp, trim_lengths,
               output_name, time_out, time_out_append, partition_count,
               input_biom_fp, save_biom):
    start = time.clock()

    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    if(input_fp is None and input_biom_fp is None):
        click.echo("No input given! See --help")
        return
    elif (input_fp is None):
        input_artifact = None
        click.echo("Loading biom table")
        input_biom = biom.load_table(input_biom_fp)
    else:
        click.echo("Importing seq data from " + input_fp)
        input_artifact = Artifact.load(input_fp)
        input_biom=None


    click.echo("{}s for importing for post_trims".format(str(time.clock() - start)))
    click.echo("partition_count: {}".format(partition_count))
    return post_trims_art(output_fp, input_artifact,
                          trim_lengths, output_name, time_out, time_out_append,
                          partition_count, input_biom, save_biom)


def post_trims_art(output_fp, input_artifact=None,
                   trim_lengths=[100,90],
                   output_name="deblurred_pt_", time_out=None,
                   time_out_append=None, partition_count=None, input_biom=None,
                   save_biom=False):
    """Post trims to various specified lengths.
    Saves qza's if specified. With naming format "deblurred_pt_<length>.qza
    eg. If input is length 100, trim_incr=10 and num_trims=5, post trims to
    lengths 100, 90, 80, 70, 60

    Parameters
    ----------
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
    trim_lengths: array_like of INT
        lengths we are trimming to
    output_name: str
        string that postrimmed qzas are named as, with length appended

    Returns
    -------
    list of length trim_lengths of post_trimmed seqs as qiime2.Artifacts of
    type FeatureTable[Frequency]
    pandas.DataFrame of collapse data
    """
    start = time.clock()

    if(input_biom is None):
        input_biom = input_artifact.view(biom.Table)

    #otus = input_biom.ids(axis="observation")
    #trim_length = len(otus[0])
    #for otu in otus:
    #    if(len(otu) != trim_length):
    #        raise ValueError("Input table reads are not all same length. Invalid")
    #        return

    pt_bioms = []
    pt_arts = []
    for l in trim_lengths:
        click.echo("Post-trimming to length {}, pc={}".format(str(l), str(partition_count)))
        pt_biom = methods.post_trim(input_biom, l, partition_count)
        #pt_bioms.append(pt_biom)
        pt_artifact = Artifact.import_data("FeatureTable[Frequency]", pt_biom)
        #pt_arts.append(pt_artifact)
        if(output_fp is not None):
            pt_artifact.save(output_fp + "/" + output_name + str(l) + ".qza")
        if(save_biom):
            with open(output_fp + "/" + output_name + str(l) + ".biom", "w") as file:
                pt_biom.to_json("deblur-testing", file)

    clps = methods.get_collapse_counts(pt_bioms)
    clps.to_csv(output_fp + "/collapse.csv", index=False)

    elapsed = time.clock() - start
    click.echo("{}s for post_trims".format(str(elapsed)))
    if time_out is not None and time_out_append is not None:
        with open(time_out, "a") as file:
            file.write("{}\t{}\n".format(time_out_append, str(elapsed)))
    return pt_arts

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              help="Path to directory holding pre, post qzas, collapse.csv . "
                   "In naming format output by post_trims and pre_trims. "
                   "Required that each post-trim has corresponding pre-trim eg:"
                   "deblurred_pre_100nt.qza")
@click.option("-ip", "--input-path-file", type=click.Path(exists=True, dir_okay=False),
              help="csv file containing paths to .biom or .qza's we are using"
                   "for analysis in format:\n"
                   "<pre_longest>,<post_longest>\n"
                   "<pre_2nd_longest>,<post_2nd_longest>\n")
@click.option("-cfp", "--clps-df", type=click.Path(exists=True, dir_okay=False),
              help="Path to collapse file")
@click.option('-o', '--output-fp',type=click.Path(file_okay=False,exists=True),
              default = None, required=False, help='Path to output csv files')
@click.option('--trim-incr', type=click.INT, default=10, required=False,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5, required=False,
              help='Number of lengths to trim to, default 5')
@click.option("-tl", "--trim-lengths", type=click.INT, multiple=True, required=False,
              help="Trim lengths")
def analysis(input_fp, input_path_file, clps_df, output_fp, trim_incr, num_trims, trim_lengths):

    start = time.clock()

    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    if(input_fp is None and input_path_file is None):
        click.echo("No inputs supplied, see --help!")
        return
    elif(input_fp is not None):
        if input_fp.endswith('/'):
            input_fp = input_fp[:-1]
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

        clps_df = pd.read_csv(input_fp + "/collapse.csv")
    else:
        if clps_df is None:
            click.echo("Supply collapse path!")
            return
        paths = pd.read_csv(input_path_file, header=None)
        pre_artifacts = [load_artifact(x) for x in paths.iloc[:,0]]
        post_artifacts = [load_artifact(x) for x in paths.iloc[:,1]]
        clps_df = pd.read_csv(clps_df)


    click.echo("{}s for loading qza's for analysis"\
               .format(str(time.clock() - start)))

    if len(trim_lengths) == 0:
        trim_lengths = None

    return analysis_art(pre_artifacts, post_artifacts, clps_df,
                        trim_incr, num_trims, output_fp, trim_lengths)

def load_artifact(fp):
    if fp.endswith(".biom"):
        return Artifact.import_data("FeatureTable[Frequency]", biom.load_table(fp))
    else:
        return Artifact.load(fp)

def analysis_art(pre_artifacts, post_artifacts, clps_df=None, trim_incr=10,
                 num_trims=5, output_fp=None, trim_lengths=None,
                 pre_bioms=None, post_bioms=None):
    """Returns analysis data on pre/post artifacts

    Parameters
    ----------
    pre_artifacts: array_like of qiime2.Artifacts type FeatureTable[Frequency]
        pre-trimmed Artifacts in descending trim length order. Should be in
        same order as post_artifacts
    post_artifacts: array_like of qiime2 artifacts
        post-trimmed Artifacts in descending trim length order. Should be in
        same order as pre_artifacts
    clps_df: pd.DataFrame
        Dataframe holding collapse data as output by
        methods.get_collapse_data()
    trim_incr: int, optional
        Percent amount to decrement by. Should correspond to input artifacts
    num_trims: int, optional
        Number of different lengths to trim to. Each trim_incr % less.
        Should correspond to input artifacts
    output_fp: str
        Path to directory which output tsv's are saved, should not end with /
    trim_lengths: array_like of type INT
        list of lenghts we are trimming to
    pre_bioms: array_like of biom.Table
        pre-trimmed bioms in descending trim length order. Should be in
        same order as post_bioms
    post_artifacts: array_like of biom.Table
        post-trimmed bioms in descending trim length order. Should be in
        same order as pre_bioms

    Returns
    -------
    pandas DataFrame of pariwise mantel test data
    pandas DataFrame of distance between pre and post matching otus and
    collapse counts
    pandas DataFrame of distance between pre and post matching samples
    pandas DataFrame of sOTU and sample counts
    pandas DataFrame of changes in reads per sample from pre to post (pre-post)
    """
    start = time.clock()
    print("Enter analysis_art")
    if pre_bioms is None and post_bioms is None:
        pre_bioms = [art.view(biom.Table) for art in pre_artifacts]
        post_bioms = [art.view(biom.Table) for art in post_artifacts]

    if(trim_lengths is None):
        otus = pre_bioms[0].ids(axis="observation")
        length = len(otus[0])
        trim_lengths, prc = calculate_trim_lengths(length, trim_incr,num_trims)

    click.echo("Sorting bioms")

    click.echo("Calculating pre-post distances")
    pre_post, pre_post_sample, pre_overlaps, post_overlaps = \
        methods.get_pre_post_distance_data(pre_bioms, post_bioms, trim_lengths)
    # merge pre_post with collapse data so we can perform correlation
    if clps_df is not None:
        pre_post = pd.merge(pre_post, clps_df[["seq","num_collapses"]], on="seq")

    if(output_fp is not None):
        pre_post.to_csv(output_fp + "/pre_post.csv", index=False)
        pre_post_sample.to_csv(output_fp + "/pre_post_sample.csv", index=False)

    click.echo("Calculating pairwise diversity")
    pairwise_mantel = methods.get_pairwise_diversity_data(pre_bioms, post_bioms,
                                                          trim_lengths)
    if(output_fp is not None):
        pairwise_mantel.to_csv(output_fp + "/pairwise_mantel.csv", index=False)

    click.echo("Calculating count info")
    counts, read_changes = methods.get_count_data(pre_bioms, pre_overlaps,
                                                  post_bioms, post_overlaps,
                                                  trim_lengths)

    if(output_fp is not None):
        counts.to_csv(output_fp + "/counts.csv", index=False)
        read_changes.to_csv(output_fp + "/read_changes.csv", index=False)

    click.echo("{}s for analysis".format(str(time.clock() - start)))
    return pairwise_mantel, pre_post, pre_post_sample, counts, read_changes

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              required=True,
              help="Path to directory holding csv files made by analysis")
@click.option('-o', '--output-fp',type=click.Path(file_okay=False, exists=True),
              default=None, required=False, help='Path to output csv files')
def do_plots(input_fp, output_fp):
    start = time.clock()
    if input_fp.endswith('/'):
        input_fp = input_fp[:-1]
    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    pre_post = pd.read_csv(input_fp + "/pre_post.csv")
    pre_post_sample = pd.read_csv(input_fp + "/pre_post_sample.csv")
    #pairwise_mantel = pd.read_csv(input_fp + "/pairwise_mantel.csv")
    #counts = pd.read_csv(input_fp + "/counts.csv")
    #read_changes = pd.read_csv(input_fp + "/read_changes.csv")
    pairwise_mantel = None
    counts = None
    read_changes = None

    click.echo("{}s for importing data for plots".format(str(time.clock() - start)))

    return plot_pd(pairwise_mantel, pre_post, pre_post_sample, counts, read_changes, output_fp)

def plot_pd(pairwise_mantel, pre_post, pre_post_sample, counts, read_changes, output_fp = None):
    """Plots data generated by analysis()

    Parameters
    ----------
    pairwise_mantel: pandas DataFrame
    pre_post: pandas DataFrame
    pre_post_sample: pandas DataFrame
    counts: pandas DataFrame
    read_changes: pandas DataFrame
        all in formart output by analysis()
    output_fp: str
        directory where to deposit plot png files
    Returns
    -------
    pyplot figure of each respective plot
    """
    start = time.clock()

    pp_plot = sns.boxplot(x="length", y="dist", data=pre_post, hue="dist_type")
    pp_plot.set(xlabel="trim_length", title="Distribution of distances between matching pre/post otus")

    if output_fp is not None:
        plt.savefig(output_fp + "/pre_post.png")
    plt.figure()

    pps_plot = sns.boxplot(x="length", y="dist", data=pre_post_sample, hue="dist_type")
    pps_plot.set(xlabel="trim_length", title="Distribution of distances between matching pre/post samples")

    if output_fp is not None:
        plt.savefig(output_fp + "/pre_post_sample.png")
    plt.figure()
    
    mantel_plot = sns.lmplot(x="trim_length", y="r_sq", row="dist_type", hue="dist_type",
                             data=pairwise_mantel, ci=None, fit_reg=False)
    mantel_plot.set(xlabel="Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/pairwise_mantel.png")
    plt.figure()

    clps_reg = sns.lmplot(x="num_collapses", y="dist", data=pre_post,
                          col="dist_type", ci=None)

    if output_fp is not None:
        plt.savefig(output_fp + "/collapse.png")
    plt.figure()

    plt.plot(counts["trim_length"], counts["sOTU_overlap_count"])
    plt.plot(counts["trim_length"], counts["sOTU_unique_pre"])
    plt.plot(counts["trim_length"], counts["sOTU_unique_post"])
    plt.plot(counts["trim_length"], counts["diff_otu"])
    plt.legend(counts.columns[1:], loc='best')
    plt.title("sOTU counts")
    plt.xlabel("Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/counts.png")
    plt.figure()

    cols = list(read_changes.columns)
    for colname in read_changes.columns[1:]:
        plt.plot(read_changes["trim_length"], read_changes[colname])
    plt.title("pre - post: reads-per-sample count")
    plt.xlabel("Trim Length")

    if output_fp is not None:
        plt.savefig(output_fp + "/read_changes.png")

    click.echo("{}s for plotting".format(str(time.clock() - start)))

@click.command()
@click.option("-i", "--input-fp", required=True,
              type=click.Path(file_okay=False, exists=True),
              help="Path to directory that holds sequences.fastq.gz and"
                   "barcodes.fastq.gz")
@click.option('-m','-metadata', required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Path to sample metadata file")
@click.option('-mbc', '--metadata_bc_col', required=True,
              help="Name of column in metadata file that holds barcodes")
@click.option('-rev_bc', is_flag=True, help="Will reverse barcodes")
@click.option('-rev_map_bc', is_flag=True, help="Will reverse mapping barcodes")
@click.option('-l', '--trim-length', type=click.INT, default=100,
              help='Trim length, default set to max possible length.')
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-nc', '--num-cores', type=click.INT, default=1,
              help="Number of cores to parallelize deblur")
@click.option("-o", "--output-fp", required=True,
              type=click.Path(file_okay=False),
              help="Directory where we will output everything, see other"
                   "functions to see output formats. demux as demux.qza")
def pre_post(input_fp, metadata, metadata_bc_col, rev_bc, rev_map_bc,
             trim_length, trim_incr, num_trims, num_cores, output_fp):
    """Runs the entire pre_post analysis pipeline from demux to plotting,
    defaults to output EVERYTHING
    """
    click.echo("pre_post() starting at " + time.strftime("[%H:%M:%S]"))
    start = time.clock()
    demux, bc_md = do_demux_art(input_fp, metadata, metadata_bc_col, rev_bc,
                            rev_map_bc, output_fp + "/demux.qza")

    pre_arts = pre_trims_art(demux, trim_length, trim_incr, num_trims,
                             output_fp, num_cores)

    pt_arts, clps = post_trims_art(output_fp, pre_arts[0], trim_incr,
                                   num_trims)

    pw_mantel, pre_post, pre_post_sample, counts, read_changes = \
        analysis_art(pre_arts, pt_arts, clps, trim_incr, num_trims, output_fp)

    plot_pd(pw_mantel, pre_post, pre_post_sample, counts, read_changes, output_fp)
    click.echo("{}s for entire pre_post()".format(str(time.clock()-start)))

@click.command()
@click.option("-i", "--input-fp", required=True,
              type=click.Path(dir_okay=False, exists=True),
              help="Paths to pre_trimmed seqs in .biom format", multiple=True)
@click.option("-o", "--output-fp", required=True,
              type=click.Path(file_okay=False),
              help="Directory where we will output everything, see other"
                   " functions to see output formats. demux as demux.qza")
@click.option("-to", "--time-out",
              help="Path to file we are appending time info to")
@click.option("-toa", "--time-out-append",
              help="Identification string to append to time out file")
@click.option("-tl", "--trim-lengths", type=click.INT, multiple=True, required=False,
              help="Trim lengths")
def biom_to_post(input_fp, output_fp, time_out, time_out_append, trim_lengths):
    """Runs the analysis pipeline starting from pre trimmed .biom files
    Do not trim to max length, is implicitly done.
    Eg. if you want to trim ot lengths 150, 100, 90 only specifpy 100 and 90
    """
    start = time.clock()

    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    get_tl = False
    if(len(trim_lengths) == 0):
        trim_lengths = []
        get_tl = True

    pre_arts = []
    pre_bioms = []
    print(input_fp)
    input_fp = list(input_fp)
    click.echo("fooooooo")
    print(input_fp)
    for fp in input_fp:
        as_biom = biom.load_table(fp)
        print(type(as_biom))
        if get_tl:
            trim_lengths.append(get_length_biom(as_biom))
        pre_bioms.append(as_biom)
        print(type(as_biom))
        as_artifact = Artifact.import_data("FeatureTable[Frequency]", as_biom)
        as_artifact.save(output_fp + "/" + os.path.basename(fp))
        pre_arts.append(as_artifact)

    # Sort in descending order by length
    pre_arts = [x for _,x in sorted(zip(trim_lengths, pre_arts), reverse=True)]

    pt_arts, clps = post_trims_art(output_fp, pre_arts[0],
                                   trim_lengths=trim_lengths[1:])

    pt_arts.insert(0, pre_arts[0])

    pw_mantel, pre_post, pre_post_sample, counts, read_changes = \
        analysis_art(pre_arts, pt_arts, clps, trim_lengths=trim_lengths,
                     output_fp=output_fp)

    plot_pd(pw_mantel, pre_post, pre_post_sample, counts, read_changes, output_fp)

    elapsed = time.clock()-start
    click.echo("{}s for entire biom_post()".format(str(elapsed)))
    if time_out is not None and time_out_append is not None:
        with open(time_out, "a") as file:
            file.write("{}\t{}\n".format(time_out_append, str(elapsed)))

@click.command()
@click.option("-i", "--input-fp", type=click.Path(dir_okay=False, exists=True),
              required=True, help="Path to .biom file we want to convert")
@click.option("-o", "--output-fp", type=click.Path(dir_okay=False),
              required=True, help="Path we will output qza to, extension optional")
def biom_to_qiime(input_fp, output_fp):
    """Converts a .biom file to .qza file"""

    as_biom = biom.load_table(input_fp)
    as_artifact = Artifact.import_data("FeatureTable[Frequency]", as_biom)
    as_artifact.save(output_fp)

@click.command()
@click.option("-i", "--input-fp", type=click.Path(dir_okay=False, exists=True),
              required=True, help="Path to .qza file we want to convert")
@click.option("-o", "--output-fp", type=click.Path(dir_okay=False),
              required=True, help="Path we will output .biom to")
def qiime_to_biom(input_fp, output_fp):
    """Converts a .qza file to .biom file"""

    as_artifact = Artifact.load(input_fp)
    as_biom = as_artifact.view(biom.Table)
    as_json = as_biom.to_json(generated_by="deblur-testing")
    with open(output_fp, "w") as f:
        f.write(as_json)

@click.command()
@click.option("-nr", "--n-row", type=click.INT, required=False, default=100,
                help="Number of rows/observations")
@click.option("-nc", "--n-col", type=click.INT, required=False, default=100,
                help="Number of columns/samples")
@click.option("-o", "--output-fp", type=click.Path(dir_okay=False),
                default=None, help="Output path for biom table as json")
@click.option("-l", "--len-seq", type=click.INT, default=150,
                help="Length of sequences")
@click.option("-sc", "--scale", type=click.INT, default=1000,
              help="Factor to scale random numbers by")
def make_rand_biom(n_row, n_col, output_fp, len_seq, scale):
    """
    Makes a random biom table
    """
    oids = []
    for j in range(n_row):
        oid = ""
        for i in range(len_seq):
            oid += choice("ACTG")
        oids.append(oid)

    vals = (np.random.rand(n_row, n_col)*scale).astype(int)
    sids = ["s"+str(x) for x in range(n_col)]
    tbl = biom.Table(vals, oids, sids)

    if output_fp is not None:
        with open(output_fp, "w") as file:
            tbl.to_json("mamaher", file)

    return tbl

@click.command()
@click.option("-mi","--max-length-biom-fp", type=click.Path(dir_okay=False, exists=True),
              required=True, help="Path to max length biom table")
@click.option("-si", "--smaller-biom-fp", type=click.Path(dir_okay=False, exists=True),
              required=True, help="Path to smaller biom table, can be multiple")
@click.option("-o", "--output-fp", type=click.Path(file_okay=False),
              help="Path to output subsetted.biom files to")
@click.option("--identical", is_flag=True, help="Have this if qiita_ids have identical samples accross tables")
def subset_emp_bioms(max_length_biom_fp, smaller_biom_fp, output_fp, identical,
                     max_biom=None, smaller_bioms=None, file_only=True):
    """
    Subsets bioms by qiita ID and filters so they hold exact same samples,
    saves them to folder specified.

    Naming format
    max length intersected with shorter:
        <qiita_id>_<max_length>nt_int_<short_length>nt.biom
    shorter intersect math length:
        <qiita_id>_<short_length>nt_int_<max_length>nt.biom
    """
    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    # Manual biom options are here so testing is easier
    if(max_biom is None and smaller_bioms is None):
        max_biom = biom.load_table(max_length_biom_fp)
        max_length = get_length_biom(max_biom)

        smaller_bioms = []
        lengths = []
        for fp in smaller_biom_fp:
            as_biom = biom.load_table(fp)
            smaller_bioms.append(as_biom)
            lengths.append(get_length_biom(as_biom))

    def partitioner(i, m):
        return i.split('.')[0]  # sample IDs are of the form: <qiita_study_id>.<sample_id>

    # dictionaries mapping qid to table
    max_part_dict = dict([x for x in max_biom.partition(partitioner)])
    smaller_dicts = dict([x for tbl in smaller_bioms for x in tbl.partition(partitioner())])

    for qid in max_part_dict:
        for length, qid_2_table in zip(lengths, smaller_dicts):
            try:
                qs_small = qid_2_table[qid]
            except KeyError:
                print("qid {} not found in table length {}, skipping".format(qid, str(length)))
            qs_max = max_part_dict[qid]

            if(not identical):
                shared_ids = set(qs_max.ids()).intersection(set(qs_small.ids()))
                qs_small.filter(shared_ids, inplace = True)
                qs_max.filter(shared_ids, inplace = True)

            qs_max_name = "{}_{}nt_int_{}nt.biom".format(qid, str(max_length), str(length))
            qs_small_name = "{}_{}nt_int_{}nt.biom".format(qid, str(length), str(max_length))

            if(output_fp is not None):
                with open(output_fp + "/" + qs_max_name) as file:
                    qs_max.to_json(generated_by="deblur_testing_subset_emp_bioms",
                                   direct_io=file)
                with open(output_fp + "/" + qs_small_name) as file:
                    qs_small.to_json(generated_by="deblur_testing_subset_emp_bioms",
                                   direct_io=file)

@click.command()
@click.option("-mi","--main-input-fp", required=True, type=click.Path(dir_okay=False, exists=True),
              help="Path to input biom table")
@click.option("-oi","--ot-input-fp", type=click.Path(dir_okay=False, exists=True),
              multiple=True,
              help="Path to input biom table")
@click.option("-s", "--start", required=True, type=click.INT,
              help="Start # of samples to subsample")
@click.option("-e", "--end", required=True, type=click.INT,
              help="End # of samples to subsample")
@click.option("-c", "--count", required=True, type=click.INT,
              help="# of different subsamples we do in between")
@click.option("-o","--output-fp", type=click.Path(file_okay=False),
              help="Path to output subsampled biom tables to. Naming format:" \
                   "<original name>_ss_<# subsample>.biom")
def subsample_biom(main_input_fp, ot_input_fp, start, end, count, output_fp):
    """Subsamples a biom table by sample many times and subsamples its
    accompanying tables as well with the same IDs

    Does NOT sort
    Does NOT filter observations so resulting bioms might be different size

    Returns
    -------
    output_bioms: array_like
        Array of lists where first element in each element
        is main biom, next n elements are corresponding ot_bioms
        There is an entry for each subsample amount
    """
    if output_fp.endswith('/'):
        output_fp = output_fp[:-1]

    all_fp = list(ot_input_fp)
    all_fp.insert(0, main_input_fp)
    basenames = [os.path.basename(fp).split(".")[0] for fp in all_fp]


    c_start = time.clock()
    main_biom = biom.load_table(main_input_fp)
    # There are 2 samples in 150nt that are not in 100,90nt, filtering them here
    not_shared_ids = ["1064.G.CV298", "2229.W2.N13.EH1.Thomas.CMB.Seaweed.lane5.NoIndex.L005"]
    main_biom.filter(not_shared_ids, invert=True, inplace=True)
    click.echo("{}s for load main biom".format(str(time.clock() - c_start)))

    c_start = time.clock()
    ot_bioms = []
    for fp in ot_input_fp:
        ot_bioms.append(biom.load_table(fp))
    click.echo("{}s for load other bioms".format(str(time.clock() - c_start)))

    output_bioms = []
    for s_count in np.linspace(start, end, count):
        c_start = time.clock()
        list_entry = []
        s_count = int(s_count)

        # Not using biom subsample function since it removes empties
        ss_ids = np.random.choice(main_biom.ids(), size=s_count, replace=False)
        ss_main_biom = main_biom.filter(ss_ids, inplace=False)
        ids = ss_main_biom.ids()

        list_entry.append(ss_main_biom)
        # Subset our other bioms
        for ot_biom in ot_bioms:
            ot_ss = ot_biom.filter(ids, inplace=False)
            list_entry.append(ot_ss)

        for tbl, fn in zip(list_entry, basenames):
            save_path = "{}/{}_ss_{}.biom".format(output_fp, fn, str(s_count))
            with open(save_path, "w") as file:
                tbl.to_json(generated_by="deblur_testing_subsample_biom",
                            direct_io=file)

        output_bioms.append(list_entry)
        click.echo("{}s for s_count {}".format(str(time.clock() - c_start), str(s_count)))

    return output_bioms

@click.command()
@click.option("-pre","--pre-fp", type=click.Path(exists=True, dir_okay=False),
              required=True, multiple=True, help="Path to input pre_trim biom")
@click.option("-pt","--post-fp", type=click.Path(exists=True, dir_okay=False),
              required=True, multiple=True, help="Path to input post_trim biom")
@click.option("-tl", "--trim-length", type=click.INT, multiple=True,
              help="Trim lengths")
@click.option('-o', '--output-fp',type=click.Path(dir_okay=False),
              default = None, required=False,
              help='Path to output pairwise_mantel.csv')
def pairwise_mantel(pre_fp, post_fp, trim_length, output_fp):
    pre_bioms = [biom.load_table(fp) for fp in pre_fp]
    pt_bioms = [biom.load_table(fp) for fp in post_fp]
    if trim_length is None:
        trim_length = np.zeros(len(pre_bioms))

    df = methods.get_pairwise_diversity_data(pre_bioms, pt_bioms, trim_length)
    df.to_csv(output_fp, index=False)

@click.command()
@click.option("-i", "--biom-fp", type=click.Path(dir_okay=False, exists=True),
              required=True, help="Path to biom")
@click.option('-o', '--output-fp', type=click.Path(dir_okay=False),
              default=None, required=True,
              help='Path to output collapse csv')
def get_collapse_count(biom_fp, output_fp):
    """
    Creates csv of collapse counts for a post-trimmed biom table
    """
    clps = methods.get_collapse_counts([biom.load_table(biom_fp)])
    clps.to_csv(output_fp, index=False)


def contains_ids(bioms, ids, axis):
    """Given a list of bioms, determines whether they all have a list of
    ids"""
    for biom in bioms:
        for id in ids:
            if (not biom.exists(id, axis)):
                return False
    return True


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

def get_length_biom(a_biom):
    """Returns length of sequences in a bio table

    Parameters
    ----------
    a_biom: biom.Table
        biom.Table whose seq length we want

    Returns
    -------
    Length of the seqs in input
    """
    otus = a_biom.ids(axis="observation")
    trim_length = len(otus[0])
    for otu in otus:
        if(len(otu) != trim_length):
            raise ValueError("Input table reads are not all same length. Invalid")
            return

    return trim_length
