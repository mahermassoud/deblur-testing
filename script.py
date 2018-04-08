import click
from qiime2 import Artifact, Metadata
from qiime2.plugins.demux.methods import emp_single
import methods
import biom
import glob
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
@click.option('-o', '--output-fp',  type=click.Path(),
              help="Path to where demuxed qza is saved. Must specify filename. qza extension optional. Does not save if not supplied")
# TEST with rdemux -seqs=test/analysis_testing_wd/mock-3/emp-single-end-sequences/
# -metadata=test/analysis_testing_wd/mock-3/sample-metadata.tsv -out=rdemux_out
# -metadata_bc_col=BarcodeSequence
def demux(input_fp, metadata, metadata_bc_col, rev_bc, rev_map_bc, output_fp = None):
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
    # TODO docstring shows up in help string and it is excessive
    if(input_fp is None or metadata is None or metadata_bc_col is None):
        click.echo("Incorrect. Run \'rdemux --help\' flag to see correct usage")
        #TODO change so it tells you whats wrong
        return

    click.echo("Importing seq data from " + input_fp)
    art = Artifact.import_data("EMPSingleEndSequences", input_fp)

    print("Loading metadata from " + metadata)
    barcode_metadata = Metadata.load(metadata) # TODO exception catching?

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
    # TODO issues testing this one!

@click.command()
@click.option('-i', '--input-fp',
              type=click.Path(exists=True),
              help="Path to demuxed qza")
@click.option('-l', '--trim-length', type=click.INT, default=100,
              help='Trim length, default set to max possible length. '
                   'Supply if possible since calcuation is expensive')
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-o', '--output-fp', type=click.Path(),
              help='Path to output deblurred qza files, optional')
@click.option('-nc', '--num-cores', type=click.INT, default=1,
              help="Number of cores to parallelize deblur")
def pre_trims(input_fp, trim_length = 100, trim_incr = 10,
              num_trims = 5, output_fp = None, num_cores = 1):
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

def pre_trims_art(input_artifact, trim_length = 100, trim_incr = 10,
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

    trim_lengths = calculate_trim_lengths(trim_length, trim_incr, num_trims)

    deblurreds = []
    for l in trim_lengths:
        db_out = methods.do_deblur(input_artifact, l, num_cores=num_cores)
        deblurreds.append(db_out)
        if(output_fp is not None):
            db_out.save("deblurred_pre_" + str(l) + ".qza")

    return deblurreds

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True), required=True,
              help="Path to deblurred qza that we are post-trimming")
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-o', '--output-fp', type=click.Path(),
              help='Path to output post-trimmed qza files, optional')
def post_trims(input_fp, trim_incr = 10,
               num_trims = 5, output_fp = None):
    """Post trims to various specified lengths.
    Saves qza's if specified. With naming format "deblurred_pt_<length>.qza
    eg. If input is length 100, trim_incr=10 and num_trims=5, post trims to
    lengths 100, 90, 80, 70, 60

    Parameters
    ----------
    input_fp: str
        Path to qza of demuxed sequences. Either this or input_artifact
        must be supplied.
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
    click.echo("Importing seq data from " + input_fp)
    input_artifact = Artifact.load(input_fp)

    return post_trims_art(input_artifact, trim_incr, num_trims, output_fp)


def post_trims_art(input_artifact = None, trim_incr = 10,
               num_trims = 5, output_fp = None):
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

    Returns
    -------
    list of length trim_lengths of post_trimmed seqs as biom tables
    """
    input_biom = input_artifact.view(biom.Table)
    otus = input_biom.ids(axis="observation")
    trim_length = len(otus[0])
    for otu in otus:
        print(len(otu))
        if(len(otu) != trim_length):
            raise ValueError("Input table reads are not all same length. Invalid")
            return

    trim_lengths, = calculate_trim_lengths(trim_length, trim_incr, num_trims)

    pt_bioms = []
    for l in trim_lengths:
        pt_biom = methods.post_trim(input_biom, l)
        pt_bioms.append(pt_biom)
        if(output_fp is not None):
            pt_artifact = Artifact.import_data("FeatureTable[Frequency]",
                                               pt_biom)
            pt_artifact.save("deblurred_pt_" + l + ".qza")

    return pt_bioms

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              required=True,
              help="Path to directory holding pre and post qza's. In naming " #TODO
              "format output by post_trims and pre_trims. Required that each "
              "post-trim has corresponding pre-trim")
@click.option('-o', '--output-fp', type=click.Path(file_okay=False),
              required=True,
              help='Path to output tsv files')
def analysis(input_fp, output_fp = None, trim_incr = 10, num_trims = 5):

    pres = dict()
    for path in glob.glob(input_fp + "/deblurred_pre_*.qza"):
        pm = path.replace("_",".")
        pm = pm.split(".")
        length = int(pm[len(pm)-2])
        artifact = Artifact.load(path)
        pres[length] = artifact
    pre_artifacts = [pres[x] for x in sorted(pres.keys())]

    posts = dict()
    for path in glob.glob(input_fp + "/deblurred_pt_*.qza"):
        pm = path.replace("_",".")
        pm = pm.split(".")
        length = int(pm[len(pm)-2])
        artifact = Artifact.load(path)
        posts[length] = artifact
    post_artifacts = [posts[x] for x in sorted(posts.keys())]

    ##Saves qza's if specified. With naming format "deblurred_pre_<length>nt.qza
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
    pre_bioms = [art.view(biom.Table,) for art in pre_artifacts]
    post_bioms = [art.view(biom.Table) for art in post_artifacts]

    otus = pre_bioms.ids(axis="observation")
    length = len(otus[0])
    trim_lengths = calculate_trim_lengths(length, trim_incr,num_trims)

    pairwise_mantel = methods.get_pairwise_diversity_data(pre_bioms, post_bioms,
                                                          trim_lengths)
    pre_post, pre_overlaps, post_overlaps = \
        methods.get_pre_post_distance_data(pre_bioms, post_bioms, trim_lengths)

    counts, read_changes = methods.get_count_data(pre_bioms, pre_overlaps,
                                                  post_bioms, post_overlaps,
                                                  trim_lengths)

    if(output_fp is not None):
        pairwise_mantel.to_csv(output_fp + "/pairwise_mantel.csv")
        pre_post.to_csv(output_fp + "/pre_post.csv")
        counts.to_csv(output_fp + "/counts.csv")
        read_changes.to_csv(output_fp + "/read_changes.csv")

    return pairwise_mantel, pre_post, counts, read_changes

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


# TODO Shell script for integration testing-- See if working TOGETHER
# MIGHT want to normalize data before bray curtis,
# some biom operations are inplace!
# Norm so each sample is vector of magnitude 1
# GO TO CODE REVIEWS!!!
    # Put on github
    # Pep-8
# TODO test saving
# TODO more testing!!!!
# TODO consistent with ' vs "
# TODO use click groups
# TODO ask dan abount seperation of concerns
# TODO do eg.s in python console format
# TODO throw exceptions?
# TODO look @ math behind mantel
# T