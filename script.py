import click
from qiime2 import Artifact, Metadata
from qiime2.plugins.demux.methods import emp_single
import methods
import biom

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
def pre_trims(input_fp = None, input_artifact = None, trim_length = 100, trim_incr = 10,
              num_trims = 5, output_fp = None, num_cores = 1):
    """Quality filters and then pre_trims sequences to various pre-trim lengths
    Saves qza's if specified. With naming format "deblurred_pre_<length>nt.qza

    Parameters
    ----------
    input_fp: path
        Path to qza of demuxed sequences
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
    if(input_fp is None and input_artifact is None):
        click.echo("No input provided. Provide qiime artifact or path to "
                   "deblurred qza.")
        return
    if(input_artifact is None):
        click.echo("Importing seq data from " + input_fp)
        input_artifact = Artifact.load(input_fp)

    if(trim_length is None):
        click.echo("Determining max possible trim length")
        trim_length = methods.get_shortest_seq(input_artifact)

    # Calculate actual trim lengths
    trim_lengths = []
    percent = 100
    for i in range(num_trims):
        trim_lengths.append(int(trim_length * (percent/100)))
        percent = percent - trim_incr

    deblurreds = []
    for l in trim_lengths:
        db_out = methods.do_deblur(input_artifact, l, num_cores=num_cores)
        deblurreds.append(db_out)
        if(output_fp is not None):
            db_out.save("deblurred_pt" + str(l) + "_nt.qza")


    return deblurreds

@click.command()
@click.option("-i","--input-fp", type=click.Path(exists=True),
              help="Path to deblurred qza that we are post-trimming")
@click.option('--trim-incr', type=click.INT, default=10,
              help='Percent increment amount for different trim lengths, default 10%')
@click.option('-n', '--num-trims', type=click.INT, default=5,
              help='Number of lengths to trim to, default 5')
@click.option('-o', '--output-fp', type=click.Path(),
              help='Path to output post-trimmed qza files, optional')
def post_trims(input_fp = None, input_artifact = None, trim_incr = 10,
               num_trims = 5, output_fp = None):
    """Post trims to various specified lengths.
    Saves qza's if specified. With naming format "deblurred_pt_<length>nt.qza
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
    if(input_fp is None and input_artifact is None):
        click.echo("No input provided. Provide qiime artifact or path to "
                   "deblurred qza.")
        return
    if(input_artifact is None):
        click.echo("Importing seq data from " + input_fp)
        input_artifact = Artifact.load(input_fp)

    input_biom = input_artifact.view(biom.Table)
    otus = input_biom.ids(axis="observation")
    trim_length = len(otus[0])
    for otu in otus:
        print(len(otu))
        if(len(otu) != trim_length):
            click.echo("Input table reads are not all same length. Invalid")
            return

    # Calculate actual trim lengths
    trim_lengths = []
    percent = 100
    for i in range(num_trims):
        trim_lengths.append(int(trim_length * (percent/100)))
        percent = percent - trim_incr

    pt_bioms = []
    for l in trim_lengths:
        pt_biom = methods.post_trim(input_biom, l)
        pt_bioms.append(pt_biom)
        if(output_fp is not None):
            pt_artifact = Artifact.import_data("FeatureTable[Frequency]",
                                               pt_biom)
            pt_artifact.save("deblurred_pt_" + l + "nt.qza")

    return pt_bioms

@click.command()
@click.option()
@click.option("-i","--input-fp", type=click.Path(exists=True, file_okay=False),
              required=True,
              help="Path to directory holding pre and post qza's. In naming " #TODO
              "format output by post_trims and pre_trims. Required that each "
              "post-trim has corresponding pre-trim")
@click.option('-o', '--output-fp', type=click.Path(file_okay=False),
              required=True,
              help='Path to output tsv files')
def analysis(input_fp, input_artifacts, output_fp):
    return
"""
#### Inputs
- all pre-trim and post-trim qzas OR python artifacts
- output dir

#### Outputs
- tsv distances between pre, post otus
- tsv for sOTU counts and pre - post, include percentages col
- tsv for mantel test
- all exported to output_dir/data

#### Pseudocode

```
if(using qzas)
    load all of them to python artifacts

# Distance pre to post plots
# Histograms of BC + Jaccard
# Box plots of distributions by trim length
# Histograms for each trim length
Make single pandas dataframe
export it to csv

# sOTU counts by trim length
# Pre - post reads per sample counts
Do as percentages, put in pandas dataframe
export as tsv


# Pairwise distance Mantel test
# R squared vs. trim length
Export a csv
```
"""

# TODO Shell script for integration testing-- See if working TOGETHER
# MIGHT want to normalize data before bray curtis,
# some biom operations are inplace!
# Norm so each sample is vector of magnitude 1
# GO TO CODE REVIEWS!!!
    # Put on github
    # Pep-8
# TODO test saving
# TODO more testing!!!!
# TODO consistent with 'vs "
# TODO use click groups
# TODO seperate concerns by having a function that acts on artifacts and a function that acts on files that calls the artifact function