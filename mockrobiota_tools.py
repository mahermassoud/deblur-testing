import wget

"""
Utilities for working with mockrobiota dataset
"""
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
