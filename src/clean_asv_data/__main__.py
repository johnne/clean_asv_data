import pandas as pd


def read_clustfile(f, sep="\t"):
    """
    Reads a cluster membership file for ASVs

    :param f:
    :param sep:
    :return:
    """
    return pd.read_csv(f, sep=sep, index_col=0, header=0)


def generate_reader(f, chunksize, nrows):
    """
    Sets up a reader with pandas. Handles both chunksize>=1 and chunksize=None

    :param f: Input file
    :param chunksize: Number of rows to read per chunk
    :param nrows: Number of total rows to read
    :return:
    """
    r = pd.read_csv(
        f, sep="\t", index_col=0, header=0, nrows=nrows, chunksize=chunksize
    )
    if chunksize is not None:
        return r
    return [r]
