#!/usr/bin/env python
import argparse
import sys
from argparse import ArgumentParser
import pandas as pd
import tqdm
from clean_asv_data.__main__ import generate_reader, read_clustfile, read_config


def count_clusters(clustdf, countsfile, clust_column, chunksize=None, nrows=None):
    """
    Calculates sums of clusters in each sample

    :param clustdf: Dataframe with ASVs as index and a column with cluster membership
    :param countsfile: Counts of ASVs in each sample
    :param cluster_column: column name of cluster designation
    :param chunksize: Number of rows to read at a time from the countsfile
    :param nrows: Number of total rows to read (development)
    :return: Dataframe with summed counts per cluster
    """
    reader = generate_reader(f=countsfile, chunksize=chunksize, nrows=nrows)
    written = 0
    with sys.stdout as fhout:
        for i, df in enumerate(
            tqdm.tqdm(reader, desc="reading counts", unit=" chunks")
        ):
            if i == 0:
                columns = df.columns
            df_merged = pd.merge(
                df, clustdf, left_index=True, right_index=True, how="inner"
            )
            _cluster_sum = df_merged.groupby(clust_column).sum(numeric_only=True)
            _cluster_sum.loc[:, columns].to_csv(fhout, sep="\t")
            written += _cluster_sum.shape[0]
    return written


def main(args):
    # Read config
    args = read_config(args.configfile, args)
    clustdf = read_clustfile(args.clustfile)
    written = count_clusters(
        clustdf,
        args.countsfile,
        args.clust_column,
        chunksize=args.chunksize,
        nrows=args.nrows,
    )
    sys.stderr.write(f"Wrote {written} clusters\n")


def main_cli():
    parser = ArgumentParser()
    parser.add_argument(
        "--countsfile",
        type=str,
        help="Tab-separated file with counts of ASVs (rows) in samples (columns)",
    )
    parser.add_argument(
        "--clustfile",
        type=str,
        help="Tab-separated file with ASV ids in first column and a column specifying "
        "the cluster it belongs to",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments.",
    )
    parser.add_argument(
        "--clust_column",
        type=str,
        help="Name of cluster column. Defaults to 'cluster'",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        help="If countsfile is very large, specify chunksize to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)
