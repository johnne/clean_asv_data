#!/usr/bin/env python
import argparse
import sys
from argparse import ArgumentParser
import pandas as pd
import tqdm
from clean_asv_data.__main__ import generate_reader, read_clustfile, read_config


def sum_clusters(clustdf, countsfile, clust_column, chunksize=None, nrows=None):
    """
    Calculates sums of clusters in each sample

    :param clustdf: Dataframe with ASVs as index and a column with cluster membership
    :param countsfile: Counts of ASVs in each sample
    :param clust_column: column name of cluster designation
    :param chunksize: Number of rows to read at a time from the countsfile
    :param nrows: Number of total rows to read (development)
    :return: Dataframe with summed counts per cluster
    """
    reader = generate_reader(f=countsfile, chunksize=chunksize, nrows=nrows)
    cluster_sum = pd.DataFrame()
    for df in tqdm.tqdm(reader, desc="reading counts", unit=" chunks"):
        merged = pd.merge(clustdf.loc[:, clust_column], df, left_index=True, right_index=True)
        _cluster_sum = merged.groupby(clust_column).sum(numeric_only=True)
        cluster_sum = cluster_sum.add(_cluster_sum, fill_value=0)
    return cluster_sum


def main(args):
    # Read config
    args = read_config(args.configfile, args)
    sys.stderr.write(f"Reading {args.clustfile}\n")
    clustdf = read_clustfile(args.clustfile)
    cluster_sum = sum_clusters(
        clustdf,
        args.countsfile,
        args.clust_column,
        chunksize=args.chunksize,
        nrows=args.nrows,
    )
    with sys.stdout as fhout:
        cluster_sum.to_csv(fhout, sep="\t")



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
