#!/usr/bin/env python
import argparse
import sys
from argparse import ArgumentParser
import pandas as pd
import tqdm
from clean_asv_data.__main__ import (
    generate_reader,
    read_clustfile,
    read_config,
    read_metadata,
)


def sum_clusters(
    clustdf,
    countsfile,
    clust_column,
    blanks=None,
    subset=None,
    chunksize=None,
    nrows=None,
):
    """
    Calculates sums of clusters in each sample

    :param clustdf: Dataframe with ASVs as index and a column with cluster membership
    :param countsfile: Counts of ASVs in each sample
    :param clust_column: column name of cluster designation
    :param chunksize: Number of rows to read at a time from the countsfile
    :param nrows: Number of total rows to read (development)
    :return: Dataframe with summed counts per cluster
    """
    if blanks is None:
        blanks = []
    if subset is None:
        subset = []
    reader = generate_reader(f=countsfile, chunksize=chunksize, nrows=nrows)
    cluster_sum = pd.DataFrame()
    for df in tqdm.tqdm(reader, desc="reading counts", unit=" chunks"):
        if len(subset) > 0:
            subset_intersection = list(set(subset).intersection(set(df.columns)))
            df = df.loc[:, subset_intersection]
        merged = pd.merge(
            clustdf.loc[:, clust_column],
            df.drop(blanks, axis=1, errors="ignore"),
            left_index=True,
            right_index=True,
        )
        _cluster_sum = merged.groupby(clust_column).sum(numeric_only=True)
        cluster_sum = cluster_sum.add(_cluster_sum, fill_value=0)
    return cluster_sum


def main(args):
    # Read config
    args = read_config(args.configfile, args)
    sys.stderr.write(f"Reading {args.clustfile}\n")
    clustdf = read_clustfile(args.clustfile)
    # Read metadata
    metadata = None
    subset = None
    if args.metadata:
        metadata = read_metadata(args.metadata, index_name=args.metadata_index_name)
        # Extract blanks from metadata
        if not args.noblanks:
            blanks = list(
                metadata.loc[metadata[args.sample_type_col].isin(args.blank_val)].index
            )
            sys.stderr.write("####\n" f"Found {len(blanks)} blanks in metadata\n")
        else:
            blanks = []
        if args.subset_val:
            subset = metadata.loc[metadata[args.subset_col] == args.subset_val].index
            sys.stderr.write(
                "####\n"
                f"Found {len(subset)} samples for {args.subset_col}:{args.subset_val}\n"
            )
    cluster_sum = sum_clusters(
        clustdf,
        args.countsfile,
        args.clust_column,
        blanks,
        subset,
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
        "--metadata", type=str, help="Tab-separated file with metadata for each sample"
    )
    parser.add_argument(
        "--metadata_index_name",
        type=str,
        help="Name of column in metadata file that contains sample ids (default: 'sampleID_NGI'))",
        default="sampleID_NGI",
    )
    parser.add_argument(
        "--sample_type_col",
        type=str,
        default="lab_sample_type",
        help="Use this column in metadata to identify sample type (default 'lab_sample_type')",
    )
    parser.add_argument(
        "--blank_val",
        type=str,
        nargs="+",
        default=["buffer_blank", "extraction_neg", "pcr_neg"],
        help="Values in <sample_type_col> that identify blanks (default 'buffer_blank', 'extraction_neg', 'pcr_neg')",
    )
    parser.add_argument(
        "--noblanks",
        action="store_true",
        help="Ignore blanks",
    )
    parser.add_argument(
        "--subset_col",
        type=str,
        default="dataset",
        help="Column in metadata to use for subsetting the counts on (default: 'dataset')",
    )
    parser.add_argument(
        "--subset_val",
        type=str,
        help="Value in subset_col to use for subsetting the counts on",
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
        default="cluster",
        help="Name of cluster column (default: 'cluster')",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        help="If countsfile is very large, specify chunksize to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)
