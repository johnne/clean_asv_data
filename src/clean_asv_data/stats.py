#!/usr/bin/env python
import argparse

import pandas as pd
import tqdm
from argparse import ArgumentParser
import sys
from clean_asv_data.__main__ import read_config, generate_reader, read_metadata


def read_counts(countsfile, blanks=None, subset=None, chunksize=None, nrows=None):
    """
    Read counts file in chunks and calculate ASV sum and ASV occurrence
    """
    if blanks is None:
        blanks = []
    if subset is None:
        subset = []
    reader = generate_reader(countsfile, chunksize, nrows)
    dataframe = pd.DataFrame()
    sys.stderr.write(f"Reading {countsfile} in chunks of {chunksize} lines\n")
    for df in tqdm.tqdm(reader, unit=" chunks"):
        if len(subset)>0:
            subset_intersection = list(set(subset).intersection(set(df.columns)))
            df = df.loc[:, subset_intersection]
        # Count number of samples where each ASV occurs
        asv_occ = pd.DataFrame(
            df.drop(blanks, axis=1, errors="ignore").gt(0).sum(axis=1), columns=["occurrence"]
        )
        # Sum counts for each ASV
        asv_sum = pd.DataFrame(df.drop(blanks, axis=1, errors="ignore").sum(axis=1), columns=["reads"])
        _dataframe = pd.merge(asv_sum, asv_occ, left_index=True, right_index=True)
        dataframe = pd.concat([dataframe, _dataframe])
    return dataframe


def main(args):
    args = read_config(args.configfile, args)
    metadata = None
    subset = None
    if args.metadata:
        metadata = read_metadata(args.metadata, index_name=args.metadata_index_name)
        # Extract blanks from metadata
        if not args.noblanks:
            blanks = list(metadata.loc[metadata[args.sample_type_col].isin(args.blank_val)].index)
            sys.stderr.write("####\n" f"Found {len(blanks)} blanks in metadata\n")
        else:
            blanks = []
        if args.subset_val:
            subset = metadata.loc[metadata[args.subset_col] == args.subset_val].index
            sys.stderr.write("####\n" f"Found {len(subset)} samples for {args.subset_col}:{args.subset_val}\n")
    dataframe = read_counts(
        args.countsfile, blanks, subset, chunksize=args.chunksize, nrows=args.nrows
    )
    sys.stderr.write(f"Writing stats for {dataframe.shape[0]} ASVs to stdout\n")
    dataframe.index.name = "ASV"
    with sys.stdout as fhout:
        dataframe.to_csv(fhout, sep="\t")


def main_cli():
    parser = ArgumentParser()
    parser.add_argument(
        "--countsfile",
        type=str,
        help="ASV counts file. Tab-separated, samples in columns, ASVs in rows",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments.",
    )
    parser.add_argument(
        "--metadata", type=str, 
        help="Tab-separated file with metadata for each sample"
    )
    parser.add_argument(
        "--metadata_index_name", type=str, help="Name of column in metadata file that contains sample ids (default: 'sampleID_NGI'))",
        default="sampleID_NGI"
    )
    parser.add_argument(
        "--sample_type_col",
        type=str,
        default="lab_sample_type",
        help="Use this column in metadata to identify sample type (default 'lab_sample_type')",
    ),
    parser.add_argument(
        "--blank_val",
        type=str,
        nargs="+",
        default=['buffer_blank', 'extraction_neg', 'pcr_neg'],
        help="Values in <sample_type_col> that identify blanks (default 'buffer_blank', 'extraction_neg', 'pcr_neg')",
    ),
    parser.add_argument(
        "--noblanks",
        action="store_true",
        help="Ignore blanks",
    ),
    parser.add_argument(
        "--subset_col", type=str, default="dataset",
        help="Column in metadata to use for subsetting the counts on (default: 'dataset')"
    )
    parser.add_argument(
        "--subset_val", type=str,
        help="Value in subset_col to use for subsetting the counts on"
    ),
    parser.add_argument(
        "--chunksize",
        type=int,
        help="Size of chunks (in lines) to read from " "countsfile",
    )
    parser.add_argument(
        "--nrows",
        type=int,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    main(args)
