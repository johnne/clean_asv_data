#!/usr/bin/env python

import pandas as pd
import tqdm
from argparse import ArgumentParser
import sys
from clean_asv_data.__main__ import read_config, generate_reader


def read_counts(countsfile, chunksize=None, nrows=None):
    """
    Read counts file in chunks and calculate ASV sum and ASV occurrence
    """
    reader = generate_reader(countsfile, chunksize, nrows)
    dataframe = pd.DataFrame()
    sys.stderr.write(f"Reading {countsfile} in chunks of {chunksize} lines\n")
    for df in tqdm.tqdm(reader, unit=" chunks"):
        asv_occ = pd.DataFrame(df.gt(0).sum(axis=1), columns=["occurrence"])
        asv_sum = pd.DataFrame(df.sum(axis=1), columns=["reads"])
        _dataframe = pd.merge(asv_sum, asv_occ, left_index=True, right_index=True)
        dataframe = pd.concat([dataframe, _dataframe])
    return dataframe


def main(args):
    args = read_config(args.configfile, args)
    dataframe = read_counts(args.countsfile, chunksize=args.chunksize, nrows=args.nrows)
    sys.stderr.write(f"Writing stats for {dataframe.shape[0]} ASVs to stdout\n")
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
