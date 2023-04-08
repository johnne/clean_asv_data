#!/usr/bin/env python

import pandas as pd
import tqdm
from argparse import ArgumentParser
import sys
from clean_asv_data.__main__ import read_config


def read_counts(f, chunksize=10000, nrows=None):
    """
    Read counts file in chunks and calculate ASV sum and ASV occurrence
    """
    dataframe = pd.DataFrame()
    sys.stderr.write(f"Reading {f} in chunks of {chunksize} lines\n")
    for df in tqdm.tqdm(
        pd.read_csv(f, sep="\t", index_col=0, chunksize=chunksize, nrows=nrows),
        unit=" chunks",
    ):
        asv_occ = pd.DataFrame(df.gt(0).sum(axis=1), columns=["occurrence"])
        asv_sum = pd.DataFrame(df.sum(axis=1), columns=["reads"])
        _dataframe = pd.merge(asv_sum, asv_occ, left_index=True, right_index=True)
        dataframe = pd.concat([dataframe, _dataframe])
    return dataframe


def main(args):
    args = read_config(args.configfile, args)
    dataframe = read_counts(args.countsfile)
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
        help="Path to a yaml-format configuration file. Can be used to set arguments."
    )
    args = parser.parse_args()
    main(args)
