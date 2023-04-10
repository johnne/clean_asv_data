#!/usr/bin/env python
import argparse

import pandas as pd
import re
from argparse import ArgumentParser
import tqdm
import sys
from clean_asv_data.__main__ import generate_reader, read_config


def generate_subs(regex, regex_split):
    subs = []
    for r in regex:
        pattern, repl = r.split(regex_split)
        subs.append((pattern, repl))
    return subs


def read_and_rename(f, regex, regex_split, chunksize=None, nrows=None):
    """
    Reads a file and renames the sample names in columns

    :param f: input file
    :param regex: list of regular expressions to use to rename the samples
    :param regex_split: character used to split the regex into pattern and replace
    :param chunksize: number of rows to read from the input file at a time
    :param nrows: number of total rows to read from the file
    :return:
    """
    subs = generate_subs(regex, regex_split)
    reader = generate_reader(f, chunksize, nrows)
    with sys.stdout as fhout:
        for i, df in enumerate(tqdm.tqdm(reader)):
            if i == 0:
                for pattern, repl in subs:
                    df.rename(columns=lambda x: re.sub(pattern, repl, x), inplace=True)
                columns = df.columns
            else:
                df.columns = columns
            df.to_csv(fhout, sep="\t")


def main(args):
    args = read_config(args.configfile, args)
    read_and_rename(
        args.input, args.regex, args.regex_split, args.chunksize, args.nrows
    )


def main_cli():
    parser = ArgumentParser()
    parser.add_argument(
        "input",
        type=str,
        help="Tab-separated input file with sample names in " "columns",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments."
    )
    parser.add_argument(
        "--regex",
        nargs="+",
        help="One or more regular expressions to use to "
        "rename column names of the input file. ",
    )
    parser.add_argument(
        "--regex-split",
        type=str,
        help="Character used to split the regular expressions "
        "into"
        "<pattern> and <repl>. For example with --regex "
        "'FL\d+_L,L the --regex-split ',' will replace "
        "'FL\d+_L' with 'L'. Default ','",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        help="If input file is very large, specify chunksize "
        "to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)


if __name__ == "__main__":
    main_cli()
