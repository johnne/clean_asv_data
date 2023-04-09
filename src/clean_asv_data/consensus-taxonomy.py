# Calculate consensus taxonomy of clusters by sorting taxonomic annotations
# by read counts and choosing the assignment at X% of reads

# Input:
# - countsfile
# - clustfile

from argparse import ArgumentParser
import pandas as pd
from clean_asv_data.__main__ import generate_reader, read_clustfile, read_config


def main(args):
    args = read_config(args.configfile, args)
    reader = generate_reader(f=args.countsfile, chunksize=args.chunksize, nrows=args.nrows)

def main_cli():
    parser = ArgumentParser()
    parser.add_argument(
        "--countsfile",
        type=str,
        help="Counts file of ASVs")
    parser.add_argument(
        "--clustfile",
        type=str,
        help="Taxonomy file for ASVs. Should also include a column with cluster designation.")
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments.",
    )
    parser.add_argument(
        "--ranks",
        nargs="+",
        help="Ranks to use for calculating consensus. Must be present in the clustfile"
    )
    parser.add_argument(
        "--consensus_threshold",
        type=int,
        help="Threshold (in %%) at which to assign taxonomy to a cluster")
    parser.add_argument(
        "--chunksize",
        type=int,
        help="If countsfile is very large, specify chunksize to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)